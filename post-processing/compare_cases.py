"""
Compare the two MOOSE cases: uniform heating (Case A) vs. OpenMC heating (Case B).

Reads (from the moose/ folder, after running both cases):
  pele_uniform_out.csv, pele_openmc_out.csv          -- peak values, energy balance, run time
  pele_<case>_out_<line>_*.csv                       -- profiles along the column height
And optionally the OpenMC statepoint, to include the neutronics run time in
the cost comparison.

Writes:
  post-processing/comparison_summary.md   -- table you can paste into the README
  images/case_comparison.png              -- profiles along the column, both cases

Usage (from the repository root):
    python post-processing/compare_cases.py
    python post-processing/compare_cases.py --moose-dir moose --statepoint neutronics/statepoint.250.h5
"""
import argparse
import csv
import glob
import os

import numpy as np

CASES = {'uniform': 'A: Uniform heat', 'openmc': 'B: OpenMC heat'}
LINES = ['fuel_centerline', 'graphite_near_fuel', 'graphite_near_coolant']


# ------------------------------------------------------------------------------
# Reading MOOSE output
# ------------------------------------------------------------------------------
def read_csv_columns(path):
    """Return {column name: numpy array} for a MOOSE CSV file."""
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f'{path} is empty')
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def read_postprocessors(moose_dir, case):
    path = os.path.join(moose_dir, f'pele_{case}_out.csv')
    if not os.path.exists(path):
        raise SystemExit(f'Missing {path}. Run the {case} case first (see README).')
    cols = read_csv_columns(path)
    # Steady runs write the initial state and the solution; keep the last row
    return {k: v[-1] for k, v in cols.items()}


def read_line(moose_dir, case, line):
    """Profile along one sampling line (latest output file), or None."""
    files = sorted(glob.glob(os.path.join(moose_dir, f'pele_{case}_out_{line}_*.csv')))
    return read_csv_columns(files[-1]) if files else None


def openmc_runtime(statepoint):
    """Total OpenMC wall time in seconds from a statepoint, or None."""
    if not statepoint or not os.path.exists(statepoint):
        return None
    import h5py
    with h5py.File(statepoint, 'r') as f:
        return float(f['runtime/total'][()])


def latest_statepoint(folder='neutronics'):
    files = glob.glob(os.path.join(folder, 'statepoint.*.h5'))
    return max(files, key=os.path.getmtime) if files else None


# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------
def build_rows(pp, t_openmc):
    """List of (label, uniform value, openmc value, unit, format)."""
    u, o = pp['uniform'], pp['openmc']

    def c(v):  # K -> C
        return v - 273.15

    def mpa(v):
        return v / 1e6

    def balance(d):
        return 100.0 * (d['power_out'] / d['power_in'] - 1.0)

    rows = [
        ('Peak fuel temperature', c(u['max_T_fuel']), c(o['max_T_fuel']), '°C', '{:.0f}'),
        ('Height of peak fuel temperature', 100 * u['z_of_max_T_fuel'],
         100 * o['z_of_max_T_fuel'], 'cm', '{:.0f}'),
        ('Peak graphite temperature', c(u['max_T_graphite']), c(o['max_T_graphite']), '°C', '{:.0f}'),
        ('Peak von Mises stress, fuel', mpa(u['max_vonmises_fuel']),
         mpa(o['max_vonmises_fuel']), 'MPa', '{:.1f}'),
        ('Peak von Mises stress, graphite', mpa(u['max_vonmises_graphite']),
         mpa(o['max_vonmises_graphite']), 'MPa', '{:.1f}'),
        ('Total heat generated', u['power_in'], o['power_in'], 'W', '{:.1f}'),
        ('Energy balance error (out vs. in)', balance(u), balance(o), '%', '{:+.2f}'),
        ('MOOSE run time', u['run_time'], o['run_time'], 's', '{:.1f}'),
    ]
    # Diagnostic: heat out from -k*grad(T) at the wall (reads low on a coarse
    # mesh). Only in CSVs from runs that have the power_out_gradient postprocessor.
    if 'power_out_gradient' in u and 'power_out_gradient' in o:
        def grad_balance(d):
            return 100.0 * (d['power_out_gradient'] / d['power_in'] - 1.0)
        rows.insert(7, ('Wall-gradient heat flux vs. in (mesh diagnostic)',
                        grad_balance(u), grad_balance(o), '%', '{:+.2f}'))
    if t_openmc is not None:
        rows.append(('OpenMC run time (needed for Case B only)', 0.0, t_openmc, 's', '{:.0f}'))
        rows.append(('Total compute time', u['run_time'], o['run_time'] + t_openmc, 's', '{:.0f}'))
    return rows


def format_table(rows):
    lines = ['| Quantity | Case A: Uniform | Case B: OpenMC | Difference (B − A) |',
             '|---|---|---|---|']
    for label, a, b, unit, fmt in rows:
        diff = b - a
        # Percent change only makes sense for quantities with a true zero
        # (not for Celsius temperatures, heights or percentages)
        rel = f' ({100 * diff / a:+.0f}%)' if a != 0 and unit in ('MPa', 'W', 's') else ''
        lines.append(f'| {label} | {fmt.format(a)} {unit} | {fmt.format(b)} {unit} | '
                     f'{fmt.format(diff)} {unit}{rel} |')
    return '\n'.join(lines)


# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
def plot_profiles(profiles, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    panels = [
        ('fuel_centerline', 'power_q', 1e-6, 'Fission power density (MW/m$^3$)', 'Heat source (fuel)'),
        ('fuel_centerline', 'temperature', None, 'Temperature (°C)', 'Fuel centerline temperature'),
        ('graphite_near_fuel', 'vonmises_stress', 1e-6, 'von Mises stress (MPa)',
         'Graphite stress next to fuel'),
        ('graphite_near_coolant', 'vonmises_stress', 1e-6, 'von Mises stress (MPa)',
         'Graphite stress next to coolant'),
    ]
    styles = {'uniform': dict(color='gray', ls='--'), 'openmc': dict(color='C0', ls='-')}

    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 5), sharey=True)
    for ax, (line, var, scale, xlabel, title) in zip(axes, panels):
        for case, label in CASES.items():
            data = profiles[case].get(line)
            if data is None or var not in data:
                continue
            x = data[var] - 273.15 if var == 'temperature' else data[var] * scale
            ax.plot(x, 100 * data['z'], label=label, **styles[case])
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel('Height in fuel column (cm)  [helium flows down]')
    axes[0].legend(loc='best', fontsize=8)
    fig.suptitle('Uniform vs. OpenMC heat source: profiles along the fuel column')
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f'Saved {out_png}')


# ------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--moose-dir', default='moose')
    p.add_argument('--statepoint', help='OpenMC statepoint (default: newest in neutronics/)')
    p.add_argument('--summary', default='post-processing/comparison_summary.md')
    p.add_argument('--plot', default='images/case_comparison.png')
    args = p.parse_args()

    pp = {case: read_postprocessors(args.moose_dir, case) for case in CASES}
    profiles = {case: {line: read_line(args.moose_dir, case, line) for line in LINES}
                for case in CASES}
    t_openmc = openmc_runtime(args.statepoint or latest_statepoint())

    table = format_table(build_rows(pp, t_openmc))
    print(table)

    for case, d in pp.items():
        err = abs(d['power_out'] / d['power_in'] - 1.0)
        if err > 0.01:
            print(f'WARNING: {case} case energy balance is off by {100 * err:.1f}% '
                  '-- check the mesh, BCs or solver convergence before trusting results.')

    os.makedirs(os.path.dirname(args.summary) or '.', exist_ok=True)
    with open(args.summary, 'w') as fh:
        fh.write('# Case comparison: uniform vs. OpenMC heat source\n\n')
        fh.write('Generated by `post-processing/compare_cases.py`.\n\n')
        fh.write(table + '\n\n')
        fh.write('![Profiles along the fuel column](../images/case_comparison.png)\n')
    print(f'Saved {args.summary}')

    if any(v is not None for d in profiles.values() for v in d.values()):
        os.makedirs(os.path.dirname(args.plot) or '.', exist_ok=True)
        plot_profiles(profiles, args.plot)
    else:
        print('No line profiles found -- skipping plot.')


if __name__ == '__main__':
    main()
