"""
Coolant sensitivity sweep: how hot does the fuel get as the helium flow changes?

The baseline (250 -> 700 C helium) needs only ~1.3 g/s per coolant channel to
carry the fuel power, so the flow is barely turbulent and h is low (~150
W/m^2-K). Most of the fuel temperature rise is the helium film. This sweep
re-runs the OpenMC-power MOOSE case for:

  1. FLOW sweep -- more helium flow per channel at the same 250 C inlet. More
     flow means a smaller temperature rise (lower outlet) AND a higher h (from
     the same correlation map_power.py uses). This is the physically
     consistent sweep.
  2. h-only sweep -- the baseline helium temperatures, but a larger h, e.g. if
     the correlation underestimates h or the channel surface is enhanced. This
     separates the film effect from the outlet temperature.

All runs use the same power (OpenMC shape, same q_avg) and the same mesh.

Usage (all from the repository root, in WSL):
    # 1. Write the coolant profiles and the run list
    python post-processing/coolant_sweep.py setup
    # 2. Run MOOSE (moose conda env; runs 3 cases at a time, 4 ranks each)
    cd moose && bash sweep/run_sweep.sh && cd ..
    # 3. Collect the results, plot, write the summary
    python post-processing/coolant_sweep.py report
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'coupling'))
from map_power import HELIUM_CP, estimate_htc  # noqa: E402  (same correlation)

MOOSE_DIR = os.path.join(ROOT, 'moose')
SWEEP_DIR = os.path.join(MOOSE_DIR, 'sweep')
BASE_COOLANT = os.path.join(MOOSE_DIR, 'openmc_coolant_temperature.csv')
BASE_CSV = os.path.join(MOOSE_DIR, 'pele_openmc_out.csv')

# Flow multipliers relative to the baseline flow (1 = 250 -> 700 C)
FLOW_FACTORS = [1.0, 1.5, 2.25, 3.0, 4.5, 6.0]
# h-only sweep, W/m^2-K, at the baseline helium temperatures
H_VALUES = [300.0, 600.0, 1200.0, 2400.0]

# TRISO reference temperatures (commonly quoted, for orientation only)
T_NORMAL_LIMIT_C = 1250.0
T_ACCIDENT_LIMIT_C = 1600.0


# ------------------------------------------------------------------------------
# Inputs: baseline settings straight from the files the baseline run used
# ------------------------------------------------------------------------------
def read_h_wall():
    text = open(os.path.join(MOOSE_DIR, 'pele_core_sim.i')).read()
    m = re.search(r'(?m)^\s*h_wall\s*=\s*([0-9.eE+-]+)', text)
    if not m:
        raise SystemExit('h_wall not found at the top of moose/pele_core_sim.i')
    return float(m.group(1))


def read_coolant(path):
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    return (np.array([float(r['z']) for r in rows]),
            np.array([float(r['T_coolant']) for r in rows]))


def read_pp(path):
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    return {k: float(v) for k, v in rows[-1].items()}


def baseline():
    if not os.path.exists(BASE_COOLANT) or not os.path.exists(BASE_CSV):
        raise SystemExit('Run the baseline OpenMC case first (see README, steps 2 and 4).')
    z, t = read_coolant(BASE_COOLANT)
    t_in, t_out = t.min(), t.max()          # coldest helium = inlet, hottest = outlet
    if not (z[np.argmin(t)] > z[np.argmax(t)]):  # flow is down: inlet at the top
        raise SystemExit('Expected downward flow (hottest helium at the bottom).')
    power = read_pp(BASE_CSV)['power_in']   # W, quarter cell
    m_dot = power / (HELIUM_CP * (t_out - t_in))
    return dict(z=z, t=t, t_in=t_in, t_out=t_out, power=power, m_dot=m_dot,
                h=read_h_wall())


def runs(b):
    """List of sweep cases: name, kind, flow factor, T_out, h, coolant profile."""
    frac = (b['t'] - b['t_in']) / (b['t_out'] - b['t_in'])   # heat picked up so far, 0..1
    out = []
    for f in FLOW_FACTORS:
        dt = (b['t_out'] - b['t_in']) / f
        t_out = b['t_in'] + dt
        _, h = estimate_htc(b['m_dot'] * f, 0.5 * (b['t_in'] + t_out))
        out.append(dict(name=f'flow_x{f:g}', kind='flow', factor=f, t_out=t_out, h=h,
                        t=b['t_in'] + dt * frac))
    for h in H_VALUES:
        out.append(dict(name=f'h_{h:g}', kind='h', factor=1.0, t_out=b['t_out'], h=h,
                        t=b['t']))
    return out


# ------------------------------------------------------------------------------
# setup: coolant CSVs + run script
# ------------------------------------------------------------------------------
def setup(args):
    b = baseline()
    os.makedirs(SWEEP_DIR, exist_ok=True)
    cmds = []
    print(f"Baseline: {b['t_in'] - 273.15:.0f} -> {b['t_out'] - 273.15:.0f} C, "
          f"{4e3 * b['m_dot']:.2f} g/s per channel, h = {b['h']:.0f} W/m^2-K "
          f"(correlation gives {estimate_htc(b['m_dot'], 0.5 * (b['t_in'] + b['t_out']))[1]:.0f})")
    for r in runs(b):
        cool = os.path.join(SWEEP_DIR, f"coolant_{r['name']}.csv")
        with open(cool, 'w') as fh:
            fh.write('z,T_coolant\n')
            for zi, ti in zip(b['z'], r['t']):
                fh.write(f'{zi:.6e},{ti:.6f}\n')
        cmds.append(
            f"run {r['name']} h_wall={r['h']:.3f} "
            f"Functions/coolant_temperature/data_file=sweep/coolant_{r['name']}.csv "
            f"Outputs/file_base=sweep/{r['name']} Outputs/exodus=false")
        print(f"  {r['name']:<12} T_out {r['t_out'] - 273.15:6.0f} C   h {r['h']:6.0f} W/m^2-K")

    script = os.path.join(SWEEP_DIR, 'run_sweep.sh')
    with open(script, 'w', newline='\n') as fh:
        fh.write('#!/bin/bash\n'
                 '# Written by post-processing/coolant_sweep.py setup. Run from moose/.\n'
                 f'# {args.jobs} cases at a time, {args.ranks} MPI ranks each.\n'
                 'EXE=${COMBINED:-combined-opt}\n'
                 'run() {\n'
                 '  name=$1; shift\n'
                 f'  mpiexec -n {args.ranks} "$EXE" -i pele_core_sim.i case_openmc.i "$@" '
                 '> sweep/$name.log 2>&1 && echo "done   $name" || echo "FAILED $name (see sweep/$name.log)"\n'
                 '}\n')
        for i, c in enumerate(cmds):
            fh.write(c + ' &\n')
            if (i + 1) % args.jobs == 0:
                fh.write('wait\n')
        fh.write('wait\n')
    print(f'Wrote {len(cmds)} coolant profiles and {os.path.relpath(script, ROOT)}')


# ------------------------------------------------------------------------------
# report: table + plot
# ------------------------------------------------------------------------------
def collect():
    b = baseline()
    res = []
    for r in runs(b):
        path = os.path.join(SWEEP_DIR, f"{r['name']}.csv")
        if not os.path.exists(path):
            print(f"  missing {os.path.relpath(path, ROOT)} -- skipped")
            continue
        pp = read_pp(path)
        r.update(T_fuel=pp['max_T_fuel'] - 273.15, T_graph=pp['max_T_graphite'] - 273.15,
                 s_fuel=pp['max_vonmises_fuel'] / 1e6, s_graph=pp['max_vonmises_graphite'] / 1e6,
                 z_peak=100 * pp['z_of_max_T_fuel'],
                 balance=100 * (pp['power_out'] / pp['power_in'] - 1),
                 flow=4e3 * b['m_dot'] * r['factor'])
        res.append(r)
    return b, res


def table(b, res):
    lines = ['| Run | He flow per channel | He outlet | h | Peak fuel T | Peak graphite T | '
             'Peak von Mises fuel / graphite | Energy balance |',
             '|---|---|---|---|---|---|---|---|']
    for r in res:
        label = (f"Flow ×{r['factor']:g}" if r['kind'] == 'flow' else f"h only = {r['h']:.0f}")
        lines.append(
            f"| {label} | {r['flow']:.2f} g/s | {r['t_out'] - 273.15:.0f} °C | {r['h']:.0f} W/m²·K | "
            f"{r['T_fuel']:.0f} °C | {r['T_graph']:.0f} °C | "
            f"{r['s_fuel']:.1f} / {r['s_graph']:.1f} MPa | {r['balance']:+.2f} % |")
    return '\n'.join(lines)


def crossing(x, y, level):
    """x where y first drops below level (linear interpolation), or None."""
    for i in range(len(y) - 1):
        if y[i] >= level > y[i + 1]:
            return x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
    return None


def plot(b, res, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    series = '#2a78d6'      # categorical slot 1
    ink2 = '#52514e'        # secondary text / reference lines
    grid = '#e4e3df'
    flow = [r for r in res if r['kind'] == 'flow']
    honly = sorted([r for r in res if r['kind'] == 'h' or (r['kind'] == 'flow' and r['factor'] == 1.0)],
                   key=lambda r: r['h'])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    panels = [
        (axes[0], [r['flow'] for r in flow], [r['T_fuel'] for r in flow],
         'Helium flow per coolant channel (g/s)',
         'More helium flow (same 250 °C inlet)', flow),
        (axes[1], [r['h'] for r in honly], [r['T_fuel'] for r in honly],
         'Heat transfer coefficient h (W/m²·K)',
         'Higher h only (baseline helium temperatures)', honly),
    ]
    for ax, x, y, xlabel, title, rows in panels:
        for level, text in ((T_ACCIDENT_LIMIT_C, '~1600 °C TRISO accident limit'),
                            (T_NORMAL_LIMIT_C, '~1250 °C typical normal-operation limit')):
            ax.axhline(level, color=ink2, lw=1, ls='--', zorder=1)
            ax.text(x[-1], level + 12, text, color=ink2, fontsize=8, ha='right', va='bottom')
        ax.plot(x, y, color=series, lw=2, marker='o', ms=7, zorder=3,
                markeredgecolor='white', markeredgewidth=1.5)
        for r, xi, yi in zip(rows, x, y):
            if ax is axes[0]:
                note, off = f"{r['t_out'] - 273.15:.0f} °C out\nh {r['h']:.0f}", (7, 4)
            else:
                note, off = f"{yi:.0f} °C", (7, -12)
            ax.annotate(note, (xi, yi), textcoords='offset points', xytext=off,
                        fontsize=7.5, color=ink2)
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=10)
        ax.grid(color=grid, lw=0.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    from matplotlib.ticker import NullFormatter, NullLocator
    axes[1].set_xscale('log')
    axes[1].xaxis.set_minor_locator(NullLocator())
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    axes[1].set_xticks([r['h'] for r in honly])
    axes[1].set_xticklabels([f"{r['h']:.0f}" for r in honly])
    axes[0].set_ylim(top=T_ACCIDENT_LIMIT_C + 120)
    axes[0].set_ylabel('Peak fuel temperature (°C)')
    fig.suptitle('Coolant sensitivity: peak fuel temperature, OpenMC power shape', fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f'Saved {os.path.relpath(out_png, ROOT)}')


def report(args):
    b, res = collect()
    if not res:
        raise SystemExit('No sweep results found -- run moose/sweep/run_sweep.sh first.')
    t = table(b, res)
    print(t)
    flow = [r for r in res if r['kind'] == 'flow']
    honly = sorted([r for r in res if r['kind'] == 'h' or r['factor'] == 1.0], key=lambda r: r['h'])
    notes = []
    for level, name in ((T_ACCIDENT_LIMIT_C, '1600 °C'), (T_NORMAL_LIMIT_C, '1250 °C')):
        if flow[0]['T_fuel'] < level:
            notes.append(f"- Peak fuel is already below {name} at the baseline flow "
                         f"({flow[0]['T_fuel']:.0f} °C), with only {level - flow[0]['T_fuel']:.0f} K margin.")
            continue
        xf = crossing([r['flow'] for r in flow], [r['T_fuel'] for r in flow], level)
        xh = crossing([r['h'] for r in honly], [r['T_fuel'] for r in honly], level)
        notes.append(f"- Peak fuel drops below {name} at about "
                     + (f"{xf:.1f} g/s per channel ({xf / flow[0]['flow']:.1f}× baseline flow)"
                        if xf else 'a flow outside the sweep')
                     + ' in the flow sweep, and at about '
                     + (f"h ≈ {xh:.0f} W/m²·K" if xh else 'an h outside the sweep')
                     + ' with the baseline helium temperatures.')
    print('\n'.join(notes))
    for r in res:
        if abs(r['balance']) > 0.5:
            print(f"WARNING: {r['name']} energy balance {r['balance']:+.2f}%")

    os.makedirs(os.path.dirname(args.summary), exist_ok=True)
    with open(args.summary, 'w', encoding='utf-8') as fh:
        fh.write('# Coolant sensitivity sweep (OpenMC power shape)\n\n'
                 'Generated by `post-processing/coolant_sweep.py`. Same power and mesh in every run; '
                 'helium inlet fixed at 250 °C, flowing down. Flow-sweep h comes from the same '
                 'correlation as `coupling/map_power.py`.\n\n')
        fh.write(t + '\n\n' + '\n'.join(notes) + '\n\n')
        fh.write('TRISO limits are commonly quoted reference values, for orientation only.\n\n')
        fh.write('![Coolant sensitivity](../images/coolant_sweep.png)\n')
    print(f'Saved {os.path.relpath(args.summary, ROOT)}')
    plot(b, res, args.plot)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('setup', help='write coolant profiles and moose/sweep/run_sweep.sh')
    s.add_argument('--jobs', type=int, default=3, help='cases run at the same time')
    s.add_argument('--ranks', type=int, default=4, help='MPI ranks per case')
    r = sub.add_parser('report', help='collect results, plot, write summary')
    r.add_argument('--summary', default=os.path.join(ROOT, 'post-processing', 'coolant_sweep_summary.md'))
    r.add_argument('--plot', default=os.path.join(ROOT, 'images', 'coolant_sweep.png'))
    args = p.parse_args()
    setup(args) if args.cmd == 'setup' else report(args)


if __name__ == '__main__':
    main()
