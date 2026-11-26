"""
Effect of OpenMC material temperatures on the power shape and fuel temperature.

Compares, along the fuel centerline:
  - Case A (uniform heat)
  - Case B with room-temperature nuclear data (the baseline OpenMC run)
  - Case B with operating-temperature nuclear data
    (python neutronics/openmc_rve_model.py --operating-temps, then
     coupling/map_power.py and the MOOSE OpenMC case into moose/hot/; see README)

Usage (from the repository root):
    python post-processing/temperature_effect.py
Writes images/temperature_effect.png and prints the peak values.
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = [
    # label, centerline csv, postprocessor csv, color, line style
    ('A: uniform heat', 'moose/pele_uniform_out_fuel_centerline_0001.csv',
     'moose/pele_uniform_out.csv', '#898781', '--'),
    ('B: OpenMC, room-temperature data', 'moose/pele_openmc_out_fuel_centerline_0001.csv',
     'moose/pele_openmc_out.csv', '#2a78d6', '-'),
    ('B: OpenMC, operating-temperature data', 'moose/hot/pele_openmc_hot_out_fuel_centerline_0001.csv',
     'moose/hot/pele_openmc_hot_out.csv', '#eb6834', '-'),
]


def read(path):
    with open(os.path.join(ROOT, path), newline='') as fh:
        rows = list(csv.DictReader(fh))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2), sharey=True)
    for label, line, pp, color, ls in RUNS:
        if not os.path.exists(os.path.join(ROOT, line)):
            print(f'missing {line} -- skipped')
            continue
        d = read(line)
        z = [100 * v for v in d['z']]
        q = d['power_q']
        q_avg = sum(q) / len(q)
        last = {k: v[-1] for k, v in read(pp).items()}
        print(f"{label:<40} axial peaking {max(q) / q_avg:.3f}   peak fuel {last['max_T_fuel'] - 273.15:6.0f} C   "
              f"peak fuel stress {last['max_vonmises_fuel'] / 1e6:5.1f} MPa")
        axes[0].plot([v / q_avg for v in q], z, color=color, ls=ls, lw=2, label=label)
        axes[1].plot([t - 273.15 for t in d['temperature']], z, color=color, ls=ls, lw=2, label=label)
    axes[0].set_xlabel('Power / column average')
    axes[0].set_title('Axial power shape (fuel centerline)', fontsize=10)
    axes[1].set_xlabel('Fuel centerline temperature (°C)')
    axes[1].set_title('Fuel temperature', fontsize=10)
    axes[0].set_ylabel('Height in fuel column (cm)  [helium flows down]')
    for ax in axes:
        ax.grid(color='#e4e3df', lw=0.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    axes[1].legend(loc='upper right', fontsize=8, frameon=False)
    fig.suptitle('Nuclear data temperature: room temperature vs. operating temperature', fontsize=11)
    fig.tight_layout()
    out = os.path.join(ROOT, 'images', 'temperature_effect.png')
    fig.savefig(out, dpi=150)
    print('Saved images/temperature_effect.png')


if __name__ == '__main__':
    main()
