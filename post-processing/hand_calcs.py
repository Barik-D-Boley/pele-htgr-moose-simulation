"""
Hand-calculation VERIFICATION of the MOOSE thermal-stress results.

Idea: at each height, heat made in the fuel has to travel outward through a
chain of "thermal resistances" to reach the helium. Each step has a simple
textbook formula, so we can estimate the temperatures WITHOUT finite elements:

    helium  --(film)-->  coolant wall  --(graphite)-->  fuel wall  --(fuel)-->  fuel center

  1. Helium temperature T_c(z): energy balance along the channel
     (same as coupling/map_power.py).
  2. Film drop:      dT_film  = q'' / h
     q'' = heat flux on the coolant wall = (power per meter) / (wall perimeter)
  3. Graphite drop:  dT_graph = q' * L / (k * W)      (flat-slab approximation)
     The real heat path between the two holes is 2D, so this one is only
     good to roughly +/-30%; it is reported as a low-high range.
  4. Fuel drop:      dT_fuel  = q * R^2 / (4 k)       (heated cylinder, exact)

Axial heat conduction along the column is ignored. That is fine where the
power changes slowly, but a sharp feature (like the OpenMC power spike at the
column ends) gets smeared over roughly one "spreading length" in MOOSE, so
MOOSE should come out somewhat COOLER than this estimate right at the ends.
The script prints that spreading length.

Thermal stresses are given as ORDER-OF-MAGNITUDE scales (sigma ~ E*alpha*dT),
not exact values: MOOSE's peaks should be the same order, not equal.

If MOOSE results exist, the script puts them side by side.

Usage (from the repository root):
    python post-processing/hand_calcs.py
"""
import csv
import glob
import os
import re

import numpy as np

# ------------------------------------------------------------------------------
# Geometry (m) -- same as the CAD / Cubit mesh
# ------------------------------------------------------------------------------
CELL_X, CELL_Y = 0.0188, 0.0094
R_COOL = 0.00794            # coolant hole radius, centered at (0, 0)
R_FUEL = 0.00635            # fuel radius, centered at (CELL_X, CELL_Y)
A_FUEL = 0.25 * np.pi * R_FUEL**2          # quarter fuel cross-section
P_COOL = 0.5 * np.pi * R_COOL               # quarter coolant wall perimeter
P_FUEL = 0.5 * np.pi * R_FUEL               # quarter fuel wall perimeter
GAP = np.hypot(CELL_X, CELL_Y) - R_COOL - R_FUEL   # shortest graphite path

# ------------------------------------------------------------------------------
# Material properties -- must match moose/pele_core_sim.i
# ------------------------------------------------------------------------------
K_FUEL = 15.0                                # W/m-K
E_FUEL, NU_FUEL, ALPHA_FUEL = 15.0e9, 0.15, 5.0e-6
E_GRAPH, NU_GRAPH, ALPHA_GRAPH = 9.8e9, 0.13, 4.0e-6
T_STRESS_FREE = 300.0


def k_graphite(t):
    return 11000.0 / t                       # W/m-K, same fit as MOOSE


# ------------------------------------------------------------------------------
# Read the case settings straight from the MOOSE inputs (so they stay in sync)
# ------------------------------------------------------------------------------
def read_number(path, pattern, default):
    try:
        m = re.search(pattern, open(path).read())
        return float(m.group(1)) if m else default
    except OSError:
        return default


def read_settings(moose_dir):
    base = os.path.join(moose_dir, 'pele_core_sim.i')
    uni = os.path.join(moose_dir, 'case_uniform.i')
    h = read_number(base, r'heat_transfer_coefficient\s*=\s*([0-9.eE+-]+)', 150.0)
    q_avg = read_number(uni, r'value\s*=\s*([0-9.eE+-]+)', 3.0e7)
    m = re.search(r"symbol_values\s*=\s*'([^']+)'", open(uni).read()) if os.path.exists(uni) else None
    t_in, t_out, height = (map(float, m.group(1).split()) if m else (523.15, 973.15, 0.8))
    return dict(h=h, q_avg=q_avg, t_in=t_in, t_out=t_out, height=height)


# ------------------------------------------------------------------------------
# Axial power shape and coolant temperature for each case
# ------------------------------------------------------------------------------
def openmc_axial_shape(path):
    """Average the MOOSE power grid (openmc_power.txt) over the fuel at each z."""
    lines = open(path).read().split('\n')
    xs = np.array(lines[1].split(), float)
    ys = np.array(lines[3].split(), float)
    zs = np.array(lines[5].split(), float)
    data = np.array(' '.join(lines[7:]).split(), float).reshape(len(zs), len(ys), len(xs))
    X, Y = np.meshgrid(xs, ys)                       # shape (ny, nx)
    inside = np.hypot(X - CELL_X, Y - CELL_Y) < R_FUEL
    q = np.array([data[k][inside].mean() for k in range(len(zs))])
    return zs, q


def coolant_profile(z, q, t_in, t_out):
    """Downward flow: helium has picked up the power produced ABOVE height z."""
    # cumulative power from the top down (trapezoid rule)
    seg = 0.5 * (q[1:] + q[:-1]) * np.diff(z)
    above = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
    return t_in + (t_out - t_in) * above / above[0]


# ------------------------------------------------------------------------------
# The resistance chain
# ------------------------------------------------------------------------------
def temperatures(z, q, t_cool, h):
    """Return dict of temperature profiles (K) and drops for one case."""
    q_line = q * A_FUEL                       # W per meter of column (quarter cell)
    dt_film = q_line / P_COOL / h
    t_wall = t_cool + dt_film

    # Graphite slab: low/high path length and width bracket the real 2D path
    w_mid = 0.5 * (P_COOL + P_FUEL)
    dt_graph = {}
    for name, length, width in [('low', GAP, w_mid), ('high', 1.5 * GAP, P_FUEL)]:
        dt = np.zeros_like(z)
        for _ in range(20):                   # k depends on temperature: iterate
            k = k_graphite(t_wall + 0.5 * dt)
            dt = q_line * length / (k * width)
        dt_graph[name] = dt
    dt_g = 0.5 * (dt_graph['low'] + dt_graph['high'])

    dt_fuel = q * R_FUEL**2 / (4.0 * K_FUEL)
    t_fuel_wall = t_wall + dt_g
    t_center = t_fuel_wall + dt_fuel
    return dict(t_cool=t_cool, dt_film=dt_film, t_wall=t_wall, dt_graph=dt_g,
                dt_graph_low=dt_graph['low'], dt_graph_high=dt_graph['high'],
                t_fuel_wall=t_fuel_wall, dt_fuel=dt_fuel, t_center=t_center,
                t_center_low=t_wall + dt_graph['low'] + dt_fuel,
                t_center_high=t_wall + dt_graph['high'] + dt_fuel)


def stress_scales(r, i):
    """Order-of-magnitude thermal stresses (Pa) at height index i."""
    # Radial temperature difference inside the fuel (parabolic profile):
    # surface stress for a free cylinder ~ E*alpha*dT / (2(1-nu))
    s_fuel = E_FUEL * ALPHA_FUEL * r['dt_fuel'][i] / (2 * (1 - NU_FUEL))
    # Temperature difference across the graphite web
    s_graph = E_GRAPH * ALPHA_GRAPH * r['dt_graph'][i] / (2 * (1 - NU_GRAPH))
    # Expansion mismatch between bonded fuel and graphite
    e_eff = 1.0 / (1.0 / E_FUEL + 1.0 / E_GRAPH)
    s_mismatch = e_eff * abs(ALPHA_FUEL - ALPHA_GRAPH) * (r['t_fuel_wall'][i] - T_STRESS_FREE)
    return s_fuel, s_graph, s_mismatch


# ------------------------------------------------------------------------------
# MOOSE results, if available
# ------------------------------------------------------------------------------
def read_moose(moose_dir, case):
    out = {}
    pp = os.path.join(moose_dir, f'pele_{case}_out.csv')
    if os.path.exists(pp):
        with open(pp, newline='') as fh:
            rows = list(csv.DictReader(fh))
        out['pp'] = {k: float(v) for k, v in rows[-1].items()}
    files = sorted(glob.glob(os.path.join(moose_dir, f'pele_{case}_out_fuel_centerline_*.csv')))
    if files:
        with open(files[-1], newline='') as fh:
            rows = list(csv.DictReader(fh))
        out['z'] = np.array([float(r['z']) for r in rows])
        out['T'] = np.array([float(r['temperature']) for r in rows])
    return out


# ------------------------------------------------------------------------------
def main(moose_dir='moose'):
    s = read_settings(moose_dir)
    print(f"Settings from MOOSE inputs: h = {s['h']:g} W/m^2-K, q_avg = {s['q_avg']:.3g} W/m^3, "
          f"helium {s['t_in'] - 273.15:.0f} -> {s['t_out'] - 273.15:.0f} C, "
          f"height {100 * s['height']:.0f} cm\n")

    # Fin-style spreading length: how far heat moves along z before the
    # coolant wall pulls it out. L = sqrt(k A / (h P)), whole cross-section.
    a_solid = CELL_X * CELL_Y - 0.25 * np.pi * R_COOL**2
    spread = np.sqrt(k_graphite(1300.0) * a_solid / (s['h'] * P_COOL))
    print(f'Axial heat spreading length ~{100 * spread:.1f} cm: power features shorter '
          'than this are smoothed out in MOOSE.\n')

    z_u = np.linspace(0.0, s['height'], 161)
    cases = {'uniform': (z_u, np.full_like(z_u, s['q_avg']))}
    power_file = os.path.join(moose_dir, 'openmc_power.txt')
    if os.path.exists(power_file):
        cases['openmc'] = openmc_axial_shape(power_file)
    else:
        print(f'({power_file} not found -- run coupling/map_power.py for the OpenMC case)\n')

    results = {}
    for case, (z, q) in cases.items():
        t_cool = coolant_profile(z, q, s['t_in'], s['t_out'])
        r = temperatures(z, q, t_cool, s['h'])
        results[case] = (z, r)
        i = int(np.argmax(r['t_center']))
        c = lambda t: t - 273.15
        sf, sg, sm = stress_scales(r, i)
        print(f'=== Case: {case} ===')
        print(f"  Hottest fuel at z = {100 * z[i]:.0f} cm: "
              f"{c(r['t_center'][i]):.0f} C  "
              f"(range {c(r['t_center_low'][i]):.0f}-{c(r['t_center_high'][i]):.0f} C)")
        print(f"    helium               {c(r['t_cool'][i]):7.0f} C")
        print(f"    + film (q''/h)       {r['dt_film'][i]:7.0f} K   <- usually the biggest step")
        print(f"    + graphite web       {r['dt_graph'][i]:7.0f} K   "
              f"(range {r['dt_graph_low'][i]:.0f}-{r['dt_graph_high'][i]:.0f})")
        print(f"    + inside the fuel    {r['dt_fuel'][i]:7.0f} K")
        print(f"  Stress scales there:  fuel ~{sf / 1e6:.1f} MPa, graphite web ~{sg / 1e6:.1f} MPa, "
              f"fuel/graphite expansion mismatch ~{sm / 1e6:.1f} MPa")

        m = read_moose(moose_dir, case)
        if 'pp' in m:
            t_m = m['pp']['max_T_fuel']
            est = r['t_center'][i]
            print(f"  MOOSE peak fuel temperature: {c(t_m):.0f} C "
                  f"(hand calc {c(est):.0f} C, difference {t_m - est:+.0f} K)")
            print(f"  MOOSE peak von Mises: fuel {m['pp']['max_vonmises_fuel'] / 1e6:.1f} MPa, "
                  f"graphite {m['pp']['max_vonmises_graphite'] / 1e6:.1f} MPa")
            lo, hi = r['t_center_low'][i] - 30, r['t_center_high'][i] + 30
            if lo <= t_m <= hi:
                print('  -> CONSISTENT with the hand calculation')
            elif t_m < lo and case == 'openmc':
                print('  -> MOOSE is cooler than the hand calc. Expected to some degree at '
                      'the column ends (axial spreading); if it is off by more than ~30% of '
                      'the rise above the helium, check the inputs.')
            else:
                print('  -> CHECK THIS: MOOSE is outside the expected range')
        print()

    plot(results, {c: read_moose(moose_dir, c) for c in results}, 'images/hand_calc_check.png')


def plot(results, moose, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    colors = {'uniform': 'gray', 'openmc': 'C0'}
    fig, ax = plt.subplots(figsize=(6.5, 6))
    for case, (z, r) in results.items():
        col = colors[case]
        ax.fill_betweenx(100 * z, r['t_center_low'] - 273.15, r['t_center_high'] - 273.15,
                         color=col, alpha=0.2)
        ax.plot(r['t_center'] - 273.15, 100 * z, color=col, lw=2,
                label=f'{case}: fuel center (hand calc)')
        ax.plot(r['t_cool'] - 273.15, 100 * z, color=col, lw=1, ls=':',
                label=f'{case}: helium')
        if 'T' in moose[case]:
            ax.plot(moose[case]['T'] - 273.15, 100 * moose[case]['z'], color=col, ls='--',
                    label=f'{case}: fuel center (MOOSE)')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Height in fuel column (cm)  [helium flows down]')
    ax.set_title('Hand-calculation check (shaded = uncertainty range)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower left')
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=150)
    print(f'Saved {out_png}')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--moose-dir', default='moose',
                    help='folder with the MOOSE inputs and results (default: moose)')
    main(ap.parse_args().moose_dir)
