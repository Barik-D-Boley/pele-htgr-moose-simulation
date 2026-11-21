"""
Map the OpenMC fuel power tally onto a grid that MOOSE can read.

What it does, step by step:
  1. Read the 'fuel_power' tally from an OpenMC statepoint file. Each value is
     the fission energy deposited in one small piece (bin) of the fuel compact.
  2. Turn those into power densities (W/m^3), scaled so the AVERAGE over the
     compact equals --q-avg. Both MOOSE cases (uniform and OpenMC) use the
     same --q-avg, so they have exactly the same total power -- only the
     SHAPE differs. That keeps the comparison fair.
  3. Re-sample the cylindrical tally onto a regular x-y-z grid (in meters,
     matching the Cubit mesh) and write it in the text format that MOOSE's
     PiecewiseMultilinear function reads.
  4. Energy-conservation check: integrate the written grid over the fuel
     compact and rescale so total power is preserved exactly.

Usage (from the repository root):
    python coupling/map_power.py
    python coupling/map_power.py --input neutronics/statepoint.150.h5 \
        --output moose/openmc_power.txt --q-avg 3.0e7 --plot
"""
import argparse
import glob
import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Geometry of the fuel compact in METERS (must match the Cubit mesh)
FUEL_CENTER_M = (0.0188, 0.0094)
FUEL_R_M = 0.00635
HEIGHT_M = 0.020

DEFAULT_Q_AVG = 3.0e7  # W/m^3, average power density in the fuel compact


# ------------------------------------------------------------------------------
# Step 1: read the tally
# ------------------------------------------------------------------------------
def read_tally(statepoint_path):
    """Return (r_edges, phi_edges, z_edges, mean, rel_err).

    Edges are in cm (OpenMC units). mean/rel_err have shape (n_r, n_phi, n_z).
    """
    import openmc  # imported here so the rest of the script works without it

    with openmc.StatePoint(statepoint_path) as sp:
        print(f'k-eff = {sp.keff}')
        tally = sp.get_tally(name='fuel_power')
        mesh = tally.find_filter(openmc.MeshFilter).mesh
        dims = tuple(mesh.dimension)
        # OpenMC orders mesh bins with the FIRST index (r) changing fastest,
        # which is Fortran ('F') order.
        mean = tally.mean.ravel().reshape(dims, order='F')
        std = tally.std_dev.ravel().reshape(dims, order='F')
        r_edges = np.asarray(mesh.r_grid)
        phi_edges = np.asarray(mesh.phi_grid)
        z_edges = np.asarray(mesh.z_grid)

    with np.errstate(divide='ignore', invalid='ignore'):
        rel_err = np.where(mean > 0, std / mean, 0.0)
    return r_edges, phi_edges, z_edges, mean, rel_err


# ------------------------------------------------------------------------------
# Step 2: energy per bin -> power density (W/m^3)
# ------------------------------------------------------------------------------
def bin_volumes(r_edges, phi_edges, z_edges):
    """Volume of each cylindrical-mesh bin (any consistent length unit)."""
    ring = 0.5 * np.diff(r_edges**2)        # area factor for each ring
    dphi = np.diff(phi_edges)
    dz = np.diff(z_edges)
    return ring[:, None, None] * dphi[None, :, None] * dz[None, None, :]


def to_power_density(tally_mean, volumes, q_avg):
    """Scale so that the volume-averaged power density equals q_avg."""
    shape = tally_mean / volumes
    avg_shape = tally_mean.sum() / volumes.sum()
    return q_avg * shape / avg_shape


# ------------------------------------------------------------------------------
# Step 3: cylindrical bins -> regular x, y, z grid in meters
# ------------------------------------------------------------------------------
def make_interpolator(r_edges_cm, phi_edges, z_edges_cm, q):
    """Return a function f(x, y, z) [meters] -> Q [W/m^3] using bin centers.

    Points outside the outermost bin centers get the nearest bin's value
    (no wild extrapolation).
    """
    r_c = 0.5 * (r_edges_cm[1:] + r_edges_cm[:-1]) / 100.0   # cm -> m
    p_c = 0.5 * (phi_edges[1:] + phi_edges[:-1])
    z_c = 0.5 * (z_edges_cm[1:] + z_edges_cm[:-1]) / 100.0

    axes = [r_c, p_c, z_c]
    # A single bin along an axis -> duplicate it so interpolation still works
    for i, ax in enumerate(axes):
        if len(ax) == 1:
            axes[i] = np.array([ax[0] - 1e-9, ax[0] + 1e-9])
            q = np.repeat(q, 2, axis=i)
    interp = RegularGridInterpolator(axes, q, method='linear')

    def f(x, y, z):
        dx = x - FUEL_CENTER_M[0]
        dy = y - FUEL_CENTER_M[1]
        r = np.hypot(dx, dy)
        phi = np.mod(np.arctan2(dy, dx), 2 * np.pi)
        pts = np.stack([np.clip(r, axes[0][0], axes[0][-1]),
                        np.clip(phi, axes[1][0], axes[1][-1]),
                        np.clip(z, axes[2][0], axes[2][-1])], axis=-1)
        return interp(pts)

    return f


def build_cartesian_grid(f, nx, ny, nz):
    """Sample f on a box that just covers the quarter fuel compact."""
    xs = np.linspace(FUEL_CENTER_M[0] - FUEL_R_M, FUEL_CENTER_M[0], nx)
    ys = np.linspace(FUEL_CENTER_M[1] - FUEL_R_M, FUEL_CENTER_M[1], ny)
    zs = np.linspace(0.0, HEIGHT_M, nz)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    return xs, ys, zs, f(X, Y, Z)


# ------------------------------------------------------------------------------
# Step 4: energy conservation check
# ------------------------------------------------------------------------------
def integrate_over_fuel(xs, ys, zs, grid, n=200):
    """Integrate the (trilinearly interpolated) grid over the quarter compact.

    This mimics what MOOSE does: it interpolates the grid trilinearly and
    integrates it over the fuel elements. Uses a fine polar midpoint rule.
    Returns watts.
    """
    interp = RegularGridInterpolator((xs, ys, zs), grid, method='linear',
                                     bounds_error=False, fill_value=None)
    dr = FUEL_R_M / n
    dphi = 0.5 * np.pi / n
    nz = max(20, len(zs) * 2)
    dz = HEIGHT_M / nz
    r = (np.arange(n) + 0.5) * dr
    phi = np.pi + (np.arange(n) + 0.5) * dphi
    z = (np.arange(nz) + 0.5) * dz
    R, P, Z = np.meshgrid(r, phi, z, indexing='ij')
    X = FUEL_CENTER_M[0] + R * np.cos(P)
    Y = FUEL_CENTER_M[1] + R * np.sin(P)
    vals = interp(np.stack([X, Y, Z], axis=-1))
    return float(np.sum(vals * R) * dr * dphi * dz)


def quarter_fuel_volume():
    return 0.25 * np.pi * FUEL_R_M**2 * HEIGHT_M


# ------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------
def write_moose_grid(path, xs, ys, zs, grid):
    """Write PiecewiseMultilinear format. X changes fastest in the DATA block."""
    with open(path, 'w') as fh:
        fh.write('AXIS X\n' + ' '.join(f'{v:.8e}' for v in xs) + '\n')
        fh.write('AXIS Y\n' + ' '.join(f'{v:.8e}' for v in ys) + '\n')
        fh.write('AXIS Z\n' + ' '.join(f'{v:.8e}' for v in zs) + '\n')
        fh.write('DATA\n')
        for k in range(len(zs)):
            for j in range(len(ys)):
                fh.write(' '.join(f'{v:.8e}' for v in grid[:, j, k]) + '\n')


def plot_profile(r_edges_cm, q, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    r_mm = 0.5 * (r_edges_cm[1:] + r_edges_cm[:-1]) * 10.0
    radial = q.mean(axis=(1, 2))
    # Area-weighted average = the value the uniform case uses
    uniform = np.average(radial, weights=np.diff(r_edges_cm**2))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(r_mm, radial / 1e6, 'o-', label='OpenMC (averaged over angle and height)')
    ax.axhline(uniform / 1e6, color='gray', ls='--',
               label='Uniform assumption (same total power)')
    ax.set_xlabel('Distance from fuel compact axis (mm)')
    ax.set_ylabel('Power density (MW/m$^3$)')
    ax.set_title('Radial power shape inside the fuel compact')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f'Saved radial profile plot to {out_png}')


# ------------------------------------------------------------------------------
def latest_statepoint(folder='neutronics'):
    files = glob.glob(os.path.join(folder, 'statepoint.*.h5'))
    if not files:
        raise SystemExit(f'No statepoint.*.h5 found in {folder}/ -- run OpenMC first.')
    return max(files, key=os.path.getmtime)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', help='OpenMC statepoint file (default: newest in neutronics/)')
    p.add_argument('--output', default='moose/openmc_power.txt')
    p.add_argument('--q-avg', type=float, default=DEFAULT_Q_AVG,
                   help='average fuel power density in W/m^3 (default %(default)g)')
    p.add_argument('--nx', type=int, default=25)
    p.add_argument('--ny', type=int, default=25)
    p.add_argument('--nz', type=int, default=11)
    p.add_argument('--plot', action='store_true',
                   help='also save images/radial_power_profile.png')
    args = p.parse_args()

    sp_path = args.input or latest_statepoint()
    print(f'Reading {sp_path}')
    r_e, p_e, z_e, mean, rel_err = read_tally(sp_path)
    print(f'Tally bins (r, phi, z): {mean.shape}; '
          f'worst statistical error = {100 * rel_err.max():.2f}%')
    if rel_err.max() > 0.05:
        print('  WARNING: some bins have >5% error. Consider more particles/batches.')

    q = to_power_density(mean, bin_volumes(r_e, p_e, z_e), args.q_avg)
    print(f'Power density in fuel: min {q.min():.3e}, max {q.max():.3e}, '
          f'peak/average = {q.max() / args.q_avg:.3f}')

    f = make_interpolator(r_e, p_e, z_e, q)
    xs, ys, zs, grid = build_cartesian_grid(f, args.nx, args.ny, args.nz)

    target = args.q_avg * quarter_fuel_volume()
    before = integrate_over_fuel(xs, ys, zs, grid)
    grid *= target / before
    after = integrate_over_fuel(xs, ys, zs, grid)
    print(f'Energy check: target {target:.4f} W, after mapping {before:.4f} W '
          f'({100 * (before / target - 1):+.2f}%), rescaled to {after:.4f} W')

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    write_moose_grid(args.output, xs, ys, zs, grid)
    print(f'Wrote {args.output}')

    if args.plot:
        os.makedirs('images', exist_ok=True)
        plot_profile(r_e, q, 'images/radial_power_profile.png')


if __name__ == '__main__':
    main()
