"""
Check that the OpenMC fission source settled before tallying started.

Plots Shannon entropy and k-eff per batch, with a line where the inactive
batches end. The entropy should be flat (just noisy) to the right of that
line. If it is still drifting there, increase settings.inactive.

Usage (from the repository root):
    python neutronics/check_convergence.py            # newest statepoint
    python neutronics/check_convergence.py neutronics/statepoint.250.h5
"""
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openmc

path = sys.argv[1] if len(sys.argv) > 1 else max(
    glob.glob('neutronics/statepoint.*.h5'), key=os.path.getmtime)

with openmc.StatePoint(path) as sp:
    entropy = np.asarray(sp.entropy)
    k_gen = np.asarray(sp.k_generation)
    n_inactive = sp.n_inactive
    keff = sp.keff

batches = np.arange(1, len(entropy) + 1)
active = entropy[n_inactive:]
half = len(active) // 2
drift = abs(active[half:].mean() - active[:half].mean())
noise = active.std()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
ax1.plot(batches, entropy, lw=1)
ax1.axvline(n_inactive, color='r', ls='--', label='start of active batches')
ax1.set_ylabel('Shannon entropy')
ax1.legend()
ax1.grid(alpha=0.3)
ax2.plot(batches, k_gen, lw=1)
ax2.axvline(n_inactive, color='r', ls='--')
ax2.axhline(keff.nominal_value, color='k', lw=0.8, label=f'k-eff = {keff}')
ax2.set_xlabel('Batch')
ax2.set_ylabel('k (per batch)')
ax2.legend()
ax2.grid(alpha=0.3)
fig.tight_layout()
os.makedirs('images', exist_ok=True)
fig.savefig('images/source_convergence.png', dpi=150)

print(f'k-eff = {keff}')
print(f'Entropy drift between 1st and 2nd half of active batches: {drift:.4f} '
      f'(batch-to-batch noise {noise:.4f})')
print('Looks converged.' if drift < noise else
      'WARNING: entropy still drifting -- increase settings.inactive.')
print('Saved images/source_convergence.png')
