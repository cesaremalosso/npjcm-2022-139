# Ab initio stress time series of liquid water

Stress-tensor time series from the Car-Parrinello molecular dynamics (CPMD) simulations of liquid water used in:

> C. Malosso, L. Zhang, R. Car, S. Baroni, D. Tisi,
> *Viscosity in water from first-principles and deep-neural-network simulations*,
> npj Computational Materials **8**, 139 (2022).
> https://doi.org/10.1038/s41524-022-00830-7

These time series were used to compute the PBE shear viscosity through the Green-Kubo relation, with cepstral analysis of the stress power spectrum ([SporTran](https://github.com/sportran/sportran)). The same data are also archived on Materials Cloud: https://doi.org/10.24435/materialscloud:x7-b0

## Simulation details

| | |
|---|---|
| Code | Quantum ESPRESSO, `cp.x` (Car-Parrinello) |
| Functional | PBE |
| Pseudopotentials | Hamann-Schlüter-Chiang-Vanderbilt norm-conserving |
| Plane-wave cutoff | 85 Ry |
| System | 64 H₂O molecules, cubic cell, edge 12.43 Å (density 1 g/cm³) |
| Time step | 3 a.u. = 0.0725665 fs |
| Ensemble | NVT equilibration (Nosé-Hoover chains), then NVE production (these data) |
| Temperatures | 454 K and 600 K |
| Length | 400 ps per temperature |

## Files

```
ab_initio-time_series/
├── stress_454K.npz      # 454 K, 2,756,093 frames
├── stress_600K.npz      # 600 K, 2,756,174 frames
└── convert_to_npz.py    # script used to create the .npz from the original cp.x .str output
```

Each `.npz` archive contains:

| key | shape / type | description |
|---|---|---|
| `stress` | `(N, 6)` float32 | stress tensor in **GPa**, columns `xx, yy, zz, xy, xz, yz` |
| `components` | `(6,)` str | names of the `stress` columns |
| `step` | `(N,)` int32 | CP step number (one frame every 2 steps) |
| `dt_ps` | scalar | CP time step in ps (3 a.u.) |
| `sampling_ps` | scalar | time between consecutive frames in ps (2 × `dt_ps` ≈ 0.145 fs) |
| `temperature_K` | scalar | average temperature (K) |
| `cell_A`, `volume_A3` | scalar | cubic cell edge (Å) and volume (Å³) |

The time in ps is `step * dt_ps`. The frames are evenly spaced, with no gaps or repeated steps. The production runs span 34.9 → 434.9 ps (454 K) and 39.3 → 439.3 ps (600 K).

## Usage

```python
import numpy as np

d = np.load('ab_initio-time_series/stress_454K.npz')
stress = d['stress']                   # (N, 6), GPa
time = d['step'] * d['dt_ps']          # ps
dt = float(d['sampling_ps'])           # ps between frames
off_diag = stress[:, 3:]               # xy, xz, yz: shear components for Green-Kubo
```

## Conversion from the original output

The original `cp.x` `.str` files print the full 3×3 tensor with 8 decimals. The tensor is symmetric to within the printed precision (1e-8 GPa), so only the 6 independent components are kept. They are stored as float32 to fit GitHub's file-size limit. The largest rounding error is 2e-6 GPa, and typical errors are ~1e-7 GPa, which is negligible for the viscosity.
