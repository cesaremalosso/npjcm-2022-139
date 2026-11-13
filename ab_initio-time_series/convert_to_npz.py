"""Convert the Car-Parrinello stress time series (.str) into compressed numpy archives.

Usage: python convert_to_npz.py stress_454K.str stress_600K.str
"""
import re
import sys

import numpy as np

AU_TIME_PS = 2.4188843265857e-5      # atomic unit of time in ps
DT_PS = 3 * AU_TIME_PS               # CP time step: 3 a.u. = 0.0725665 fs
CELL_A = 12.43                       # cubic cell edge (Angstrom)
COMPONENTS = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
IDX = ([0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2])


def load_str(fname):
    """Return (step, time_ps, stress[N,3,3]) from a cp.x-style .str file."""
    with open(fname) as fh:
        lines = [l for l in fh if not l.startswith('#')]
    assert len(lines) % 4 == 0, 'file must contain blocks of 4 lines'
    assert all(len(l.split()) == 2 for l in lines[0::4]), 'bad block header'
    a = np.array(''.join(lines).split(), dtype=float).reshape(-1, 11)
    return a[:, 0].astype(np.int64), a[:, 1], a[:, 2:].reshape(-1, 3, 3)


def convert(fname):
    step, time, s = load_str(fname)
    assert np.all(np.diff(step) == 2), 'steps are not uniformly spaced'
    assert np.abs(s - s.transpose(0, 2, 1)).max() < 1e-7, 'stress is not symmetric'
    stress = s[:, IDX[0], IDX[1]].astype(np.float32)
    temperature = int(re.search(r'(\d+)K', fname).group(1))
    out = fname.rsplit('.', 1)[0] + '.npz'
    np.savez_compressed(
        out,
        stress=stress,
        components=np.array(COMPONENTS),
        step=step.astype(np.int32),
        dt_ps=DT_PS,
        sampling_ps=2 * DT_PS,
        temperature_K=temperature,
        cell_A=CELL_A,
        volume_A3=CELL_A**3,
    )
    print(f'{fname} -> {out}: {len(step)} frames, '
          f'max float32 rounding error {np.abs(stress - s[:, IDX[0], IDX[1]]).max():.1e} GPa')


if __name__ == '__main__':
    for f in sys.argv[1:]:
        convert(f)
