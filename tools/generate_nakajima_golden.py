"""Write the golden bulk melt values of the Nakajima et al. (2021) reference code.

The reference is https://github.com/mikinakajima/MeltScalingLaw at commit
2eca30f22073dfe089407d523f1cf76aaeebaad7 (MIT). Its source is patched in memory only:
``np.float`` becomes ``float`` (NumPy 1.24 and later), the colour maps and LaTeX are not
loaded, ``GG`` is set to the G of the table, the merging branch can be forced, and the bulk
values are returned right after their 1-sigma, before the clip and the melt-pool geometry.

Usage::

    python tools/generate_nakajima_golden.py <MeltScalingLaw checkout> tests/data/nakajima2021_golden.json
"""

from __future__ import annotations

import itertools
import json
import os
import sys
import types
from pathlib import Path

REFERENCE_COMMIT = '2eca30f22073dfe089407d523f1cf76aaeebaad7'
REFERENCE_G = 6.67408e-11
M_MARS = 6.4171e23
EDITS = [
    ("rc('text', usetex=True)\n", '\n', 1),
    ('.astype(np.float)', '.astype(float)', 2),
    (
        '        self.cm_data = np.loadtxt("vik/vik.txt")\n'
        "        self.vik_map = LinearSegmentedColormap.from_list('vik', self.cm_data)\n"
        '        self.cm_data2 = np.loadtxt("turku/turku.txt")\n'
        "        self.turku_map = LinearSegmentedColormap.from_list('turku', self.cm_data2)\n",
        '',
        1,
    ),
    (
        '        if self.vel <= critical_velocity:  # merging',
        '        if FORCE_MERGE or self.vel <= critical_velocity:  # merging',
        1,
    ),
    (
        '        f_model_std = f_model * dz\n',
        '        f_model_std = f_model * dz\n'
        '        return dict(f_raw=f_model, sigma=f_model_std, h=h_model, IE=IE_model,'
        ' m=Mantle_mass_model, dIE=IE_model * (dPE + dKE), Mt=Mt, Mi=Mi, Rt=Rt, Ri=Ri,'
        ' Rti=Rti, vesc=vesc, vcr=critical_velocity,'
        ' merging=bool(FORCE_MERGE or self.vel <= critical_velocity))\n',
        1,
    ),
]


def patched_model(reference: Path, force_merge: bool, gravity: float):
    """Return the reference ``Model`` class with the patch applied.

    Parameters
    ----------
    reference : Path
        Checkout of MeltScalingLaw at ``REFERENCE_COMMIT``.
    force_merge : bool
        Use the merging coefficients whatever the impact speed.
    gravity : float
        Gravitational constant [m^3 kg^-1 s^-2] the model uses.

    Returns
    -------
    type
        A ``Model`` subclass whose ``run_model`` returns the bulk values.

    Raises
    ------
    ValueError
        If the source does not hold an edited line exactly as often as expected.
    """
    src = (reference / 'melt_model.py').read_text()
    for old, new, count in EDITS:
        if src.count(old) != count:
            raise ValueError(f'melt_model.py holds {old[:40]!r} {src.count(old)} times')
        src = src.replace(old, new)
    module = types.ModuleType('melt_model_patched')
    module.FORCE_MERGE = force_merge
    exec(compile(src, str(reference / 'melt_model.py'), 'exec'), module.__dict__)

    class Model(module.Model):
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self.GG = gravity

    return Model


def rows(reference: Path, gravity: float, cases, force_merge: bool) -> list[dict]:
    """Run the patched reference on ``cases`` from inside the checkout."""
    model = patched_model(reference, force_merge, gravity)
    here = os.getcwd()
    os.chdir(reference)  # the reference reads its data files from the working directory
    try:
        out = []
        for theta, gamma, vel, mtotal in cases:
            d = model(Mtotal=mtotal, gamma=gamma, vel=vel, entropy0=3160, impact_angle=theta)
            values = d.run_model()
            out.append(
                {'theta': theta, 'gamma': gamma, 'vel': vel, 'Mtotal_Mmars': mtotal}
                | {k: float(v) for k, v in values.items()}
                | {'forced_merge': force_merge}
            )
        return out
    finally:
        os.chdir(here)


def main(argv: list[str]) -> int:
    """Write the golden table: the grid with Morrigan's G, and the issue table both ways.

    A grid case the reference sends to its hit-and-run branch is kept only with the
    merging branch forced, the branch the port uses for every merger.
    """
    from morrigan.constants import G, M_earth

    reference, target = Path(argv[1]), Path(argv[2])
    grid = list(
        itertools.product(
            (0, 30, 45, 60, 90), (0.03, 0.1, 0.3, 0.5), (1.0, 1.5), (2.0, 9.31, 40.0)
        )
    )
    natural = rows(reference, G, grid, force_merge=False)
    hit_and_run = [
        (r['theta'], r['gamma'], r['vel'], r['Mtotal_Mmars'])
        for r in natural
        if not r['merging']
    ]
    issue = list(itertools.product((0, 45, 90), (0.03, 0.1, 0.3), (1.0,), (M_earth / M_MARS,)))
    table = {
        'reference': f'https://github.com/mikinakajima/MeltScalingLaw/tree/{REFERENCE_COMMIT}',
        'M_Mars': M_MARS,
        'grid': [r for r in natural if r['merging']]
        + rows(reference, G, hit_and_run, force_merge=True),
        'issue_morrigan_G': rows(reference, G, issue, force_merge=False),
        'issue_reference_G': rows(reference, REFERENCE_G, issue, force_merge=False),
    }
    target.write_text(json.dumps(table, indent=1) + '\n')
    print(
        f'{len(table["grid"])} grid rows ({len(hit_and_run)} forced), {len(issue)} issue rows'
    )
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
