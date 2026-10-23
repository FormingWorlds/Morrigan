"""Tests for ``src/morrigan/impact_melt.py``: the bulk mantle melt fraction of a merger.

Invariants exercised: the port equals the Nakajima et al. (2021) reference code to 1e-12
on a grid of angles, mass ratios, speeds and masses (golden table from
``tools/generate_nakajima_golden.py``); the melt fraction stays in [0, 1] and falls from
head-on to grazing; the fit-range and critical-velocity flags switch at each boundary;
every argument outside its domain is refused by name.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from morrigan.constants import G, M_earth
from morrigan.impact_melt import M_MARS, melt_fraction

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

GOLDEN = json.loads((Path(__file__).parent / 'data' / 'nakajima2021_golden.json').read_text())


def _call(row):
    """Run the port on one golden row, fed the reference code's own radii and G."""
    return melt_fraction(
        row['Mt'],
        row['Mi'],
        row['Rt'],
        row['Ri'],
        row['vel'] * row['vesc'],
        math.sin(math.radians(row['theta'])),
        R_merged=row['Rti'],
    )


def _case(theta, gamma=0.1, vel=1.0, m_total=M_earth, density=5500.0):
    """Inputs for a merger of total mass ``m_total`` at one shared density."""
    m_i, m_t = gamma * m_total, (1.0 - gamma) * m_total
    r_t, r_i = ((3.0 * m / (4.0 * math.pi * density)) ** (1.0 / 3.0) for m in (m_t, m_i))
    v_esc = math.sqrt(2.0 * G * m_total / (r_t + r_i))
    return m_t, m_i, r_t, r_i, vel * v_esc, math.sin(math.radians(theta))


_NAMES = 'M_target M_impactor R_target R_impactor v_impact impact_parameter'.split()
_ARGS = dict(zip(_NAMES, _case(45.0)))


@pytest.mark.reference_pinned
@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'row', GOLDEN['grid'], ids=lambda r: '{theta}-{gamma}-{vel}-{Mtotal_Mmars}'.format(**r)
)
def test_port_equals_the_reference_code(row):
    """Cross-implementation check: f, its 1-sigma, dIE and h equal MeltScalingLaw at
    2eca30f2 (merging coefficients; rows the reference sends to the hit-and-run branch are
    compared with its merging branch forced, and the port flags them above v_cr)."""
    res = _call(row)
    assert res.f_melt_raw == pytest.approx(row['f_raw'], rel=1e-12, abs=0.0)
    assert res.f_melt == pytest.approx(min(row['f_raw'], 1.0), rel=1e-12, abs=0.0)
    assert res.f_melt_sigma == pytest.approx(row['sigma'], rel=1e-12, abs=0.0)
    assert res.dIE == pytest.approx(row['dIE'], rel=1e-12, abs=0.0)
    assert res.h == pytest.approx(row['h'], rel=1e-12, abs=0.0)
    assert res.above_v_cr is row['forced_merge']
    # Discrimination: the latent heat in place of E_M, the 0.70 mantle fraction left out,
    # and dIE without its angle fit all land far from the pinned value.
    if row['f_raw'] < 1.0:
        assert abs(res.f_melt - row['f_raw'] * 5.2e6 / 7.18e5) > 1e-3
        assert abs(res.f_melt - row['f_raw'] * 0.70) > 1e-3
        assert abs(res.dIE - row['dIE'] / row['IE']) > 1e-3 * res.dIE


@pytest.mark.reference_pinned
@pytest.mark.parametrize(
    'row', GOLDEN['grid'], ids=lambda r: '{theta}-{gamma}-{vel}-{Mtotal_Mmars}'.format(**r)
)
def test_above_v_cr_flips_at_the_reference_critical_speed(row):
    """Just above the reference code's critical speed (in units of the mutual escape
    speed, recomputed here with Morrigan's G) the merger is flagged, just below it is not,
    so every term of the Genda et al. (2012) fit is pinned at every golden row."""
    v_esc = math.sqrt(2.0 * G * (row['Mt'] + row['Mi']) / (row['Rt'] + row['Ri']))
    b = math.sin(math.radians(row['theta']))
    for factor, expected in ((1.0 + 1e-6, True), (1.0 - 1e-6, False)):
        res = melt_fraction(
            row['Mt'],
            row['Mi'],
            row['Rt'],
            row['Ri'],
            row['vcr'] * factor * v_esc,
            b,
            R_merged=row['Rti'],
        )
        assert res.above_v_cr is expected, (factor, row['vcr'])


@pytest.mark.reference_pinned
def test_issue_table_at_one_earth_mass():
    """The table of Morrigan issue 17 (1 M_E, v_imp = v_esc, reference radii) to its two
    decimals; the 90 deg, gamma 0.3 sigma is 0.3844 in the reference code itself, so that
    cell is pinned to the code (0.38), not to the issue text (0.39)."""
    text = {
        0: ('0.36 ± 0.02', '1.00 ± 0.07', '1.00 ± 0.16'),
        45: ('0.29 ± 0.04', '0.87 ± 0.12', '1.00 ± 0.27'),
        90: ('0.09 ± 0.06', '0.27 ± 0.17', '0.60 ± 0.38'),
    }
    for row in GOLDEN['issue_morrigan_G']:
        res = _call(row)
        assert res.f_melt == pytest.approx(min(row['f_raw'], 1.0), rel=1e-12, abs=0.0)
        cell = f'{res.f_melt:.2f} ± {res.f_melt_sigma:.2f}'
        assert cell == text[row['theta']][(0.03, 0.1, 0.3).index(row['gamma'])]
    assert len(GOLDEN['issue_morrigan_G']) == 9


@pytest.mark.physics_invariant
@pytest.mark.parametrize('gamma', [0.03, 0.1, 0.3, 0.5])
def test_melt_falls_from_head_on_to_grazing(gamma):
    """At fixed masses and speed the melt fraction never rises with the impact angle, over
    the continuous angle and not only at the simulated ones; it stays in [0, 1]."""
    f = [melt_fraction(*_case(theta, gamma)).f_melt for theta in np.linspace(0.0, 90.0, 181)]
    assert all(0.0 <= v <= 1.0 for v in f)
    assert all(b <= a for a, b in zip(f, f[1:]))
    # Strictly lower at grazing than head-on unless both are clipped at 1.
    assert f[-1] < f[0] or f[0] == pytest.approx(1.0, abs=0.0)


@pytest.mark.physics_invariant
@pytest.mark.parametrize('tie', [15.0, 37.5, 52.5, 75.0])
def test_sigma_at_a_tie_uses_the_lower_simulated_angle(tie):
    """Between the simulated angles the 1-sigma uses the errors of the nearest one; at an
    exact tie, with the impact parameter formed as the driver forms it, the lower one. The
    relative error is continuous from below at the tie and jumps just above it."""

    def rel(theta):
        args = _case(theta, 0.03)[:5]
        res = melt_fraction(*args, float(np.sin(np.deg2rad(theta))))
        assert res.f_melt < 1.0  # unclipped, so sigma / f is the relative error
        return res.f_melt_sigma / res.f_melt

    assert rel(tie) == pytest.approx(rel(tie - 1e-7), rel=1e-6)
    assert abs(rel(tie + 1e-7) / rel(tie) - 1.0) > 0.01
    # 1e-11 deg above the tie is past the round-off of asin(sin), so the upper angle
    assert abs(rel(tie + 1e-11) / rel(tie) - 1.0) > 0.01


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    ('kwargs', 'expected'),
    [
        ({'m_total': 1.0001 * M_MARS}, True),
        ({'m_total': 0.9999 * M_MARS}, False),
        ({'m_total': 52.99 * M_MARS}, True),
        ({'m_total': 53.01 * M_MARS}, False),
        ({'gamma': 0.0301}, True),
        ({'gamma': 0.0299}, False),
        ({'gamma': 0.5}, True),
        ({'vel': 1.0}, True),
        ({'vel': 1.0 - 1e-12}, True),
        ({'vel': 1.0 - 1e-6}, False),
        ({'vel': 1.9999}, True),
        ({'vel': 2.0001}, False),
    ],
    ids=lambda v: str(v),
)
def test_in_range_flag_at_each_boundary(kwargs, expected):
    """The fit range (1 to 53 M_Mars, gamma 0.03 to 0.5, v_esc to 2 v_esc) flags a merger
    just outside each bound and accepts it just inside; the lower speed bound carries a
    1e-9 relative margin, so a contact speed that equals the recomputed v_esc stays in."""
    res = melt_fraction(*_case(45.0, **({'gamma': 0.1} | kwargs)))
    assert res.in_range is expected
    assert 0.0 <= res.f_melt <= 1.0


@pytest.mark.physics_invariant
def test_contact_speed_at_the_mutual_escape_speed_is_in_range():
    """A merger whose speed is the mutual escape speed, computed the way merge_embryo does,
    is in range: the recomputed v_esc must not flag it by round-off."""
    m_t, m_i, r_t, r_i = 0.8 * M_earth, 0.2 * M_earth, 6.0e6, 3.7e6
    v_esc = np.sqrt(2 * G * (m_t + m_i) / (r_t + r_i))
    res = melt_fraction(m_t, m_i, r_t, r_i, v_esc, 0.5)
    assert res.in_range is True
    assert res.above_v_cr is False


@pytest.mark.physics_invariant
def test_above_v_cr_on_both_sides_of_the_critical_speed():
    """Grazing at 1.5 v_esc is above the Genda et al. (2012) critical speed (1.047 v_esc
    at gamma 0.1), head-on at 1.5 v_esc is below it (4.463 v_esc)."""
    assert melt_fraction(*_case(90.0, vel=1.5)).above_v_cr is True
    assert melt_fraction(*_case(0.0, vel=1.5)).above_v_cr is False
    assert melt_fraction(*_case(90.0, vel=1.0)).above_v_cr is False


@pytest.mark.physics_invariant
def test_merged_radius_defaults_to_the_target_density():
    """Without R_merged the merged body keeps the target's density, as Morrigan's own
    radius update does; at one shared density this is the volume sum."""
    m_t, m_i, r_t, r_i, v, b = _case(30.0, 0.2)
    shared = melt_fraction(m_t, m_i, r_t, r_i, v, b)
    explicit = melt_fraction(m_t, m_i, r_t, r_i, v, b, R_merged=(r_t**3 + r_i**3) ** (1 / 3))
    assert shared.dIE == pytest.approx(explicit.dIE, rel=1e-12)
    # A merged radius 1 percent larger releases visibly less energy.
    bigger = melt_fraction(
        m_t, m_i, r_t, r_i, v, b, R_merged=1.01 * (r_t**3 + r_i**3) ** (1 / 3)
    )
    assert bigger.dIE < shared.dIE * (1 - 1e-3)


@pytest.mark.parametrize(
    ('name', 'bad'),
    [
        ('M_target', 0.0),
        ('M_target', -1.0),
        ('M_target', float('nan')),
        ('M_target', float('inf')),
        ('M_target', True),
        ('M_impactor', 0.0),
        ('R_target', -1.0),
        ('R_impactor', float('nan')),
        ('v_impact', 0.0),
        ('v_impact', float('inf')),
        ('impact_parameter', -0.1),
        ('impact_parameter', 1.1),
        ('impact_parameter', float('nan')),
        ('impact_parameter', False),
        ('R_merged', 0.0),
    ],
)
def test_argument_outside_its_domain_is_refused_by_name(name, bad):
    """A non-finite, non-positive, boolean or out-of-range argument raises ValueError naming
    it, so a broken merge record stops the run instead of giving a melt fraction."""
    with pytest.raises(ValueError, match=f'{name} '):
        melt_fraction(**(_ARGS | {name: bad}))
    assert 0.0 <= melt_fraction(**_ARGS).f_melt <= 1.0


def test_merged_radius_that_releases_no_energy_is_refused():
    """A merged radius so large that the merger would absorb energy is a caller error:
    ValueError names R_merged, and the default merged radius releases energy."""
    m_t, m_i, r_t, r_i, v, b = _case(45.0, 0.2)
    with pytest.raises(ValueError, match='R_merged .* release no finite positive energy'):
        melt_fraction(m_t, m_i, r_t, r_i, v, b, R_merged=100.0 * r_t)
    assert melt_fraction(m_t, m_i, r_t, r_i, v, b).dIE > 0.0


@pytest.mark.parametrize(
    ('override', 'phrase'),
    [
        ({'v_impact': 1e160}, 'overflow or underflow the energy or the escape speed'),
        ({'M_target': 10**400}, 'M_target must be'),
        (
            {'M_target': 1e160, 'M_impactor': 1e159},
            'overflow or underflow the energy or the escape speed',
        ),
        ({'R_merged': 1e-320}, 'overflow or underflow the energy or the escape speed'),
        (
            {'R_target': 1e308, 'R_impactor': 1e308},
            'overflow or underflow the energy or the escape speed',
        ),
        (
            {
                'M_target': 1e-4,
                'M_impactor': 1e-4,
                'R_target': 5e-324,
                'R_impactor': 5e-324,
                'v_impact': 1.0,
            },
            'overflow or underflow the melt fraction',
        ),
        (
            {
                'M_target': 1.0,
                'M_impactor': 1.0,
                'R_target': 1e307,
                'R_impactor': 1e307,
                'v_impact': 1e-160,
            },
            'overflow or underflow the melt fraction',
        ),
    ],
    ids=[
        'speed',
        'int-mass',
        'float-mass',
        'tiny-merged-radius',
        'huge-radii',
        'f-raw',
        'f-raw-underflow',
    ],
)
def test_inputs_that_overflow_are_refused_by_name(override, phrase):
    """Finite inputs whose energy, escape speed or melt fraction overflows, underflows or
    divides by zero in floating point raise ValueError with the phrase of their branch,
    not an arithmetic error or the RuntimeError kept for broken coefficient files."""
    with pytest.raises(ValueError, match=phrase):
        melt_fraction(**(_ARGS | override))
    assert melt_fraction(**_ARGS).dIE > 0.0


def test_merged_radius_is_named_only_when_given():
    """The overflow message names R_merged only when the caller passed one."""
    with pytest.raises(ValueError) as given:
        melt_fraction(**(_ARGS | {'R_merged': 1e-320}))
    with pytest.raises(ValueError) as default:
        melt_fraction(**(_ARGS | {'v_impact': 1e160}))
    assert 'R_merged' in str(given.value) and 'R_merged' not in str(default.value)


@pytest.mark.physics_invariant
def test_impactor_heavier_than_target_is_refused():
    """Morrigan's target is the larger body; swapped masses are refused, not silently used."""
    m_t, m_i, r_t, r_i, v, b = _case(45.0, 0.3)
    with pytest.raises(ValueError, match='M_impactor'):
        melt_fraction(m_i, m_t, r_i, r_t, v, b)
    assert melt_fraction(m_t, m_i, r_t, r_i, v, b).f_melt > 0.0


@pytest.mark.physics_invariant
def test_coefficient_files_of_the_wrong_shape_are_refused(tmp_path, monkeypatch):
    """A parameter.txt or Model_sigma.txt of the wrong shape stops the run with both
    shapes named, instead of fitting with misplaced coefficients."""
    from morrigan import impact_melt as im

    data = tmp_path / 'data' / 'nakajima2021'
    data.mkdir(parents=True)
    (data / 'parameter.txt').write_text('1.0 2.0 3.0\n')
    (data / 'Model_sigma.txt').write_text('#0 30 45 60 90\n1 2 3 4 5\n1 2 3 4 5\n1 2 3 4 5\n')
    monkeypatch.setattr(im.resources, 'files', lambda package: tmp_path)
    with pytest.raises(RuntimeError, match=r'shapes \(3,\), \(3, 5\)'):
        im._coefficients.__wrapped__()
    # The shipped files have the shapes the fit needs.
    monkeypatch.undo()
    para, sigma = im._coefficients.__wrapped__()
    assert para.shape == (12,) and sigma.shape == (3, 5)


@pytest.mark.physics_invariant
@pytest.mark.parametrize(
    'broken', ['negative h', 'negative IE', 'negative mantle', 'infinite sigma']
)
def test_broken_coefficients_stop_the_run(broken, monkeypatch):
    """A non-positive fit factor or a non-finite 1-sigma raises RuntimeError rather than
    reporting a melt fraction; the shipped coefficients never reach either guard."""
    from morrigan import impact_melt as im

    para, sigma = (a.copy() for a in im._coefficients())
    if broken != 'infinite sigma':
        para[{'negative h': 0, 'negative IE': 3, 'negative mantle': 10}[broken]] = -5.0
        match = 'non-positive fit factor'
    else:
        sigma[:] = np.inf
        match = 'leaves its bounds'
    monkeypatch.setattr(im, '_coefficients', lambda: (para, sigma))
    with pytest.raises(RuntimeError, match=match):
        melt_fraction(*_case(45.0))
    monkeypatch.undo()
    assert 0.0 < melt_fraction(*_case(45.0)).f_melt <= 1.0
