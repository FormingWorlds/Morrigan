"""Tests for ``src/morrigan/orbit_cross_K25.py``.

``orbit_cross_K25`` resolves a scheduled crossing into one of four
outcomes: nothing (a stable pair), a merger, a scattering, or an
ejection. Each path is driven deterministically here through the global
numpy seed and pinned to its conservation law: the merger to its record
schema and mass books, the scattering to the exact conservation of the
mass-weighted semi-major axis (eqs. 18-20 of Kimura et al. 2025), and
the ejection to removing exactly the body excited past e = 1 while the
survivor stays bound.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import G, M_earth, M_sun, au2m
from morrigan.helper_functions import planet_radius
from morrigan.orbit_cross_K25 import orbit_cross_K25

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0


def _pair(sep_au, ecc, a1_au=1.0, m_earths=1.0, f_atm=0.01):
    """Build the mutable state arrays for a two-planet crossing call."""
    ap = np.array([a1_au, a1_au + sep_au]) * au2m
    mp = np.array([m_earths, m_earths]) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    e = np.array([ecc, ecc])
    f = np.array([f_atm, f_atm])
    interact = np.array([True, True])
    live = np.array([True, True])
    pid = np.arange(2)
    return ap, mp, rp, f, e, interact, live, pid


@pytest.mark.physics_invariant
def test_merger_kills_one_body_and_closes_its_record():
    """A forced merger removes exactly one body and its record closes
    the mass and speed books.

    Seed 9 on a nearly touching equal-mass pair lands in the collision
    branch (the collision probability is close to one there). The
    record must carry the pre-merge masses it was given, a merged mass
    equal to the sum minus the atmosphere actually lost, a loss
    fraction inside [0, 1], and a contact speed at or above the mutual
    escape speed of the recorded bodies, which discriminates any
    mass / radius mix-up in the bookkeeping.
    """
    np.random.seed(9)
    ap, mp, rp, f, e, interact, live, pid = _pair(0.005, 0.05)
    rec = orbit_cross_K25(ap, mp, rp, M_sun, f, 0.5, e, interact, live, 2, pid, 0)

    assert rec is not None
    assert int(np.sum(live)) == 1  # a merger removes exactly one body

    assert rec['M_target_before'] == pytest.approx(M_earth, rel=1e-12)
    assert rec['M_impactor_before'] == pytest.approx(M_earth, rel=1e-12)
    # Mass books: merged mass = sum - lost atmosphere, exactly.
    atm_combined = 2.0 * M_earth * 0.01
    expected = 2.0 * M_earth - rec['atm_mass_loss_frac'] * atm_combined
    assert rec['M_merged_after'] == pytest.approx(expected, rel=1e-12)
    assert 0.0 <= rec['atm_mass_loss_frac'] <= 1.0

    # Contact speed floors at the mutual escape speed of the recorded pair.
    v_esc = np.sqrt(2.0 * G * 2.0 * M_earth / (rec['R_target_before'] + rec['R_impactor']))
    assert rec['v_c'] >= v_esc * (1.0 - 1e-9)

    # The survivor's post-merge eccentricity is recorded and bound.
    assert 0.0 <= rec['e_after'] < 1.0


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_scattering_conserves_the_mass_weighted_orbit_sum():
    """A scattering event conserves sum(M a) exactly and pushes the pair
    apart, following eqs. 18-20 of Kimura et al. (2025).

    Analytical limit: eqs. 19-20 shift the two orbits by amounts
    weighted with the opposite body's mass fraction, so M1 da1 + M2 da2
    = 0 identically and the mass-weighted semi-major-axis sum is
    conserved to machine precision. Seed 0 on a well-separated pair
    lands in the scattering branch (collision probability near zero).
    The inner orbit must move in and the outer out, both bodies stay
    alive with unchanged masses, and no merge record is produced. A
    formula scattering both orbits the same way, the plausible wrong
    transcription, would break the conservation at the percent level.
    """
    np.random.seed(0)
    ap, mp, rp, f, e, interact, live, pid = _pair(0.05, 0.05)
    a_before = ap.copy()
    m_before = mp.copy()
    rec = orbit_cross_K25(ap, mp, rp, M_sun, f, 0.5, e, interact, live, 2, pid, 0)

    assert rec is None
    assert bool(np.all(live))
    np.testing.assert_array_equal(mp, m_before)  # scattering moves orbits, not mass

    # Conservation: the mass-weighted orbit sum is exact.
    assert float(np.sum(mp * ap)) == pytest.approx(float(np.sum(m_before * a_before)), rel=1e-12)
    # Direction: inner in, outer out, by a nonzero amount.
    assert ap[0] < a_before[0]
    assert ap[1] > a_before[1]


@pytest.mark.physics_invariant
def test_ejection_removes_the_escaping_body_and_binds_the_survivor():
    """An ejection kills exactly the body excited past e = 1 and leaves
    the survivor on a bound orbit that closes the pair's energy books.

    A compact massive pair at 30 au has a mutual escape eccentricity
    far above one, so the scattering excitation drives an ejection
    (seed 0). Exactly one body must die; the survivor must carry an
    eccentricity inside [0, 1), sit inside its own original orbit
    (it absorbs the escaper's binding energy), and satisfy the exact
    energy closure M_s/a_new = M_s/a_s + M_l/a_l of the survivor
    formula. This is the regression pin for the ejection bookkeeping:
    flagging the wrong body would leave a hyperbolic e >= 1 orbit
    alive, which the bound below rejects. The unequal-mass case
    (1 and 50 Earth masses at 1.0 and 1.5 au, seed 0) pins the
    ordering: the lighter body takes the larger kick and escapes, and
    the heavy survivor can land between the two original orbits, so
    only its own original orbit bounds it.
    """
    np.random.seed(0)
    ap, mp, rp, f, e, interact, live, pid = _pair(0.3, 0.3, a1_au=30.0, m_earths=100.0, f_atm=0.0)
    a_before = ap.copy()
    m_before = mp.copy()
    rec = orbit_cross_K25(ap, mp, rp, M_sun, f, 0.5, e, interact, live, 2, pid, 0)

    assert rec is None
    assert int(np.sum(live)) == 1  # exactly one body escapes

    survivor = int(np.argmax(live))
    ejected = 1 - survivor
    assert 0.0 <= e[survivor] < 1.0
    assert 0.0 < ap[survivor] < a_before[survivor]
    # Exact energy closure of the survivor formula.
    assert mp[survivor] / ap[survivor] == pytest.approx(
        m_before[survivor] / a_before[survivor] + m_before[ejected] / a_before[ejected],
        rel=1e-12,
    )

    # Unequal masses: the lighter body escapes, the heavy survivor stays
    # bound between the original orbits, and the energy books still close.
    np.random.seed(0)
    ap2 = np.array([1.0, 1.5]) * au2m
    mp2 = np.array([1.0, 50.0]) * M_earth
    rp2 = np.array([planet_radius(m, _RHO) for m in mp2])
    e2 = np.array([0.3, 0.3])
    live2 = np.array([True, True])
    rec2 = orbit_cross_K25(ap2, mp2, rp2, M_sun, np.zeros(2), 0.5, e2,
                           np.ones(2, dtype=bool), live2, 2, np.arange(2), 0)
    assert rec2 is None
    assert list(live2) == [False, True]  # the 1 Earth-mass body is ejected
    assert 0.0 <= e2[1] < 1.0
    assert ap2[1] == pytest.approx(
        50.0 * M_earth / (50.0 * M_earth / (1.5 * au2m) + 1.0 * M_earth / (1.0 * au2m)),
        rel=1e-12,
    )
    # The survivor sits inside its own orbit but outside the escaper's,
    # which discriminates the true bound from the min-of-both wrong bound.
    assert 1.0 * au2m < ap2[1] < 1.5 * au2m


@pytest.mark.physics_invariant
def test_a_stable_wide_pair_is_left_untouched():
    """A Jacobi-stable, non-overlapping two-planet system passes through
    the crossing resolver unchanged.

    Guard branch: a wide near-circular pair fails both trigger
    conditions, so the resolver must return early with no record and
    no mutation of any array, which is the do-no-harm contract the
    driver relies on when a scheduled event turns out benign. Seeded
    with 5 to pin that not even a random draw is consumed on this path
    (the follow-up draw matches a fresh seed-5 draw exactly).
    """
    np.random.seed(5)
    ap, mp, rp, f, e, interact, live, pid = _pair(0.5, 0.001)
    before = [ap.copy(), mp.copy(), e.copy(), live.copy(), f.copy()]
    rec = orbit_cross_K25(ap, mp, rp, M_sun, f, 0.5, e, interact, live, 2, pid, 0)

    assert rec is None
    for now, then in zip([ap, mp, e, live, f], before):
        np.testing.assert_array_equal(now, then)

    # No random draw was consumed: the next draw equals a fresh one.
    after_draw = np.random.uniform()
    np.random.seed(5)
    fresh_draw = np.random.uniform()
    assert after_draw == pytest.approx(fresh_draw, rel=0.0, abs=0.0)
