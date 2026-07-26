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


def _pair(sep_au, ecc, a1_au=1.0, m_earths=1.0):
    """Build the mutable state arrays for a two-planet crossing call."""
    ap = np.array([a1_au, a1_au + sep_au]) * au2m
    mp = np.array([m_earths, m_earths]) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    e = np.array([ecc, ecc])
    interact = np.array([True, True])
    live = np.array([True, True])
    pid = np.arange(2)
    return ap, mp, rp, e, interact, live, pid


@pytest.mark.physics_invariant
def test_merger_kills_one_body_and_closes_its_record():
    """A forced merger removes exactly one body and its record closes
    the mass and speed books.

    Seed 9 on a nearly touching equal-mass pair lands in the collision
    branch (the collision probability is close to one there). The
    record must carry the pre-merge masses it was given, a merged mass
    equal to their plain sum, and a contact speed at or above the
    mutual escape speed of the recorded bodies, which discriminates any
    mass / radius mix-up in the bookkeeping.
    """
    np.random.seed(9)
    ap, mp, rp, e, interact, live, pid = _pair(0.005, 0.05)
    rec = orbit_cross_K25(ap, mp, rp, M_sun, 0.5, e, interact, live, 2, pid, 0)

    assert rec is not None
    assert int(np.sum(live)) == 1  # a merger removes exactly one body

    assert rec['M_target_before'] == pytest.approx(M_earth, rel=1e-12)
    assert rec['M_impactor_before'] == pytest.approx(M_earth, rel=1e-12)
    # Mass books: a merger is perfect, so the merged mass is the plain sum.
    assert rec['M_merged_after'] == pytest.approx(2.0 * M_earth, rel=1e-12)
    # Discrimination: reporting either body's own mass, or their mean, would
    # land a whole Earth mass away from the sum.
    assert abs(rec['M_merged_after'] - M_earth) > 0.5 * M_earth

    # Contact speed floors at the mutual escape speed of the recorded pair.
    v_esc = np.sqrt(2.0 * G * 2.0 * M_earth / (rec['R_target_before'] + rec['R_impactor']))
    assert rec['v_c'] >= v_esc * (1.0 - 1e-9)

    # The survivor's eccentricity is recorded on both sides of the merger and
    # both are bound. Reporting both is what lets a consumer following a planet
    # on a different orbit apply the change this collision made rather than
    # transplanting an absolute value that belongs to this body.
    assert 0.0 <= rec['e_after'] < 1.0
    assert 0.0 <= rec['e_before'] < 1.0

    # e_before is the target's eccentricity at the moment of the merger, read
    # from the state before merge_embryo overwrites it. That is not the value
    # the pair went in with: viscous stirring excites the orbits first, so the
    # pair entered at 0.05 and collides well above it. Both facts matter, since
    # recording the input value instead would describe a collision that never
    # happened at that geometry.
    assert rec['e_before'] > 0.05
    # Reading the same array one line too late would report the merged value,
    # which is the easy mistake because both come from `ecc[target_idx]`.
    assert rec['e_before'] != pytest.approx(rec['e_after'], rel=1e-6)


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
    ap, mp, rp, e, interact, live, pid = _pair(0.05, 0.05)
    a_before = ap.copy()
    m_before = mp.copy()
    rec = orbit_cross_K25(ap, mp, rp, M_sun, 0.5, e, interact, live, 2, pid, 0)

    assert rec is None
    assert bool(np.all(live))
    np.testing.assert_array_equal(mp, m_before)  # scattering moves orbits, not mass

    # Conservation: the mass-weighted orbit sum is exact.
    assert float(np.sum(mp * ap)) == pytest.approx(
        float(np.sum(m_before * a_before)), rel=1e-12
    )
    # Direction: inner in, outer out, by a nonzero amount.
    assert ap[0] < a_before[0]
    assert ap[1] > a_before[1]


@pytest.mark.physics_invariant
def test_ejection_removes_the_escaping_body_and_binds_the_survivor():
    """An ejection removes the escaper and leaves the survivor on the
    orbit that still reaches the encounter, or removes both if no such
    orbit is bound.

    The survivor absorbs the escaper's binding energy, so its orbit
    tightens to satisfy M_s/a_new = M_s/a_s + M_l/a_l exactly, and it
    must still pass through the place the encounter happened. That place
    now lies outside the tightened orbit, so it is the apocentre and the
    eccentricity follows as a_old/a_new - 1, equivalently
    (M_ejected/M_survivor)(a_survivor,0/a_escaper,0).

    That closed form says when nothing stays bound: the survivor must
    start outside the escaper by more than their mass ratio. Equipartition
    gives the lighter body the larger kick, so the survivor is the heavier
    one and must also be the outer one.

    Both regimes are pinned. An equal-mass pair at 30.0 and 30.3 au
    (seed 0) leaves a bound survivor, because the tie is broken by
    ejecting the outer body and the inner survivor is then inside the
    escaper. A 1.0 inner, 1.2 outer pair at 30 and 40 au (seed 0) crosses
    the threshold, 1.111 by the closed form, and both bodies go.
    """
    np.random.seed(0)
    ap, mp, rp, e, interact, live, pid = _pair(0.3, 0.3, a1_au=30.0, m_earths=100.0)
    a_before = ap.copy()
    rec = orbit_cross_K25(ap, mp, rp, M_sun, 0.5, e, interact, live, 2, pid, 0)

    assert rec is None
    # The tie sends the outer body out, so exactly one survives.
    assert list(live) == [True, False]
    # The survivor's orbit closes the pair's energy books exactly.
    assert ap[0] == pytest.approx(
        100.0 * M_earth / (100.0 * M_earth / a_before[0] + 100.0 * M_earth / a_before[1]),
        rel=1e-12,
    )
    # And its eccentricity is the apocentre closure, 30.0/30.3 short of 1.
    assert e[0] == pytest.approx(a_before[0] / ap[0] - 1.0, rel=1e-12)
    assert e[0] == pytest.approx(0.990099, rel=1e-5)
    assert 0.0 <= e[0] < 1.0
    # Discrimination: the pre-change form gives 0.502, a factor of two out.
    assert abs((1.0 - ap[0] / ((a_before[0] + a_before[1]) / 2)) - e[0]) > 1e-2

    # Nothing stays bound when the survivor starts outside the escaper by
    # more than the mass ratio: 1.0 inner and 1.2 outer at 30 and 40 au
    # gives (1.0/1.2)(40/30) = 1.111, so the survivor is unbound too.
    np.random.seed(0)
    ap2 = np.array([30.0, 40.0]) * au2m
    mp2 = np.array([1.0, 1.2]) * M_earth
    rp2 = np.array([planet_radius(m, _RHO) for m in mp2])
    e2 = np.array([0.3, 0.3])
    live2 = np.array([True, True])
    rec2 = orbit_cross_K25(
        ap2, mp2, rp2, M_sun, 0.5, e2, np.ones(2, dtype=bool), live2, 2, np.arange(2), 0
    )
    assert rec2 is None
    assert list(live2) == [False, False]
    assert e2[1] == pytest.approx((1.0 / 1.2) * (40.0 / 30.0), rel=1e-3)
    assert e2[1] >= 1.0


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
    ap, mp, rp, e, interact, live, pid = _pair(0.5, 0.001)
    before = [ap.copy(), mp.copy(), e.copy(), live.copy()]
    rec = orbit_cross_K25(ap, mp, rp, M_sun, 0.5, e, interact, live, 2, pid, 0)

    assert rec is None
    for now, then in zip([ap, mp, e, live], before):
        np.testing.assert_array_equal(now, then)

    # No random draw was consumed: the next draw equals a fresh one.
    after_draw = np.random.uniform()
    np.random.seed(5)
    fresh_draw = np.random.uniform()
    assert after_draw == pytest.approx(fresh_draw, rel=0.0, abs=0.0)


@pytest.mark.physics_invariant
def test_the_capped_fallback_keeps_the_eccentricity_a_body_already_carries():
    """When no draw can make the orbits overlap, the default does not
    discard a larger secular eccentricity.

    The rejection loop gives up when the crossing eccentricity exceeds
    twice the escape eccentricity, a condition fixed before the loop
    starts, and falls back to a geometric default of 1.001 times the
    crossing value. Two 1e-3 Earth-mass bodies at 1.0 and 1.5 au with
    e = 0.4 each (seed 0) put that ratio at 7.4, so the branch fires on
    the first iteration, and the carried 0.4 is the larger of the two
    against a default of 0.2002.

    Taking the larger is what every other site in the file does. Simply
    overwriting reports a contact speed of 7598.968 m/s against the
    15081.809 m/s the incoming orbits imply, a clean factor of two that
    feeds the consumer's erosion law directly.
    """
    np.random.seed(0)
    ap = np.array([1.0, 1.5]) * au2m
    mp = np.array([1.0e-3, 1.0e-3]) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    ecc = np.array([0.4, 0.4])
    live = np.array([True, True])
    rec = orbit_cross_K25(
        ap, mp, rp, M_sun, 0.5, ecc, np.ones(2, dtype=bool), live, 2, np.arange(2), 0
    )

    assert rec is not None
    # The carried eccentricity survives the fallback rather than being
    # replaced by the smaller geometric default.
    assert ecc[1] == pytest.approx(0.4, rel=1e-12)
    assert rec['v_c'] == pytest.approx(15081.809, rel=1e-6)
    # Discrimination: overwriting instead of taking the larger halves it.
    assert abs(7598.968 - rec['v_c']) > 100 * 1e-6 * rec['v_c']
    # The contact speed still floors at the mutual escape speed.
    v_esc = np.sqrt(2.0 * G * mp.sum() / (rec['R_target_before'] + rec['R_impactor']))
    assert rec['v_c'] >= v_esc * (1.0 - 1e-9)


@pytest.mark.physics_invariant
def test_the_recorded_pre_impact_eccentricity_belongs_to_the_target():
    """`e_before` is the surviving body's eccentricity, not the impactor's.

    Both bodies' eccentricities are excited before the collision and both live
    in the same array, so the record is one index away from describing the
    wrong body. On an equal-mass pair the two are identical and the mistake is
    invisible, which is why this uses a 1.0 and 0.3 Earth-mass pair: the
    rejection sampling scales each body's excitation by the square root of the
    other's mass, so the heavier target ends up on a different eccentricity
    from its impactor and the two can be told apart.

    Seed 9 puts this pair in the collision branch with body 0, the heavier one,
    as the target.
    """
    np.random.seed(9)
    ap = np.array([1.0, 1.005]) * au2m
    mp = np.array([1.0, 0.3]) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    e = np.array([0.05, 0.05])
    rec = orbit_cross_K25(
        ap,
        mp,
        rp,
        M_sun,
        0.5,
        e,
        np.array([True, True]),
        np.array([True, True]),
        2,
        np.arange(2),
        0,
    )

    assert rec is not None
    # The heavier body survives, so it is the one whose eccentricity is reported.
    assert rec['id_target'] == 0
    assert rec['M_target_before'] > rec['M_impactor_before']

    # Pinned: the target's own excited eccentricity at the moment of contact.
    # Reading the impactor's slot instead lands on a different value, which is
    # what an index mix-up would report.
    assert rec['e_before'] == pytest.approx(0.125086, rel=1e-4)
    assert 0.0 <= rec['e_before'] < 1.0

    # It is also not the post-merge value, which is what reading the array one
    # line later would give.
    assert rec['e_after'] == pytest.approx(0.043559, rel=1e-4)
    assert rec['e_before'] != pytest.approx(rec['e_after'], rel=1e-6)


@pytest.mark.physics_invariant
def test_a_merge_record_never_reports_an_unbound_orbit():
    """A record describes a closed orbit, and says so when it had to cap one.

    The excitation applied before a collision is drawn without an upper bound,
    and the fallback taken when no draw satisfies the overlap condition is not
    bounded either, so a body can reach the collision branch already past e = 1
    while the ejection branch treats that same condition as grounds for removing
    it. A consumer cannot use a record describing an open orbit: periapsis, the
    time-averaged separation and the Hill radius all assume a closed one, and
    the coupled framework refuses a whole impact history over it.

    The cap is therefore applied where the record is written, not to the model
    state, and it warns, because a capped record no longer describes the
    geometry the model held and an impact reported that way is not trustworthy.
    """
    from morrigan.orbit_cross_K25 import MAX_RECORDED_ECC, _bounded_eccentricity

    # An already-closed orbit passes through untouched, including the boundary.
    assert _bounded_eccentricity(0.0, 'before', 0) == 0.0
    assert _bounded_eccentricity(0.5, 'before', 0) == pytest.approx(0.5, rel=1e-12)
    assert _bounded_eccentricity(MAX_RECORDED_ECC, 'after', 0) == pytest.approx(
        MAX_RECORDED_ECC, rel=1e-12
    )

    # An open orbit is capped and reported, rather than passed to a consumer
    # that would reject the whole timeline over it.
    with pytest.warns(UserWarning, match='not a closed orbit'):
        capped = _bounded_eccentricity(1.815, 'before', 3)
    assert capped == pytest.approx(MAX_RECORDED_ECC, rel=1e-12)
    assert 0.0 <= capped < 1.0

    # The warning names the side of the collision and the body, so an ensemble
    # run can be traced back to the impact that produced it.
    with pytest.warns(UserWarning, match='Body 7'):
        _bounded_eccentricity(2.0, 'after', 7)
    with pytest.warns(UserWarning, match='after the collision'):
        _bounded_eccentricity(2.0, 'after', 7)

    # Discrimination: capping at 1.0 rather than below it would still fail the
    # consumer's half-open bound, which is the whole point of the cap.
    assert MAX_RECORDED_ECC < 1.0
