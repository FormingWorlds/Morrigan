"""Tests for ``src/morrigan/merge_embryo.py``.

``merge_embryo`` is the merger bookkeeping: it must close the mass
books (mass conserved exactly on the plain sum), place
the merged body on the eq. 16 orbit of Kimura et al. (2025), bound the
eq. 17 eccentricity, kill exactly the smaller body, and clamp the
merged eccentricity inside its composition bounds. ``collision_velocity`` must floor
at the mutual escape speed, the analytic zero-eccentricity limit of the
sqrt(v_inf^2 + v_esc^2) contact-speed convention.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import G, M_earth, M_sun, au2m
from morrigan.helper_functions import planet_radius
from morrigan.merge_embryo import collision_velocity, merge_embryo

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0


def _pair(m1=2.0, m2=1.0, a1=1.0, a2=1.2, e1=0.05, e2=0.08):
    """Return the mutable arrays merge_embryo works on, in model units."""
    ap = np.array([a1, a2]) * au2m
    mp = np.array([m1, m2]) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    ecc = np.array([e1, e2])
    live = np.array([True, True])
    return ap, mp, rp, ecc, live


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_collision_velocity_floors_at_the_mutual_escape_speed():
    """The contact speed is sqrt(v_inf^2 + v_esc^2) and so equals the
    mutual escape speed exactly for circular orbits.

    Analytical limit: with both eccentricities zero the relative speed
    at infinity vanishes and the contact speed must equal
    sqrt(2 G (M1 + M2) / (R1 + R2)) to machine precision; this is the
    convention the Kegerreis et al. (2020) loss law expects. The
    eccentric case is pinned by hand through v_inf = e_ij v_kep:
    11945.7 m/s for the reference pair, which discriminates the wrong
    linear combination v_inf + v_esc (14317.5 m/s there). Speed grows
    monotonically with eccentricity.
    """
    ap, mp, rp, _, _ = _pair()
    v_esc = np.sqrt(2.0 * G * (mp[0] + mp[1]) / (rp[0] + rp[1]))

    v_circ = collision_velocity(list(ap), list(mp), list(rp), M_sun, [0.0, 0.0])
    assert v_circ == pytest.approx(v_esc, rel=1e-12)

    v_ecc = collision_velocity(list(ap), list(mp), list(rp), M_sun, [0.05, 0.08])
    assert v_ecc == pytest.approx(11945.706, rel=1e-6)
    # Discrimination: the linear-sum wrong form v_inf + v_esc evaluates to
    # 14317.5 m/s here, twenty percent above the quadrature value.
    e_ij = np.sqrt(0.05**2 + 0.08**2)
    v_inf = e_ij * np.sqrt(G * M_sun / (1.1 * au2m))
    assert v_inf + v_esc == pytest.approx(14317.5, rel=1e-4)
    assert abs((v_inf + v_esc) - v_ecc) > 100 * 1e-6 * v_ecc

    # Monotone in eccentricity, and always at or above the escape floor.
    v_more = collision_velocity(list(ap), list(mp), list(rp), M_sun, [0.1, 0.16])
    assert v_more > v_ecc > v_circ


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_merger_closes_its_mass_and_orbit_books():
    """A merger conserves mass exactly and lands on the eq. 16 orbit.

    Seeded with 7 (the pericentre-alignment draw). The merged
    semi-major axis is the orbital-energy-weighted form of eq. 16 of
    Kimura et al. (2025): a = (M1 + M2) / (M1/a1 + M2/a2), by hand
    1.05882 au for 2 M_e at 1.0 au absorbing 1 M_e at 1.2 au. Both the
    arithmetic mean (1.1 au) and the mass-weighted mean (1.0667 au),
    the two plausible wrong forms, sit far outside the pin. The merger
    is perfect, so the surviving body carries the plain sum of the two
    masses and nothing is shed.
    """
    ap, mp, rp, ecc, live = _pair()
    # merge_embryo mutates these arrays in place, so capture the pre-merge
    # masses before the call rather than reading them back afterwards.
    mass_before = float(mp.sum())
    target_before = float(mp[0])

    np.random.seed(7)
    v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
    ap_n, mp_n, ecc_n, live_n = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5)

    a_expected = 3.0 * M_earth / (2.0 * M_earth / (1.0 * au2m) + 1.0 * M_earth / (1.2 * au2m))
    assert ap_n[0] == pytest.approx(a_expected, rel=1e-12)
    # Discrimination: arithmetic and mass-weighted means both fail the pin.
    assert abs(1.1 * au2m - a_expected) > 1e6
    assert abs((2.0 * 1.0 + 1.0 * 1.2) / 3.0 * au2m - a_expected) > 1e6

    # Mass is conserved exactly: the survivor is the plain sum.
    assert float(mp_n[0]) == pytest.approx(mass_before, rel=1e-12)
    # Discrimination: the survivor gains the whole impactor, so the merged
    # mass differs from the pre-merge target by a third of the total here.
    # Reporting the target's own mass unchanged would fail this.
    assert abs(float(mp_n[0]) - target_before) > 0.3 * mass_before

    # The smaller body is dead.
    assert list(live_n) == [True, False]


@pytest.mark.physics_invariant
def test_survivor_selection_is_by_mass_not_position():
    """The heavier body survives regardless of array order, and an
    equal-mass tie keeps the first slot.

    Swapping the pair so the heavier body sits in slot 1 must flip the
    dead flag to slot 0; an exact tie resolves to slot 0 by the >=
    convention. Seeded with 11. This guards the target / impactor
    indexing against position-based selection, which would pass any
    test that always puts the heavy body first.
    """
    np.random.seed(11)
    ap, mp, rp, ecc, live = _pair(m1=1.0, m2=2.0)
    mass_before = float(mp.sum())
    v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
    _, mp_n, _, live_n = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5)
    assert list(live_n) == [False, True]
    # The survivor, whichever slot it is in, carries the whole pair's mass.
    assert float(mp_n[1]) == pytest.approx(mass_before, rel=1e-12)

    np.random.seed(11)
    ap, mp, rp, ecc, live = _pair(m1=1.0, m2=1.0)
    v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
    _, _, _, live_tie = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5)
    assert list(live_tie) == [True, False]


@pytest.mark.physics_invariant
def test_the_merged_eccentricity_stays_within_its_composition_bounds():
    """A far-above-escape impact still conserves mass, and the merged
    eccentricity stays inside its vector-composition bounds.

    The eq. 17 eccentricity is a vector composition, so across seeds it
    must stay within [|M1 e1 - M2 e2|, M1 e1 + M2 e2] / (M1 + M2). The
    sweep uses a close pair (1.0 and 1.05 au) whose pericentre-alignment
    cosine is interior to (-1, 1), so the drawn alignment genuinely
    varies with the seed; seeds 0-9 must therefore produce a
    non-degenerate spread of merged eccentricities inside the bounds. The
    edge cases run on inputs the function actually reads: a pair of
    circular orbits, where eq. 17 must compose exactly zero, and a
    million-to-one mass ratio, where the merged orbit and eccentricity
    must collapse onto the heavy body's own.
    """
    # Circular pair: the vector composition of two zero eccentricities is
    # exactly zero, whatever alignment is drawn.
    ap, mp, rp, ecc, live = _pair(m1=1.0, m2=1.0, e1=0.0, e2=0.0)
    mass_before = float(mp.sum())
    np.random.seed(0)
    _, mp_n, ecc_n, _ = merge_embryo(ap, mp, rp, M_sun, ecc, 1.0e4, live, 0.5)
    assert float(ecc_n[0]) == pytest.approx(0.0, abs=1e-15)
    assert float(mp_n[0]) == pytest.approx(mass_before, rel=1e-12)

    # Extreme mass ratio: the heavy body dominates both the eq. 16 orbit and
    # the eq. 17 eccentricity, so the merged values sit on its own.
    np.random.seed(1)
    ap, mp, rp, ecc, live = _pair(m1=1.0e6, m2=1.0, a1=1.0, a2=1.4, e1=0.02, e2=0.30)
    a_heavy, e_heavy = float(ap[0]), float(ecc[0])
    ap_n, mp_n, ecc_n, _ = merge_embryo(ap, mp, rp, M_sun, ecc, 1.0e4, live, 0.5)
    assert float(ap_n[0]) == pytest.approx(a_heavy, rel=1e-4)
    assert float(ecc_n[0]) == pytest.approx(e_heavy, rel=1e-4)
    assert float(mp_n[0]) == pytest.approx((1.0e6 + 1.0) * M_earth, rel=1e-12)

    lo = abs(2.0 * 0.05 - 1.0 * 0.08) / 3.0
    hi = (2.0 * 0.05 + 1.0 * 0.08) / 3.0
    merged_eccs = []
    for seed in range(10):
        np.random.seed(seed)
        ap, mp, rp, ecc, live = _pair(a2=1.05)
        v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
        _, _, ecc_n, _ = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5)
        assert lo * (1.0 - 1e-9) <= ecc_n[0] <= hi * (1.0 + 1e-9)
        merged_eccs.append(float(ecc_n[0]))
    # The alignment draw genuinely varies: the sweep is not degenerate.
    assert max(merged_eccs) - min(merged_eccs) > 0.0
