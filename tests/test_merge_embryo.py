"""Tests for ``src/morrigan/merge_embryo.py``.

``merge_embryo`` is the merger bookkeeping: it must close the mass
books (rock conserved, atmosphere split into retained and lost), place
the merged body on the eq. 16 orbit of Kimura et al. (2025), bound the
eq. 17 eccentricity, kill exactly the smaller body, and clamp the
Kegerreis loss fraction into [0, 1]. ``collision_velocity`` must floor
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


def _pair(m1=2.0, m2=1.0, a1=1.0, a2=1.2, e1=0.05, e2=0.08, f1=0.01, f2=0.02):
    """Return the mutable arrays merge_embryo works on, in model units."""
    ap = np.array([a1, a2]) * au2m
    mp = np.array([m1, m2]) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    ecc = np.array([e1, e2])
    live = np.array([True, True])
    f_atm = np.array([f1, f2])
    return ap, mp, rp, ecc, live, f_atm


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
    ap, mp, rp, _, _, _ = _pair()
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
    """A merger conserves rock, splits the atmosphere into retained plus
    lost, and lands on the eq. 16 orbit.

    Seeded with 7 (the pericentre-alignment draw). The merged
    semi-major axis is the orbital-energy-weighted form of eq. 16 of
    Kimura et al. (2025): a = (M1 + M2) / (M1/a1 + M2/a2), by hand
    1.05882 au for 2 M_e at 1.0 au absorbing 1 M_e at 1.2 au. Both the
    arithmetic mean (1.1 au) and the mass-weighted mean (1.0667 au),
    the two plausible wrong forms, sit far outside the pin. Rock
    closure is exact: total mass minus atmosphere is unchanged by the
    merger, and the surviving atmosphere fraction is stored over the
    new total mass, guarding the mass-into-fraction-array slip.
    """
    ap, mp, rp, ecc, live, f_atm = _pair()
    atm_before = float(mp[0] * f_atm[0] + mp[1] * f_atm[1])
    rock_before = float(mp.sum() - atm_before)

    np.random.seed(7)
    v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
    ap_n, mp_n, ecc_n, live_n, f_n, frac_lost = merge_embryo(
        ap, mp, rp, M_sun, ecc, v_c, live, 0.5, f_atm
    )

    a_expected = 3.0 * M_earth / (2.0 * M_earth / (1.0 * au2m) + 1.0 * M_earth / (1.2 * au2m))
    assert ap_n[0] == pytest.approx(a_expected, rel=1e-12)
    # Discrimination: arithmetic and mass-weighted means both fail the pin.
    assert abs(1.1 * au2m - a_expected) > 1e6
    assert abs((2.0 * 1.0 + 1.0 * 1.2) / 3.0 * au2m - a_expected) > 1e6

    # Rock is conserved through the atmosphere loss, exactly.
    rock_after = float(mp_n[0] * (1.0 - f_n[0]))
    assert rock_after == pytest.approx(rock_before, rel=1e-12)
    # Atmosphere books close: retained + lost = combined.
    atm_after = float(mp_n[0] * f_n[0])
    assert atm_after + frac_lost * atm_before == pytest.approx(atm_before, rel=1e-9)
    assert 0.0 <= frac_lost <= 1.0

    # The smaller body is dead with no atmosphere left on its slot.
    assert list(live_n) == [True, False]
    assert f_n[1] == pytest.approx(0.0, abs=0.0)
    # The stored fraction refers to the new total mass, not the old one.
    assert f_n[0] == pytest.approx(atm_after / mp_n[0], rel=1e-12)


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
    ap, mp, rp, ecc, live, f_atm = _pair(m1=1.0, m2=2.0)
    v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
    _, _, _, live_n, f_n, _ = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5, f_atm)
    assert list(live_n) == [False, True]
    assert f_n[0] == pytest.approx(0.0, abs=0.0)

    np.random.seed(11)
    ap, mp, rp, ecc, live, f_atm = _pair(m1=1.0, m2=1.0)
    v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
    _, _, _, live_tie, _, _ = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5, f_atm)
    assert list(live_tie) == [True, False]


@pytest.mark.physics_invariant
def test_extreme_impacts_clamp_the_loss_and_bound_the_eccentricity():
    """A far-above-escape impact clamps the loss fraction at one, and
    the merged eccentricity stays inside its vector-composition bounds.

    The raw Kegerreis law exceeds one for a head-on equal-mass impact
    at five times the contact floor, so the returned fraction must be
    exactly one and the merged body left with zero atmosphere: the
    call-site clamp contract. The eq. 17 eccentricity is a vector
    composition, so across seeds it must stay within
    [|M1 e1 - M2 e2|, M1 e1 + M2 e2] / (M1 + M2). The sweep uses a
    close pair (1.0 and 1.05 au) whose pericentre-alignment cosine is
    interior to (-1, 1), so the drawn alignment genuinely varies with
    the seed; seeds 0-9 must therefore produce a non-degenerate spread
    of merged eccentricities inside the bounds.
    """
    ap, mp, rp, ecc, live, f_atm = _pair(m1=1.0, m2=1.0, e1=0.05, e2=0.05)
    v_esc = np.sqrt(2.0 * G * (mp[0] + mp[1]) / (rp[0] + rp[1]))
    np.random.seed(0)
    _, mp_n, _, _, f_n, frac_lost = merge_embryo(
        ap, mp, rp, M_sun, ecc, 5.0 * v_esc, live, 0.0, f_atm
    )
    assert frac_lost == pytest.approx(1.0, rel=0.0, abs=0.0)
    assert f_n[0] == pytest.approx(0.0, abs=1e-15)
    # The lost atmosphere came off the merged mass.
    assert mp_n[0] < 2.0 * M_earth

    lo = abs(2.0 * 0.05 - 1.0 * 0.08) / 3.0
    hi = (2.0 * 0.05 + 1.0 * 0.08) / 3.0
    merged_eccs = []
    for seed in range(10):
        np.random.seed(seed)
        ap, mp, rp, ecc, live, f_atm = _pair(a2=1.05)
        v_c = collision_velocity(list(ap), list(mp), list(rp), M_sun, list(ecc))
        _, _, ecc_n, _, _, _ = merge_embryo(ap, mp, rp, M_sun, ecc, v_c, live, 0.5, f_atm)
        assert lo * (1.0 - 1e-9) <= ecc_n[0] <= hi * (1.0 + 1e-9)
        merged_eccs.append(float(ecc_n[0]))
    # The alignment draw genuinely varies: the sweep is not degenerate.
    assert max(merged_eccs) - min(merged_eccs) > 0.0
