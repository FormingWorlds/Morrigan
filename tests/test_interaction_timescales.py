"""Tests for ``src/morrigan/interaction_timescales.py``.

The three timescales transcribe Kimura et al. (2025): the viscous
relaxation and collision timescales (their eqs. 6, 8-12, 23) and the
Petit et al. (2020) three-body crossing time with the eq. 28 stability
gate. The tests pin the two pair timescales at hand-checked reference
values, pin the analytic location of the stability boundary, and
exercise the guard branches: the overlapping-orbit zero, the
infinitely stable plateau, and the triplet-only contract of the
wrapper.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import M_earth, M_sun, au2m
from morrigan.helper_functions import hill_sphere, planet_radius
from morrigan.interaction_timescales import (
    interaction_wrapper,
    tau_col,
    tau_cross_petit,
    tau_vis,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0
_STABLE = 1e20  # the sentinel the module returns for an indefinitely stable system


def _pair():
    """A reference interacting pair: two Earths at 0.95 and 1.05 au."""
    ap = [0.95 * au2m, 1.05 * au2m]
    mp = [M_earth, M_earth]
    rp = [planet_radius(M_earth, _RHO)] * 2
    return ap, mp, rp


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_tau_vis_pins_its_closed_form():
    """The viscous relaxation timescale reproduces its transcription of
    Kimura et al. (2025) eqs. 6, 8-10 and 23 at a hand-checked point.

    Two Earths at 0.95 / 1.05 au with e = 0.05 give 1.74771e11 s
    (5538 yr), evaluated by hand through the chain rep_e -> random
    speed -> surface density -> stirring time. At this point the input
    eccentricities set rep_e: eq. 23 takes the larger of sum(e_cross)
    and sum(e), and for an equal-mass pair the two are equal to within
    one bit (0.09999999999999999 against 0.1), so the crossing
    eccentricities are not what the pin resolves. The unequal-mass case
    at the end of this test is what exercises eq. 6. The two guards are
    the transcription slips this chain invites: the factor-3 variant of
    the final expression (the paper-vs-Fortran discrepancy the source
    resolves in favour of the Fortran) and a dropped 2 pi in the surface
    number density; each moves the result by its own factor and lands
    far outside tolerance. The rep_e**4 scaling pins the eccentricity
    dependence exactly.
    """
    ap, mp, rp = _pair()
    t_vis = tau_vis(ap, mp, rp, M_sun, [0.05, 0.05])
    assert t_vis == pytest.approx(1.7477101e11, rel=1e-4)

    # Discrimination: the factor-3 variant of the final expression evaluates
    # to 5.2431303e11 s and the dropped-2pi variant to 2.7815670e10 s at this
    # reference point; both sit far outside the pin's tolerance band.
    assert abs(5.2431303e11 - t_vis) > 100 * 1e-4 * t_vis
    assert abs(2.7815670e10 - t_vis) > 100 * 1e-4 * t_vis

    # rep_e**4 scaling: doubling both eccentricities (above the geometric
    # crossing floor) scales the timescale by exactly 2**4 = 16.
    ratio = tau_vis(ap, mp, rp, M_sun, [0.1, 0.1]) / t_vis
    assert ratio == pytest.approx(16.0, rel=1e-9)

    # Unequal masses resolve the eq. 6 denominator: with a 1 and 10 Earth-mass
    # pair at near-zero input eccentricity, the crossing eccentricities set
    # rep_e and the shared-denominator form gives 6.4189361e9 s, hand-derived.
    # The variant that builds ecross_j with its own (sqrt(M_j) a_i swapped)
    # denominator gives 6.1122437e9 s, 4.8 percent away, so this pin
    # discriminates the two forms where the equal-mass point cannot.
    mp_uneq = [M_earth, 10 * M_earth]
    rp_uneq = [planet_radius(m, _RHO) for m in mp_uneq]
    t_uneq = tau_vis(ap, mp_uneq, rp_uneq, M_sun, [0.001, 0.001])
    assert t_uneq == pytest.approx(6.4189361e9, rel=1e-4)
    assert abs(6.1122437e9 - t_uneq) > 100 * 1e-4 * t_uneq

    # Edge: the timescale is positive and finite for a near-circular pair.
    t_circ = tau_vis(ap, mp, rp, M_sun, [1e-6, 1e-6])
    assert np.isfinite(t_circ) and t_circ > 0.0


@pytest.mark.physics_invariant
def test_tau_col_pins_and_returns_a_time_not_a_rate():
    """The collision timescale reproduces its hand-checked value and is
    returned as a time, not the collision rate it is built from.

    The same reference pair with e = 0.05 gives 1.97470e12 s. The
    sharpest guard is the inversion: the function assembles a collision
    rate (eq. 11) and returns its reciprocal, and forgetting the final
    inversion would return ~5e-13 s, twenty-five orders of magnitude
    away. Gravitational focusing is pinned directionally: enlarging both
    radii shortens the timescale.
    """
    ap, mp, rp = _pair()
    t_col = tau_col(ap, mp, rp, M_sun, [0.05, 0.05])
    assert t_col == pytest.approx(1.9746995e12, rel=1e-4)
    # Discrimination: the non-inverted rate would be smaller than 1 s.
    assert t_col > 1.0

    # Focusing and cross-section: bigger bodies collide sooner.
    rp_big = [2 * r for r in rp]
    assert tau_col(ap, mp, rp_big, M_sun, [0.05, 0.05]) < t_col

    # The eq. 6 denominator appears twice in the module, once per timescale,
    # so it needs pinning on both. An unequal 1 and 10 Earth-mass pair at
    # near-zero input eccentricity puts sum(e_cross) = 0.1027 far above
    # sum(e) = 0.002, so eq. 23 takes the crossing branch and the shared
    # denominator is what sets the answer: 2.4354420e11 s, hand-derived.
    # Building e_cross_j on its own swapped denominator instead gives
    # 2.3769219e11 s, 2.4 percent away and outside the tolerance band, which
    # the equal-mass reference point cannot distinguish at all.
    mp_uneq = [M_earth, 10 * M_earth]
    rp_uneq = [planet_radius(m, _RHO) for m in mp_uneq]
    t_uneq = tau_col(ap, mp_uneq, rp_uneq, M_sun, [0.001, 0.001])
    assert t_uneq == pytest.approx(2.4354420e11, rel=1e-4)
    assert abs(2.3769219e11 - t_uneq) > 100 * 1e-4 * t_uneq

    # Edge: a near-circular pair still returns a positive, finite time
    # (the geometric crossing eccentricity keeps rep_e above zero).
    t_circ = tau_col(ap, mp, rp, M_sun, [1e-6, 1e-6])
    assert np.isfinite(t_circ) and t_circ > 0.0


@pytest.mark.physics_invariant
def test_tau_cross_petit_orders_by_separation_and_guards_overlap():
    """The three-body crossing time lengthens with separation and its
    two guard branches return their sentinels.

    Equal-mass Earth triplets at 4, 6 and 8 mutual Hill radii must give
    strictly increasing, finite crossing times below the stable
    sentinel: packing controls instability. A triplet whose orbits
    already overlap (negative gap after the eccentricity terms) returns
    exactly 0.0, the immediate-instability guard, and a wide triplet
    beyond the resonance-overlap boundary returns the 1e20 stable
    sentinel.
    """
    a0 = 1.0 * au2m
    r_hill = hill_sphere(a0, 2 * M_earth, M_sun)
    times = []
    for sep in (4.0, 6.0, 8.0):
        a = [a0, a0 + sep * r_hill, a0 + 2 * sep * r_hill]
        times.append(tau_cross_petit(a, [M_earth] * 3, M_sun, [0.01] * 3, 3))
    assert all(np.isfinite(t) and 0.0 < t < _STABLE for t in times)
    assert times[0] < times[1] < times[2]

    # Guard: already-crossing orbits are instantaneously unstable.
    a_overlap = [1.0 * au2m, 1.02 * au2m, 1.04 * au2m]
    assert tau_cross_petit(a_overlap, [M_earth] * 3, M_sun, [0.3] * 3, 3) == pytest.approx(
        0.0, abs=0.0
    )

    # Guard: a wide triplet is stable indefinitely.
    a_wide = [0.5 * au2m, 1.0 * au2m, 2.0 * au2m]
    assert tau_cross_petit(a_wide, [M_earth] * 3, M_sun, [0.01] * 3, 3) == pytest.approx(
        _STABLE
    )


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_interaction_wrapper_pins_the_analytic_stability_boundary():
    """The two-planet stability gate flips at the separation eq. 28 of
    Kimura et al. (2025) predicts, and the wrapper is defined on
    triplets only.

    Analytical limit: the gate EJ = (5/8) sum(e^2)/h^2 - (3/8) x^2 + 4.5
    with x the pair separation in units of h * aM changes sign, for both
    eccentricities equal to 2h, at x_crit = sqrt((5 + 4.5) * 8 / 3) =
    5.0332. The test brackets that root: at x = 4.95 the wrapper falls
    through to the Petit crossing time (finite), at x = 5.10 it returns
    the stable sentinel. The guards show the bracket discriminates the
    coefficients: replacing 4.5 by 3.0 moves the root to 4.62, below the
    bracket, and dropping the 5/8 moves it to 5.77, above it. A bare
    unstable pair, without the third body, fails loudly with IndexError:
    the wrapper's contract is triplets.
    """
    a1 = 1.0 * au2m
    h = hill_sphere(a1, 2 * M_earth, M_sun) / a1
    ecc = 2.0 * h

    def gate(x):
        # place the pair at exactly x mutual-Hill-scaled separations of aM,
        # with a third body 100 Hill radii out, too remote to matter
        da = x * h * a1 / (1.0 - x * h / 2.0)
        return interaction_wrapper(
            [a1, a1 + da, a1 + da + 100.0 * h * a1],
            [M_earth] * 3,
            M_sun,
            [ecc, ecc, 0.0],
            2,
        )

    inside = gate(4.95)
    outside = gate(5.10)
    assert np.isfinite(inside) and 0.0 < inside < _STABLE
    assert outside == pytest.approx(_STABLE)

    # The bracket pins the eq. 28 coefficients: wrong-constant variants
    # put the root outside [4.95, 5.10].
    x_crit = np.sqrt((5.0 / 8.0 * 8.0 + 4.5) * 8.0 / 3.0)
    assert 4.95 < x_crit < 5.10
    x_wrong_const = np.sqrt((5.0 / 8.0 * 8.0 + 3.0) * 8.0 / 3.0)  # 4.5 -> 3.0
    x_wrong_coeff = np.sqrt((8.0 + 4.5) * 8.0 / 3.0)  # dropped 5/8
    assert x_wrong_const < 4.95
    assert x_wrong_coeff > 5.10

    # Contract: the wrapper needs a triplet once the pair is unstable.
    da = 3.0 * h * a1
    with pytest.raises(IndexError):
        interaction_wrapper([a1, a1 + da], [M_earth] * 2, M_sun, [ecc, ecc], 2)
