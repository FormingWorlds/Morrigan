"""Tests for ``src/morrigan/crossing_pair.py``.

``crossing_pair`` schedules the next orbit-crossing event and names the
pair that will cross. The tests pin the two-body branches (a stable pair
is deferred to the 1e20 sentinel, an unstable pair is scheduled at
1.5 t), and pin the eq. 4 pair-selection geometry of Kimura et al.
(2025): the interacting pair is the one with the smaller
perihelion-to-aphelion gap, which eccentricity can flip away from the
plain semi-major-axis gap ordering.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import M_earth, M_sun, au2m
from morrigan.crossing_pair import crossing_pair
from morrigan.helper_functions import planet_radius

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0
_STABLE = 1e20


def _call(a_au, ecc, t, n=None, masses=None):
    """Call crossing_pair on a system given in au, with zero secular modes.

    ``masses`` is in Earth masses and defaults to an equal-mass system.
    Pass unequal masses to resolve the eq. 6 denominator, which is
    symmetric under any mass-index error when every body weighs the same.
    """
    n = n if n is not None else len(a_au)
    ap = np.asarray(a_au, dtype=float) * au2m
    mp = np.full(n, M_earth) if masses is None else np.asarray(masses, dtype=float) * M_earth
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    ecc = np.asarray(ecc, dtype=float)
    ecc_vec = np.zeros((n, n))
    g = np.zeros(n)
    beta = np.zeros(n)
    interact = np.ones(n, dtype=bool)
    return crossing_pair(ap, mp, rp, M_sun, ecc, ecc_vec, g, beta, interact, n, t, t)


@pytest.mark.physics_invariant
def test_two_body_branches_defer_or_schedule():
    """A stable two-planet system is deferred indefinitely; an unstable
    one is scheduled promptly at one and a half times the current time.

    The stable branch (wide, near-circular pair) must return the 1e20
    sentinel through the empty triplet loop, and the unstable branch
    (overlapping pair) returns exactly 1.5 t, the prompt-event schedule.
    The edge case is t = 0, where the prompt schedule degenerates to an
    immediate event at time zero.
    """
    icross, t_event = _call([1.0, 1.5], [0.001, 0.001], t=8e10)
    assert icross == 0
    assert t_event == pytest.approx(_STABLE)

    icross, t_event = _call([1.0, 1.02], [0.05, 0.05], t=8e10)
    assert icross == 0
    assert t_event == pytest.approx(1.5 * 8e10, rel=1e-12)
    # Discrimination: the two branches are far apart, not near-degenerate.
    assert t_event < _STABLE / 1e6

    # Edge: at t = 0 the unstable schedule degenerates to an immediate event.
    _, t_zero = _call([1.0, 1.02], [0.05, 0.05], t=0.0)
    assert t_zero == pytest.approx(0.0, abs=0.0)


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_pair_selection_follows_the_perihelion_gap_geometry():
    """The interacting pair is chosen by the closest physical approach,
    eq. 4 of Kimura et al. (2025), not by the plain semi-major-axis gap.

    Analytical geometry: the gap between neighbours is the outer body's
    perihelion minus the inner body's aphelion, d = (1 - e_out) a_out -
    (1 + e_in) a_in. The triplet (1.0, 1.06, 1.13) au is packed tightly
    enough (about five mutual Hill radii) to be genuinely unstable, so
    the scheduled pair is the eq. 4 choice and not the untouched
    initialiser slot. Circular, the plain gaps order inner < outer
    (0.06 < 0.07 au) and the inner pair is selected; giving the
    outermost body e = 0.05 collapses the outer physical gap to
    0.95 * 1.13 - 1.06 = 0.0135 au, hand-derived, and must flip the
    selection to the outer pair. A selection built on plain
    semi-major-axis differences could never respond to that
    eccentricity, which is the discrimination this flip provides.
    """
    icross_circular, t_circ = _call([1.0, 1.06, 1.13], [0.0, 0.0, 0.0], t=1e9)
    assert icross_circular == 0
    assert np.isfinite(t_circ) and t_circ < _STABLE

    icross_eccentric, t_ecc = _call([1.0, 1.06, 1.13], [0.0, 0.0, 0.05], t=1e9)
    assert icross_eccentric == 1

    # Both schedules run strictly forward from the current time.
    assert t_circ > 1e9 and t_ecc > 1e9

    # Edge: with the middle body eccentric instead, both physical gaps
    # shrink together and the inner pair keeps the smaller one
    # (0.95 * 1.06 - 1.0 = 0.007 < 1.13 - 1.05 * 1.06 = 0.017), hand-derived.
    icross_middle, _ = _call([1.0, 1.06, 1.13], [0.0, 0.05, 0.0], t=1e9)
    assert icross_middle == 0


@pytest.mark.physics_invariant
def test_event_times_run_forward_for_a_packed_system():
    """A tightly packed multi-planet system schedules a finite event in
    the future, never in the past.

    Five equal-mass planets spaced a few Hill radii apart are deep in
    the unstable regime, so the scheduled event must be finite, sit
    strictly after the current time, and name an interior pair index.
    The guard input is the same system widened tenfold, which must push
    the schedule out by many orders of magnitude: the scheduler responds
    to packing, not just to planet count.
    """
    a_packed = [1.0, 1.03, 1.06, 1.09, 1.12]
    icross, t_event = _call(a_packed, [0.01] * 5, t=1e8)
    assert 0 <= icross <= 3
    assert t_event > 1e8
    assert np.isfinite(t_event)

    a_wide = [1.0, 1.3, 1.6, 1.9, 2.2]
    _, t_wide = _call(a_wide, [0.01] * 5, t=1e8)
    assert t_wide > 100.0 * t_event


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_unequal_masses_resolve_the_crossing_eccentricity_denominator():
    """The scheduled event time carries the shared eq. 6 denominator.

    Eq. 6 gives both bodies of a pair the same denominator,
    sqrt(M_j) a_i + sqrt(M_i) a_j, pairing each mass with the OTHER
    body's semi-major axis. Building each body's crossing eccentricity
    on its own denominator instead is the transcription slip the form
    invites, and it is invisible in an equal-mass system, where the two
    expressions are algebraically identical. A 1, 4 and 0.5 Earth-mass
    triplet at 1.00, 1.06 and 1.13 au resolves them.

    The eccentricity matters here: eq. 23 takes the larger of the input
    and crossing eccentricities, and at e = 0.02 the correctly built
    crossing value wins while the swapped one is clamped away by the
    input, so the two forms schedule measurably different events. The
    pinned separation is 2.2e-4 relative, which is why the tolerance is
    tight; at lower eccentricities the two forms coincide and the pin
    would not discriminate at all.
    """
    icross, t_event = _call([1.0, 1.06, 1.13], [0.02, 0.02, 0.02], 0.0,
                            masses=[1.0, 4.0, 0.5])

    assert icross == 0
    assert t_event == pytest.approx(2.5199888664e11, rel=1e-6)
    # Discrimination: the per-body denominator schedules 2.5194397394e11 s,
    # 2.2e-4 away, which is 200 times the tolerance above.
    assert abs(2.5194397394e11 - t_event) > 100 * 1e-6 * t_event
    # The event is scheduled, not deferred: the stable sentinel would pass
    # the bounds above by accident.
    assert t_event < _STABLE
