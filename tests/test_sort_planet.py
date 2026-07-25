"""Tests for ``src/morrigan/sort_planet.py``.

``sort_planet`` is the post-event cleanup: it drops dead bodies and
re-sorts every state array by semi-major axis. The contract is that all
nine arrays are filtered and permuted together, so a planet keeps its
own mass, radius, and identity through the shuffle; the tests build
systems where any de-synchronisation between arrays is visible.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.sort_planet import sort_planet

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _system(a, live):
    """Build a system whose per-planet values encode their identity.

    Each planet's mass, radius, eccentricity, and density are simple
    functions of its id, so after any filter-and-sort the cross-array
    consistency can be checked per body.
    """
    n = len(a)
    pid = np.arange(n)
    return dict(
        ap=np.asarray(a, dtype=float),
        Mp=1.0e24 * (1.0 + pid.astype(float)),
        ecc=0.01 * (1.0 + pid.astype(float)),
        Rp=1.0e6 * (1.0 + pid.astype(float)),
        live_status=np.asarray(live, dtype=bool),
        interact=np.ones(n, dtype=bool),
        densities=3000.0 + 100.0 * pid.astype(float),
        planet_id=pid,
    )


def test_dead_bodies_are_dropped_and_survivors_sorted_together():
    """Dead planets vanish and every array follows the same new order.

    Three of five bodies survive, listed out of semi-major-axis order,
    so the result must be length three, ascending in a, and each
    surviving id must still carry its own identity-encoded mass,
    radius, eccentricity, and density: any array left out of the
    shared permutation breaks at least one of these equalities.
    """
    s = _system([3.0, 1.0, 5.0, 2.0, 4.0], [True, False, True, True, False])
    ap, mp, ecc, rp, live, interact, rho, pid = sort_planet(**s)

    assert len(ap) == 3
    assert list(pid) == [3, 0, 2]  # ascending a: 2.0, 3.0, 5.0
    assert np.all(np.diff(ap) > 0.0)
    assert bool(np.all(live)) and bool(np.all(interact))

    for k, planet in enumerate(pid):
        assert mp[k] == pytest.approx(1.0e24 * (1.0 + planet), rel=1e-15)
        assert rp[k] == pytest.approx(1.0e6 * (1.0 + planet), rel=1e-15)
        assert ecc[k] == pytest.approx(0.01 * (1.0 + planet), rel=1e-15)
        assert rho[k] == pytest.approx(3000.0 + 100.0 * planet, rel=1e-15)


def test_edge_systems_survive_the_cleanup():
    """An all-dead system empties cleanly and a sorted one is unchanged.

    The all-dead edge must return eight empty arrays rather than raise,
    since the driver reaches this state when the last event kills the
    final interacting body. An already-sorted all-alive system must
    pass through as the identity, which guards against a spurious
    re-ordering of equal or ordered semi-major axes.
    """
    s = _system([1.0, 2.0], [False, False])
    out = sort_planet(**s)
    assert all(len(arr) == 0 for arr in out)

    s = _system([1.0, 2.0, 3.0], [True, True, True])
    ap, mp, ecc, rp, live, interact, rho, pid = sort_planet(**s)
    assert list(pid) == [0, 1, 2]
    np.testing.assert_array_equal(ap, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_array_equal(mp, 1.0e24 * np.array([1.0, 2.0, 3.0]))
