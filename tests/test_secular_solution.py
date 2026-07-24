"""Tests for ``src/morrigan/secular_solution.py``.

The module solves the classical Laplace-Lagrange secular problem for the
system's eccentricity evolution. The defining property of that
eigen-decomposition is that the returned modes reconstruct the input
eccentricities exactly at the reference epoch, which is the analytical
anchor here (the linear secular solution of celestial mechanics, e.g.
Murray & Dermott 1999, ch. 7). The tests also cover the guard branch
that decouples radially overlapping orbits, determinism under the
global numpy seed, and the amplitude envelope bound.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import M_earth, M_sun, au2m
from morrigan.helper_functions import planet_radius
from morrigan.secular_solution import secular_solution

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0


def _system():
    """A three-planet system with distinct masses and eccentricities."""
    ap = np.array([0.8, 1.0, 1.3]) * au2m
    mp = np.array([1.0, 2.0, 0.5]) * M_earth
    ecc = np.array([0.03, 0.05, 0.02])
    rp = np.array([planet_radius(m, _RHO) for m in mp])
    return ap, mp, ecc, rp


def _reconstruct(ecc_vec, beta):
    """Rebuild each planet's epoch eccentricity from the returned modes."""
    n = ecc_vec.shape[0]
    h0 = np.array([np.sum(ecc_vec[i, :] * np.sin(beta)) for i in range(n)])
    k0 = np.array([np.sum(ecc_vec[i, :] * np.cos(beta)) for i in range(n)])
    return np.sqrt(h0**2 + k0**2)


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_recovers_the_initial_eccentricities_at_epoch():
    """The secular modes reconstruct the input eccentricities exactly at
    the reference epoch.

    Analytical limit: in the Laplace-Lagrange solution the integration
    constants are fitted to the epoch state, so summing the scaled
    eigenvectors at phase beta must return the input eccentricities to
    machine precision (the classical linear secular solution, Murray &
    Dermott 1999, ch. 7). Seeded with 42; the pericentre angles are
    random but the recovery must hold for every draw, so the assertion
    discriminates any error in the eigen-solve, the constant fit, or
    the amplitude scaling at once. The edge case is a near-circular
    system, where the recovery must survive eccentricities of 1e-8.
    """
    ap, mp, ecc, rp = _system()

    np.random.seed(42)
    ecc_vec, g, beta = secular_solution(ap, mp, ecc, rp, M_sun, 3)
    np.testing.assert_allclose(_reconstruct(ecc_vec, beta), ecc, rtol=1e-9)
    assert np.all(np.isfinite(g))
    # A coupled, non-overlapping system has at least one nonzero frequency.
    assert np.any(np.abs(g) > 0.0)

    # Edge: near-circular input is recovered at the same precision.
    tiny = np.full(3, 1e-8)
    np.random.seed(42)
    ecc_vec_t, _, beta_t = secular_solution(ap, mp, tiny, rp, M_sun, 3)
    np.testing.assert_allclose(_reconstruct(ecc_vec_t, beta_t), tiny, rtol=1e-6)


@pytest.mark.physics_invariant
def test_overlapping_orbits_are_decoupled():
    """Radially overlapping neighbours are excluded from the secular
    coupling, leaving the interaction matrix empty.

    Guard branch: when one body's perihelion dips inside its
    neighbour's aphelion the linear expansion is invalid, so the matrix
    entries are skipped. For a two-planet system whose orbits fully
    overlap (a = 1.0 and 1.02 au, e = 0.3 both) every coupling is
    skipped, both eigenfrequencies must be exactly zero, and the epoch
    recovery still holds through the then-trivial eigenbasis. Seeded
    with 1.
    """
    ap = np.array([1.0, 1.02]) * au2m
    mp = np.array([1.0, 1.0]) * M_earth
    ecc = np.array([0.3, 0.3])
    rp = np.array([planet_radius(M_earth, _RHO)] * 2)

    np.random.seed(1)
    ecc_vec, g, beta = secular_solution(ap, mp, ecc, rp, M_sun, 2)
    np.testing.assert_allclose(g, 0.0, atol=0.0)
    np.testing.assert_allclose(_reconstruct(ecc_vec, beta), ecc, rtol=1e-9)

    # Contrast: separating the same pair restores a nonzero frequency,
    # so the zero above discriminates the guard, not a degenerate solve.
    ap_wide = np.array([1.0, 1.4]) * au2m
    np.random.seed(1)
    _, g_wide, _ = secular_solution(ap_wide, mp, np.array([0.02, 0.02]), rp, M_sun, 2)
    assert np.any(np.abs(g_wide) > 0.0)


@pytest.mark.physics_invariant
def test_seeded_draws_make_the_solution_deterministic():
    """The same global seed reproduces the same secular solution.

    The pericentre angles are drawn from numpy's global state, so two
    calls under the same seed (7 here) must agree bit for bit in all
    three outputs, and a different seed (8) must change the phase
    angles: the determinism contract of the model's Monte Carlo layer.
    The amplitude envelope is also bounded here: the sum of mode
    amplitudes per planet can never fall below its epoch eccentricity,
    since the modes must at least reach the epoch state.
    """
    ap, mp, ecc, rp = _system()

    np.random.seed(7)
    v1, g1, b1 = secular_solution(ap, mp, ecc, rp, M_sun, 3)
    np.random.seed(7)
    v2, g2, b2 = secular_solution(ap, mp, ecc, rp, M_sun, 3)
    np.testing.assert_array_equal(v1, v2)
    np.testing.assert_array_equal(g1, g2)
    np.testing.assert_array_equal(b1, b2)

    np.random.seed(8)
    _, _, b3 = secular_solution(ap, mp, ecc, rp, M_sun, 3)
    assert np.any(np.abs(b3 - b1) > 0.0)

    # Envelope bound: total mode amplitude per planet covers the epoch value.
    envelope = np.sum(np.abs(v1), axis=1)
    assert np.all(envelope >= ecc * (1.0 - 1e-12))
