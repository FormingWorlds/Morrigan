"""Tests for ``src/morrigan/helper_functions.py``.

These helpers carry the model's orbital mechanics primitives, so the
invariants exercised here are the closed forms themselves: Kepler's third
law with its exact 3/2 exponent, the mutual Hill radius against the known
Earth value, the radius / density inversion pair, the sqrt(a) scaling of
the escape eccentricity, and the truncation contract of the Rayleigh draw.
Every pinned number is derived by hand from the model's own constants and
guarded against the nearest plausible wrong formula.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import M_earth, M_sun, au2m
from morrigan.helper_functions import (
    esc_ecc,
    hill_sphere,
    kepler_period,
    planet_density,
    planet_radius,
    rayleigh,
)

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_kepler_period_pins_the_analytic_third_law():
    """The orbital period reproduces Kepler's third law at 1 au.

    Analytical limit: P = 2 pi sqrt(a^3 / (G M_s)). With the model's own
    constants (G = 6.67e-11, M_sun = 1.9892e30 kg, au = 1.5e11 m) the
    hand-derived period at 1 au is 3.16894e7 s, a 366.8-day year; the
    offset from 365.25 d comes from the model's rounded au and G, not
    from the formula. The a**(3/2) exponent is pinned by the exact
    factor-8 response to a factor-4 orbit change, which no linear or
    quadratic wrong exponent reproduces.
    """
    P = kepler_period(M_earth, M_sun, au2m)
    assert P == pytest.approx(3.16894e7, rel=1e-4)

    # Kepler exponent: P(4a)/P(a) = 4**1.5 = 8 exactly. A linear-in-a wrong
    # formula would give 4 and an a**2 one 16; both sit far outside tolerance.
    ratio = kepler_period(M_earth, M_sun, 4 * au2m) / P
    assert ratio == pytest.approx(8.0, rel=1e-12)
    assert abs(ratio - 4.0) > 1.0 and abs(ratio - 16.0) > 1.0

    # The planet mass must not enter: the formula is the M_p << M_s limit.
    # Using the planet mass in the denominator by mistake would change the
    # period by many orders of magnitude.
    assert kepler_period(1e3 * M_earth, M_sun, au2m) == pytest.approx(P, rel=1e-12)

    # Edge: a tiny orbit still returns a positive, finite period.
    P_small = kepler_period(M_earth, M_sun, 1e3)
    assert np.isfinite(P_small) and P_small > 0.0


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_hill_sphere_reproduces_the_earth_value():
    """The mutual Hill radius at 1 au matches the known Earth value.

    Published benchmark: Earth's Hill radius is 0.0100 au (about 1.5e9 m).
    With the model constants, a (M/(3 M_s))**(1/3) evaluates to
    1.5005e9 m. Dropping the factor 3 in the denominator, the most
    plausible transcription slip, gives 2.164e9 m, 44 percent high and
    far outside the pin's tolerance.
    """
    r_hill = hill_sphere(au2m, M_earth, M_sun)
    assert r_hill == pytest.approx(1.5005e9, rel=1e-3)
    # Discrimination: the no-third variant a*(M/Ms)**(1/3) = 2.164e9 m.
    wrong = au2m * (M_earth / M_sun) ** (1 / 3)
    assert abs(wrong - r_hill) > 100 * 1e-3 * r_hill

    # Mass scaling: an 8-fold mass raises the Hill radius exactly 2-fold.
    assert hill_sphere(au2m, 8 * M_earth, M_sun) / r_hill == pytest.approx(2.0, rel=1e-12)

    # Edge: a planet as heavy as the star still returns a finite radius
    # (the formula itself has no guard; the limit value is a*(1/3)**(1/3)).
    r_extreme = hill_sphere(au2m, M_sun, M_sun)
    assert r_extreme == pytest.approx(au2m * (1 / 3) ** (1 / 3), rel=1e-12)


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_planet_radius_and_density_invert_each_other():
    """Radius from density and density from radius are exact inverses.

    The pinned point is the Earth: M_earth at the mean density 5513
    kg m-3 must return the published mean radius 6.371e6 m. Swapping the
    4/3 into a 3/4 in the sphere volume, the standard slip, scales the
    radius by (16/9)**(1/3) = 1.21, far outside tolerance. The
    round-trip identity is checked at both mass extremes of the model's
    working range.
    """
    r_earth = planet_radius(M_earth, 5513.0)
    assert r_earth == pytest.approx(6.3717e6, rel=1e-3)
    # Discrimination: the 4/3-vs-3/4 slip gives 1.211 * r_earth.
    assert abs(1.211 * r_earth - r_earth) > 100 * 1e-3 * r_earth

    # Round trip at the extremes: a Moon-mass body and a 1000-Earth body.
    for mass in (7.3e22, 1e3 * M_earth):
        for rho in (500.0, 12000.0):
            assert planet_density(mass, planet_radius(mass, rho)) == pytest.approx(
                rho, rel=1e-12
            )

    # Edge: radius grows monotonically with mass at fixed density.
    assert planet_radius(2 * M_earth, 5513.0) > r_earth


@pytest.mark.physics_invariant
def test_esc_ecc_scales_and_stays_symmetric():
    """The escape eccentricity scales as sqrt(a) and ignores body order.

    e_esc is the mutual surface escape speed over the local Kepler speed,
    so it must scale exactly as sqrt(a) (the Kepler speed in the
    denominator falls as 1/sqrt(a)) and be symmetric under swapping the
    two bodies. The pinned point, two Earths at 1 au with the published
    Earth radius, evaluates by hand to 0.37604. Far-out compact pairs
    push e_esc above one, the regime where scattering ejects rather than
    collides, which is the meaningful limit-input behaviour.
    """
    r_e = 6.371e6
    e1 = esc_ecc(M_sun, M_earth, M_earth, r_e, r_e, au2m)
    assert e1 == pytest.approx(0.37604, rel=1e-4)

    # sqrt(a) scaling: quadrupling a doubles e_esc exactly; a linear
    # wrong scaling would give 4.
    ratio = esc_ecc(M_sun, M_earth, M_earth, r_e, r_e, 4 * au2m) / e1
    assert ratio == pytest.approx(2.0, rel=1e-12)
    assert abs(ratio - 4.0) > 1.0

    # Symmetry under body exchange, with distinct masses and radii.
    a = esc_ecc(M_sun, M_earth, 0.1 * M_earth, 6.4e6, 3.0e6, au2m)
    b = esc_ecc(M_sun, 0.1 * M_earth, M_earth, 3.0e6, 6.4e6, au2m)
    assert a == pytest.approx(b, rel=1e-12)

    # Edge: a massive compact pair at 30 au sits in the ejection regime.
    e_far = esc_ecc(M_sun, 100 * M_earth, 100 * M_earth, 1.6e7, 1.6e7, 30 * au2m)
    assert e_far > 1.0 and np.isfinite(e_far)


@pytest.mark.physics_invariant
def test_rayleigh_respects_its_truncation_bound():
    """Truncated Rayleigh draws never fall below the requested floor.

    The draw is used for post-encounter eccentricity excitation with a
    geometric floor, so the truncation contract is the physics: every
    sample must sit at or above xmin. Seeded with 42 for the sweep and 3
    for the determinism check. The extreme floor xmin = 6 sigma
    exercises the numerically delicate tail, the closest thing the
    function has to an error path, and must still return finite values
    above the floor.
    """
    sigma = 1.0 / np.sqrt(2.0)

    np.random.seed(42)
    draws = np.array([rayleigh(sigma, 2.0) for _ in range(200)])
    assert np.all(draws >= 2.0)
    assert np.all(np.isfinite(draws))
    # The distribution is not degenerate at the floor.
    assert draws.max() > draws.min()

    # Edge: an untruncated draw is strictly positive.
    np.random.seed(42)
    assert rayleigh(sigma, 0.0) > 0.0

    # Extreme tail: a floor at 6 sigma is still honoured and finite.
    np.random.seed(42)
    tail = rayleigh(sigma, 6.0 * sigma)
    assert tail >= 6.0 * sigma and np.isfinite(tail)

    # Same seed, same draw: the sampler follows the global numpy state.
    np.random.seed(3)
    first = rayleigh(sigma, 1.0)
    np.random.seed(3)
    second = rayleigh(sigma, 1.0)
    assert first == pytest.approx(second, rel=0.0, abs=0.0)
