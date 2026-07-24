"""Tests for ``src/morrigan/mass_loss.py``.

The file transcribes the Kegerreis et al. (2020) giant-impact
atmosphere-loss scaling, so the tests pin its closed form against a
hand-derived value, guard the two historically likely transcription
slips (the mass-ratio denominator and the outer exponent), and exercise
the algebraic limits: the grazing zero at b = 1 and the unclamped
head-on regime above one that the ``merge_embryo`` call site clamps.
The equal-density cross-check against the ZEPHYRUS implementation lives
in ``tests/test_mass_loss_properties.py``, the optional-dependency
companion of this file.
"""

from __future__ import annotations

import numpy as np
import pytest

from morrigan.constants import G, M_earth
from morrigan.helper_functions import planet_radius
from morrigan.mass_loss import mass_loss

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0  # bulk density shared by both bodies [kg m-3]


def _pair(q):
    """Return (M_t, M_i, R_t, R_i, v_esc) for an impactor of mass ratio q."""
    m_t, m_i = M_earth, q * M_earth
    r_t, r_i = planet_radius(m_t, _RHO), planet_radius(m_i, _RHO)
    v_esc = np.sqrt(2.0 * G * (m_t + m_i) / (r_t + r_i))
    return m_t, m_i, r_t, r_i, v_esc


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_pins_the_kegerreis_scaling_closed_form():
    """The loss fraction reproduces the published scaling, eq. (1) of
    Kegerreis et al. (2020), at a hand-derived reference point.

    Reference point: a 0.1 mass-ratio impactor at 1.5 v_esc, equal bulk
    densities, b = 0.7. Evaluating X = 0.64 [ (v_c/v_esc)^2
    (M_i/M_tot)^(1/2) (rho_i/rho_t)^(1/2) f_V(b) ]^0.65 by hand gives
    0.147467. The two discrimination guards are the transcription slips
    this formula invites: dividing by M_t instead of M_tot gives
    0.152106, and an outer exponent of 0.5 instead of 0.65 gives
    0.206922; both sit far outside the pin's tolerance.
    """
    m_t, m_i, r_t, r_i, v_esc = _pair(0.1)
    x = mass_loss(1.5 * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, 0.7)
    assert x == pytest.approx(0.1474666, rel=1e-5)

    # Discrimination guards: the nearest wrong formulas land >100 tolerances away.
    assert abs(0.1521060 - x) > 100 * 1e-5 * x  # M_i/M_t denominator slip
    assert abs(0.2069216 - x) > 100 * 1e-5 * x  # 0.5 outer exponent slip

    # The loss fraction weakens monotonically toward grazing geometry.
    grid = [mass_loss(1.5 * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, b)
            for b in (0.0, 0.3, 0.7, 0.99)]
    assert all(a > b for a, b in zip(grid, grid[1:]))


@pytest.mark.physics_invariant
def test_grazing_limit_is_an_algebraic_zero():
    """A fully grazing impact, b = 1, removes no atmosphere at all.

    The interacting-volume factor carries (1 - b)^2, so b = 1 is an
    exact algebraic zero of the law, not a small number; this is the
    limit-input contract of the geometry. Just inside the limit the loss
    is already tiny but strictly positive, so the zero is approached
    continuously rather than by a jump.
    """
    m_t, m_i, r_t, r_i, v_esc = _pair(0.5)
    assert mass_loss(2.0 * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, 1.0) == pytest.approx(
        0.0, abs=1e-15
    )
    near = mass_loss(2.0 * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, 1.0 - 1e-6)
    assert 0.0 < near < 1e-6


@pytest.mark.physics_invariant
def test_head_on_energetic_impacts_exceed_one_unclamped():
    """The raw law exceeds one for energetic head-on mergers and the
    clamp is the caller's responsibility.

    An equal-mass head-on collision at twice the mutual escape speed
    evaluates to 1.25801 by hand, above the physical ceiling: the fitted
    power law is only meaningful on [0, 1] and ``merge_embryo`` clamps
    it at the call site. Pinning the raw value here documents that this
    function does not clamp, so a silent clamp added inside it would be
    caught as a behaviour change.
    """
    m = M_earth
    r = planet_radius(m, _RHO)
    v_esc = np.sqrt(2.0 * G * 2 * m / (2 * r))
    x = mass_loss(2.0 * v_esc, m, m, _RHO, _RHO, r, r, 0.0)
    assert x == pytest.approx(1.2580104, rel=1e-5)
    assert x > 1.0

    # At the mutual escape speed itself, the slowest possible contact,
    # the same geometry stays below one.
    x_slow = mass_loss(v_esc, m, m, _RHO, _RHO, r, r, 0.0)
    assert 0.0 < x_slow < 1.0


@pytest.mark.physics_invariant
def test_loss_grows_with_speed_and_impactor_share():
    """More speed and a larger impactor share both strip more atmosphere.

    Monotonicity in v_c at fixed geometry and in the impactor mass ratio
    at fixed speed-to-escape ratio are the two directional dependences
    of the published fit. The edge case is the vanishing impactor,
    where the loss tends to zero from above.
    """
    m_t, m_i, r_t, r_i, v_esc = _pair(0.3)
    speeds = [mass_loss(f * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, 0.3)
              for f in (1.0, 1.5, 2.0, 3.0)]
    assert all(a < b for a, b in zip(speeds, speeds[1:]))

    ratios = []
    for q in (1e-4, 0.01, 0.1, 1.0):
        m_t, m_i, r_t, r_i, v_esc = _pair(q)
        ratios.append(mass_loss(1.5 * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, 0.3))
    assert all(a < b for a, b in zip(ratios, ratios[1:]))
    # Edge: a vanishing impactor strips little, from above; the tail is the
    # slow q**(0.5*0.65) power of the mass ratio, so q = 1e-4 still leaves
    # a ~2 percent loss rather than an exponential cutoff.
    assert 0.0 < ratios[0] < 0.05
