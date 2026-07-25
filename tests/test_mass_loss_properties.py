"""Property-based and cross-implementation tests for
``src/morrigan/mass_loss.py``.

This is the optional-dependency companion of ``tests/test_mass_loss.py``:
it needs ``hypothesis`` for the domain sweeps and ``zephyrus`` for the
cross-implementation check, so both are import-or-skipped at module top
and the mandatory closed-form pins stay in the main file, which always
runs. The cross-check pins the equal-density identity between Morrigan's
interacting-volume form and the ZEPHYRUS cap-geometry form of the
Kegerreis et al. (2020) law: the b-dependence ratios agree to machine
precision, and the absolute values agree once the two packages' G
constants (6.67e-11 here, CODATA in ZEPHYRUS) are allowed for.
"""

from __future__ import annotations

import numpy as np
import pytest

hypothesis = pytest.importorskip('hypothesis')
zephyrus_collision = pytest.importorskip('zephyrus.collision')

from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

from morrigan.constants import G, M_earth  # noqa: E402
from morrigan.helper_functions import planet_radius  # noqa: E402
from morrigan.mass_loss import mass_loss  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

_RHO = 3000.0


def _geometry(q, m_t):
    """Return (m_i, r_t, r_i, v_esc) for a mass-ratio q impactor on m_t."""
    m_i = q * m_t
    r_t, r_i = planet_radius(m_t, _RHO), planet_radius(m_i, _RHO)
    v_esc = np.sqrt(2.0 * G * (m_t + m_i) / (r_t + r_i))
    return m_i, r_t, r_i, v_esc


@pytest.mark.physics_invariant
@settings(derandomize=True, max_examples=60, deadline=None)
@given(
    q=st.floats(min_value=1e-3, max_value=1.0),
    m_t=st.floats(min_value=1e22, max_value=1e26),
    b=st.floats(min_value=0.0, max_value=1.0),
    v_lo=st.floats(min_value=1.0, max_value=2.0),
    v_step=st.floats(min_value=0.05, max_value=1.0),
)
def test_loss_is_nonnegative_finite_and_monotone_in_speed(q, m_t, b, v_lo, v_step):
    """Across the fitted domain the raw loss is nonnegative, finite, and
    grows with collision speed.

    The impactor mass is co-drawn as a ratio q <= 1 of the target so the
    sweep never leaves the impactor-lighter-than-target regime the fit
    was built for. The b = 1 edge is inside the drawn range, where the
    loss must be exactly zero, which is the guard-limit behaviour of the
    geometry factor. Derandomized, so the sweep is reproducible.
    """
    m_i, r_t, r_i, v_esc = _geometry(q, m_t)
    x_lo = mass_loss(v_lo * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, b)
    x_hi = mass_loss((v_lo + v_step) * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, b)

    assert np.isfinite(x_lo) and np.isfinite(x_hi)
    assert x_lo >= 0.0
    if b >= 1.0:  # the grazing boundary of the drawn range
        assert x_lo == pytest.approx(0.0, abs=1e-15)
    else:
        assert x_hi > x_lo  # faster contact always strips more


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
@settings(derandomize=True, max_examples=60, deadline=None)
@given(
    q=st.floats(min_value=0.02, max_value=1.0),
    b=st.floats(min_value=0.0, max_value=0.95),
    v_fac=st.floats(min_value=1.0, max_value=1.6),
)
def test_equal_density_agreement_with_zephyrus(q, b, v_fac):
    """At equal bulk densities Morrigan's form matches the ZEPHYRUS
    implementation of the same published law.

    Cross-implementation cross-check: ``zephyrus.collision.mass_loss``
    implements the Kegerreis et al. (2020) appendix-B cap geometry,
    which at equal densities reduces algebraically to Morrigan's
    interacting-volume form. The b-dependence ratio X(b)/X(0) cancels
    each package's gravitational constant and must agree to 1e-9; the
    absolute values agree to 1e-3, the size of the residual left by
    Morrigan's G = 6.67e-11 against ZEPHYRUS's CODATA G, a factor
    (G_m/G_z)**0.65 = 0.99958 on X. Draws are derandomized and clamped
    below the X = 1 ceiling so the ZEPHYRUS-side clamp never engages.
    """
    m_t = M_earth
    m_i, r_t, r_i, v_esc = _geometry(q, m_t)

    x_m0 = mass_loss(v_fac * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, 0.0)
    hypothesis.assume(x_m0 < 0.95)  # keep clear of the [0, 1] clamp in ZEPHYRUS

    x_mb = mass_loss(v_fac * v_esc, m_i, m_t, _RHO, _RHO, r_i, r_t, b)
    x_z0 = zephyrus_collision.mass_loss(
        v_c=v_fac * v_esc, M_i=m_i, M_t=m_t, rho_i=_RHO, rho_t=_RHO,
        R_i=r_i, R_t=r_t, b=0.0,
    )
    x_zb = zephyrus_collision.mass_loss(
        v_c=v_fac * v_esc, M_i=m_i, M_t=m_t, rho_i=_RHO, rho_t=_RHO,
        R_i=r_i, R_t=r_t, b=b,
    )

    # Geometry identity, G-free: the two b-dependences are the same law.
    assert x_mb / x_m0 == pytest.approx(x_zb / x_z0, rel=1e-9)
    # Absolute agreement, bounded by the G-constant residual only.
    assert x_mb == pytest.approx(x_zb, rel=1.5e-3)
