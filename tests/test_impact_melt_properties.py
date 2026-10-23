"""Property tests for ``src/morrigan/impact_melt.py`` over each argument's domain.

Kept apart from ``tests/test_impact_melt.py`` so that ``pytest.importorskip('hypothesis')``
skips only these tests where hypothesis is missing, not the reference cross-check.
"""

from __future__ import annotations

import math

import pytest

hypothesis = pytest.importorskip('hypothesis')
from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

from morrigan.constants import G, M_earth  # noqa: E402
from morrigan.impact_melt import melt_fraction  # noqa: E402

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


def _case(theta, gamma=0.1, vel=1.0, m_total=M_earth, density=5500.0):
    """Inputs for a merger of total mass ``m_total`` at one shared density."""
    m_i, m_t = gamma * m_total, (1.0 - gamma) * m_total
    r_t, r_i = ((3.0 * m / (4.0 * math.pi * density)) ** (1.0 / 3.0) for m in (m_t, m_i))
    v_esc = math.sqrt(2.0 * G * m_total / (r_t + r_i))
    return m_t, m_i, r_t, r_i, vel * v_esc, math.sin(math.radians(theta))


_NAMES = 'M_target M_impactor R_target R_impactor v_impact impact_parameter'.split()
_ARGS = dict(zip(_NAMES, _case(45.0)))


_BAD = st.one_of(
    st.just(0.0),
    st.floats(max_value=-1e-300, allow_nan=False),
    st.just(float('nan')),
    st.just(float('inf')),
    st.just(float('-inf')),
    st.booleans(),
)


@settings(derandomize=True, max_examples=200, deadline=None)
@given(
    which=st.sampled_from(['M_target', 'M_impactor', 'R_target', 'R_impactor', 'v_impact']),
    bad=_BAD,
)
def test_property_bad_positive_argument_is_refused(which, bad):
    """Zero, negative, NaN, inf and booleans are refused for every positive argument."""
    with pytest.raises(ValueError, match=which):
        melt_fraction(**(_ARGS | {which: bad}))


@settings(derandomize=True, max_examples=300, deadline=None)
@given(
    m_total=st.floats(0.01 * M_earth, 20.0 * M_earth),
    gamma=st.floats(1e-4, 0.5),
    vel=st.floats(1.0, 3.0),
    b=st.floats(0.0, 1.0),
)
def test_property_valid_merger_gives_a_bounded_result(m_total, gamma, vel, b):
    """Over masses, mass ratios, speeds and impact parameters well past the fit range, the
    melt fraction is in [0, 1], its sigma and dIE are finite and positive, and the flags
    are booleans."""
    m_t, m_i, r_t, r_i, v, _ = _case(0.0, gamma, vel, m_total)
    res = melt_fraction(m_t, m_i, r_t, r_i, v, b)
    assert 0.0 <= res.f_melt <= 1.0
    assert res.f_melt_sigma >= 0.0 and math.isfinite(res.f_melt_sigma)
    assert res.dIE > 0.0 and math.isfinite(res.dIE)
    assert isinstance(res.in_range, bool) and isinstance(res.above_v_cr, bool)
