"""
!!! info "`impact_melt.py`"
    Bulk mantle melt fraction of a merger from the Nakajima et al. (2021) scaling law.
    Author(s): Tim Lichtenberg

The law fits SPH giant impacts (Nakajima et al. 2021, EPSL 568, 116983,
doi:10.1016/j.epsl.2021.116983); its merging coefficients come from the runs at
v_imp = v_esc, the campaign reaches 2 v_esc. This module ports the bulk part, Eqs. 2-4,
6, 7 and 14, with the coefficients of the reference code MeltScalingLaw at commit
2eca30f2 (https://github.com/mikinakajima/MeltScalingLaw, MIT), which add 45 degree
impacts and differ from the paper's Tables S.5 and S.6; ``data/nakajima2021/`` holds them
unchanged.
"""

from __future__ import annotations

import math
from functools import lru_cache
from importlib import resources
from numbers import Real
from typing import NamedTuple

import numpy as np
from numpy.polynomial import legendre

from morrigan.constants import G

M_MARS = 6.4171e23  # mass unit of the fit [kg] (Table 1)
E_M = 5.2e6  # specific energy that melts the mantle after decompression [J kg-1] (Table 1)
F_MANTLE = 0.70  # mantle mass fraction of the simulated bodies
FIT_ANGLES = (0.0, 30.0, 45.0, 60.0, 90.0)  # simulated impact angles [deg]
MASS_RANGE = (1.0, 53.0)  # total mass [M_Mars]
GAMMA_RANGE = (0.03, 0.5)  # impactor-to-total mass ratio
SPEED_RANGE = (1.0, 2.0)  # impact speed over the mutual escape speed
SPEED_MARGIN = 1e-9  # relative, so a contact speed equal to the recomputed v_esc is in range
_P6_AS_FITTED = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0])


class MeltResult(NamedTuple):
    """Bulk melt of one merger.

    Attributes
    ----------
    f_melt : float
        Melt mass fraction of the fitted post-impact mantle, 0.70 m(theta) M with m the
        Eq. 7 fit (0.924 to 1.0015), Eq. 14, clipped at 1.
    f_melt_raw : float
        The Eq. 14 fit before the clip, the value ``f_melt_sigma`` refers to; above 1 the
        released energy exceeds what melts the whole mantle.
    f_melt_sigma : float
        1-sigma of the unclipped fraction, from the fit errors at the nearest simulated
        angle; where ``f_melt`` is clipped at 1 it belongs to the unclipped value and can
        exceed 1.
    dIE : float
        Internal energy gain of the merged body, Eq. 4 [J].
    h : float
        Share of ``dIE`` deposited in the mantle, Eq. 6.
    in_range : bool
        Whether the merger lies in the fitted mass, mass-ratio and speed range.
    above_v_cr : bool
        Whether the impact speed exceeds the Genda et al. (2012) critical speed, where the
        reference code switches to its hit-and-run fit.
    """

    f_melt: float
    f_melt_raw: float
    f_melt_sigma: float
    dIE: float
    h: float
    in_range: bool
    above_v_cr: bool


@lru_cache(maxsize=1)
def _coefficients() -> tuple[np.ndarray, np.ndarray]:
    """Return the merging-branch fit coefficients and the 1-sigma rows at v_imp = v_esc."""
    data = resources.files('morrigan') / 'data' / 'nakajima2021'
    para = np.array((data / 'parameter.txt').read_text().splitlines()[0].split(), dtype=float)
    rows = [
        line.split('#')[0].split()
        for line in (data / 'Model_sigma.txt').read_text().splitlines()
    ]
    sigma = np.array(rows[1:4], dtype=float)  # IE, h and mantle mass, one column per angle
    if para.shape != (12,) or sigma.shape != (3, len(FIT_ANGLES)):
        raise RuntimeError(
            f'Nakajima coefficient files have shapes {para.shape}, {sigma.shape}'
        )
    return para, sigma


def _check(name: str, value, low: float = 0.0, high: float | None = None) -> float:
    """Return ``value`` as a float, or raise ValueError naming the argument.

    A value must be a finite real number (not a bool) above ``low``, or in
    ``[low, high]`` when ``high`` is given.
    """
    try:
        ok = isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)
    except OverflowError:  # an integer beyond the float range
        ok = False
    ok = ok and (low <= value <= high if high is not None else value > low)
    if not ok:
        bound = f'in [{low}, {high}]' if high is not None else f'finite and above {low}'
        raise ValueError(f'{name} must be a real number {bound}, got {value!r}')
    return float(value)


def melt_fraction(
    M_target,
    M_impactor,
    R_target,
    R_impactor,
    v_impact,
    impact_parameter,
    R_merged=None,
) -> MeltResult:
    """Return the bulk mantle melt fraction of a perfect merger (Nakajima et al. 2021).

    The released energy is the kinetic energy at contact, Eq. 2, plus the potential energy
    of the two bodies minus that of the merged body, Eq. 3. Legendre fits in
    cos(theta), with sin(theta) the impact parameter, give the share of it that becomes
    internal energy (Eq. 4), the mantle's share of that (Eq. 6) and the post-impact mantle
    mass (Eq. 7); the melt fraction is the mantle's energy gain over the mantle mass and
    E_M (Eq. 14). Every merger uses the merging coefficients, as Morrigan merges every
    collision; ``above_v_cr`` flags those the reference code would fit as hit-and-run.
    Outside the fitted range the value is the fit's, flagged by ``in_range``. The 1-sigma
    uses the fit errors of the simulated angle nearest to theta (ties to the lower one),
    a rule of this port: the reference code accepts only the five simulated angles.

    Parameters
    ----------
    M_target, M_impactor : float
        Masses of the target (the larger body) and the impactor [kg].
    R_target, R_impactor : float
        Their radii [m].
    v_impact : float
        Speed at first contact [m s-1].
    impact_parameter : float
        sin(theta) in [0, 1]; 0 is head-on, 1 is grazing.
    R_merged : float, optional
        Radius of the merged body [m]; by default the target's density at the merged mass,
        as Morrigan's own radius update gives it.

    Returns
    -------
    MeltResult
        The melt fraction, its 1-sigma, the internal energy gain, the mantle share and the
        two flags.

    Raises
    ------
    ValueError
        If an argument is not a finite real number in its domain, the impactor is the
        heavier body, or the inputs release no finite positive energy (a merged radius
        too large, or a speed whose energy overflows).
    RuntimeError
        If a fitted factor is not positive or the result leaves [0, 1]; the shipped fits
        are positive over all angles, so this signals broken coefficient files.
    """
    m_t = _check('M_target', M_target)
    m_i = _check('M_impactor', M_impactor)
    r_t = _check('R_target', R_target)
    r_i = _check('R_impactor', R_impactor)
    v = _check('v_impact', v_impact)
    b = _check('impact_parameter', impact_parameter, 0.0, 1.0)
    if m_i > m_t:
        raise ValueError(
            f'M_impactor {m_i!r} exceeds M_target {m_t!r}; the target is the larger body'
        )
    m = m_t + m_i
    r_m = r_t * (m / m_t) ** (1.0 / 3.0) if R_merged is None else _check('R_merged', R_merged)

    q = m_i / m_t
    scale = G * m_t**2 / r_t
    d_pe = (
        -0.6
        - 0.6 * q**2 / (r_i / r_t)
        - q / (1.0 + r_i / r_t)
        + 0.6 * (1.0 + q) ** 2 / (r_m / r_t)
    ) * scale
    try:
        ke = 0.5 * m_t * m_i / m * v**2
    except OverflowError:
        raise ValueError(f'v_impact {v!r} overflows the kinetic energy') from None
    energy = d_pe + ke
    if not (math.isfinite(energy) and energy > 0.0):
        raise ValueError(
            f'the merged radius R_merged {r_m!r} and v_impact {v!r} release no finite '
            f'positive energy ({energy} J)'
        )
    para, sigma = _coefficients()
    x = math.sqrt(1.0 - b * b)  # cos(theta)
    h = float(legendre.legval(x, para[0:3]))
    # The reference code's sixth-order term is twice the standard P6 (1/8, not 1/16); its
    # coefficients come with that basis, so the last IE coefficient is doubled here.
    ie = float(legendre.legval(x, para[3:10] * _P6_AS_FITTED))
    mantle = float(legendre.legval(x, para[10:12]))
    if not (h > 0.0 and ie > 0.0 and mantle > 0.0):
        raise RuntimeError(f'non-positive fit factor: h {h}, IE {ie}, mantle {mantle}')
    d_ie = ie * energy
    f_raw = h * ie * energy / (F_MANTLE * mantle * m) / E_M

    # nearest simulated angle; a tie, within the round-off of asin(sin), goes to the lower
    dist = [abs(math.degrees(math.asin(b)) - a) for a in FIT_ANGLES]
    k = next(j for j, d in enumerate(dist) if d <= min(dist) + 1e-9)
    dz = math.sqrt(
        (sigma[0, k] / ie) ** 2 + (sigma[1, k] / h) ** 2 + (sigma[2, k] / mantle) ** 2
    )
    f_sigma = f_raw * dz
    f = min(f_raw, 1.0)
    if not (0.0 <= f <= 1.0 and f_sigma >= 0.0 and math.isfinite(f_sigma)):
        raise RuntimeError(f'melt fraction {f} with sigma {f_sigma} leaves its bounds')

    ratio = v / math.sqrt(2.0 * G * m / (r_t + r_i))
    big_gamma, grazing = (m_t - m_i) / m, 1.0 - b
    v_cr = 2.43 * big_gamma**2 * grazing**2.5 - 0.0408 * big_gamma + 1.86 * grazing**2.5 + 1.08
    in_range = (
        MASS_RANGE[0] <= m / M_MARS <= MASS_RANGE[1]
        and GAMMA_RANGE[0] <= m_i / m <= GAMMA_RANGE[1]
        and SPEED_RANGE[0] * (1.0 - SPEED_MARGIN) <= ratio <= SPEED_RANGE[1]
    )
    return MeltResult(f, f_raw, f_sigma, d_ie, h, in_range, ratio > v_cr)
