"""Cross-cutting ensemble statistics of the full model.

This file is the documented exception to the one-test-file-per-source
rule: its subject is the statistical behaviour of whole seed ensembles,
which no single source file owns. The anchor is Kimura et al. (2025),
whose reference model S0 (Table 1: 15 embryos, 2.43 Earth masses in
total, from 0.1 au at 10 mutual-Hill-radii spacing, bulk density
3 g cm-3, solar-mass host) is reproduced here in an equal-mass variant
and evolved over 20 seeds. Two anchors are pinned: the closed-form
eq. 37 coefficient relating the crossing eccentricity of a pair to its
escape eccentricity, and the Hill-scaled endstate statistics (eqs.
33-34) of the evolved ensembles.
"""

from __future__ import annotations

import numpy as np
import pytest

import morrigan
from morrigan.constants import M_earth, M_sun, au2m
from morrigan.driver import run_once
from morrigan.helper_functions import esc_ecc, planet_radius

pytestmark = [pytest.mark.integration, pytest.mark.timeout(300)]

_N_INI = 15
_M_EACH = 2.43 / _N_INI * M_earth  # S0 total mass split into equal embryos
_RHO = 3000.0
_SEEDS = range(20)


def _s0_config(seed):
    """The equal-mass S0-like settings dictionary for one seed."""
    return {
        'run_simulation': {
            't': 0.0, 't_ref': 0.0, 't_event': 0.0, 'flag_event': 1,
            'a_min': 0.02, 'max_time': 0.03, 'random_seed': seed,
            'save_directory': '',
        },
        'init_par': {
            'N': _N_INI, 'e': 0.01, 'impact_angle': 45.0,
            'Mp': [_M_EACH / M_earth] * _N_INI,
            'atm_mass_fraction': [0.0] * _N_INI,
            'Ms': 1.0, 'rho_p': _RHO, 'inner_edge': 0.1, 'spacing': 10.0,
        },
    }


def _hill_statistics(a, m, e):
    """Mean Hill-scaled separation and eccentricity, eqs. 33-34 of
    Kimura et al. (2025), over adjacent surviving pairs."""
    order = np.argsort(a)
    a, m, e = a[order], m[order], e[order]
    b_h, e_h = [], []
    for i in range(len(a) - 1):
        r_hill = ((m[i] + m[i + 1]) / (3.0 * M_sun)) ** (1 / 3) * (a[i] + a[i + 1]) / 2.0
        b_h.append((a[i + 1] - a[i]) / r_hill)
        e_h.append((e[i] * a[i] + e[i + 1] * a[i + 1]) / (2.0 * r_hill))
    return float(np.mean(b_h)), float(np.mean(e_h))


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_eq37_coefficient_matches_the_published_value():
    """The crossing-to-escape eccentricity ratio reproduces the 0.07
    coefficient of eq. 37 in Kimura et al. (2025).

    Published benchmark: for equal-mass neighbours spaced b_H mutual
    Hill radii apart, eq. 37 states e_ij^2 / e_esc^2 = 0.07 (b_H/10)^2
    at 1 au, bulk density 3 g cm-3, and a solar-mass host. Building the
    same ratio from the model's own primitives (the eq. 36 crossing
    eccentricity over the mutual escape eccentricity) gives 0.0688,
    matching the published coefficient at its stated precision. The
    ratio must also be independent of the embryo mass, which cancels
    analytically; a wrong mass exponent anywhere in the chain would
    break that cancellation across the hundredfold mass sweep.
    """
    b_hill = 10.0
    ratios = []
    for mass in (0.1 * M_earth, M_earth, 10.0 * M_earth):
        a = au2m
        radius = planet_radius(mass, _RHO)
        r_hill = ((2.0 * mass) / (3.0 * M_sun)) ** (1 / 3) * a  # eq. 34, a_i = a_j
        separation = b_hill * r_hill
        e_ij = np.sqrt(2.0) * separation / (2.0 * a)  # eq. 36
        e_escape = esc_ecc(M_sun, mass, mass, radius, radius, a)
        ratios.append(e_ij**2 / e_escape**2)

    assert ratios[1] == pytest.approx(0.0688, rel=1e-3)
    # Consistency with the published rounded coefficient 0.07.
    assert ratios[1] == pytest.approx(0.07, rel=0.03)
    # Analytic mass cancellation across two orders of magnitude.
    assert ratios[0] == pytest.approx(ratios[1], rel=1e-12)
    assert ratios[2] == pytest.approx(ratios[1], rel=1e-12)
    # Discrimination: the (b_H/10)^2 dependence, not linear in spacing.
    e_ij_double = np.sqrt(2.0) * (2.0 * b_hill * r_hill) / (2.0 * au2m)
    assert (e_ij_double**2 / e_escape**2) == pytest.approx(4.0 * ratios[2], rel=1e-12)


@pytest.mark.physics_invariant
def test_s0_like_ensembles_reach_a_widened_hill_stable_endstate():
    """Twenty S0-like seeds all accrete and settle into widened,
    moderately excited endstates.

    Every seeded system must lose bodies to mergers (final count
    strictly below 15 and at least 1) while never gaining mass overall.
    The endstate statistics of eqs. 33-34 pin the physics: the mean
    Hill-scaled separation must exceed the initial spacing of 10, since
    giant-impact evolution runs until neighbours are Hill-stable, and
    the mean Hill-scaled eccentricity must sit in the moderately
    excited band (0.5 to 15, ensemble mean 1 to 8) set below by the
    escape-eccentricity scale of the region and above zero by the
    excitation every scattering leaves behind. A model that stopped
    merging, over-damped, or over-excited would leave the band.
    """
    n_final, b_h_means, e_h_means = [], [], []
    total_initial = _N_INI * _M_EACH

    for seed in _SEEDS:
        raw = run_once(0, _s0_config(seed), collect=True)
        n_final.append(len(raw['masses']))

        assert 1 <= len(raw['masses']) < _N_INI
        assert float(np.sum(raw['masses'])) <= total_initial * (1.0 + 1e-12)
        assert np.all(raw['masses'] > 0.0) and np.all(raw['a'] > 0.0)

        if len(raw['masses']) >= 2:
            b_h, e_h = _hill_statistics(raw['a'], raw['masses'], raw['ecc'])
            b_h_means.append(b_h)
            e_h_means.append(e_h)
            assert b_h > 10.0  # endstates widen beyond the initial spacing
            assert 0.5 < e_h < 15.0

    # The ensemble as a whole accretes into a few final bodies.
    assert 3.0 <= float(np.mean(n_final)) <= 8.0
    assert 1.0 < float(np.mean(e_h_means)) < 8.0
    # The endstate is not degenerate across seeds.
    assert max(n_final) > min(n_final)


@pytest.mark.physics_invariant
def test_dry_impact_chains_hand_over_masses_across_the_ensemble():
    """Every impact chain in every seeded system hands its mass over
    exactly and runs strictly forward in time.

    Run dry, the target mass of each next impact must equal the
    previous perfect-merger sum to machine precision, for every
    survivor of every seed: the chain-continuity contract the consuming
    framework validates per impact. Five seeds keep the sweep light;
    each must produce at least one multi-impact chain for the handover
    clause to bite, which the S0-like packing guarantees.
    """
    checked_chains = 0
    for seed in list(_SEEDS)[:5]:
        out = morrigan.run_system(
            seed=seed, masses=[_M_EACH] * _N_INI, eccentricity=0.01,
            inner_edge=0.1 * au2m, spacing=10.0, density=_RHO,
            impact_angle=45.0, evolution_time=0.03,
            inner_cutoff=0.02 * au2m, stellar_mass=1.0,
        )
        for chain in out['impacts'].values():
            times = [r['time'] for r in chain]
            assert times == sorted(times)
            for earlier, later in zip(chain, chain[1:]):
                assert later['M_target_before'] == pytest.approx(
                    earlier['M_merged_after'], rel=1e-12
                )
                checked_chains += 1
    assert checked_chains >= 5, 'the ensemble must exercise real multi-impact chains'
