"""Tests for the in-memory entry point in ``morrigan.driver``.

``run_system`` is the interface another program drives the model through,
so what it must guarantee is the shape and the physics of the records it
returns: a per-impact schema whose masses close, whose collision speed
sits above the mutual escape speed, whose geometry is in range, and whose
successive impacts along one body describe that body growing. These are
the same invariants the consuming framework re-checks, pinned here so a
schema regression is caught in this repository rather than downstream.
"""

from __future__ import annotations

import numpy as np
import pytest

import morrigan
from morrigan.constants import G, M_earth, au2m

pytestmark = [pytest.mark.unit, pytest.mark.timeout(60)]

#a system that produces several mergers, a body hit more than once, and a body
#left untouched, so the chain and empty-history assertions are all exercised
_MASSES = [0.2, 0.9, 0.4, 1.1, 0.3, 0.7, 1.3, 0.5]
_SCHEMA = (
    'time', 'M_target_before', 'M_impactor', 'M_merged_after', 'v_impact',
    'v_esc', 'impact_parameter', 'R_target_before', 'R_impactor', 'rho_target',
    'rho_impactor', 'a_before', 'a_after', 'e_after', 'id_target', 'id_impactor',
)


def _run(seed=7, atm_mass_fraction=0.0, impact_angle=20.0):
    """Evolve one reference system dry, unless an atmosphere is asked for."""
    return morrigan.run_system(
        seed=seed,
        masses=[m * M_earth for m in _MASSES],
        eccentricity=0.05,
        inner_edge=0.05 * au2m,
        spacing=10,
        density=5500.0,
        impact_angle=impact_angle,
        evolution_time=1.0,
        inner_cutoff=0.005 * au2m,
        stellar_mass=1.0,
        atm_mass_fraction=atm_mass_fraction,
    )


@pytest.mark.unit
@pytest.mark.physics_invariant
def test_every_impact_record_is_physically_self_consistent():
    """Each returned impact obeys the invariants the consumer will enforce.

    The record is the whole interface, so every one must close its mass as
    a perfect merger, keep its collision speed at or above the mutual
    escape speed, and carry a geometry in range. The escape speed is the
    discriminating check: it is rebuilt here from the masses and radii and
    must match the value implied by the reported collision speed, so a
    record that mixed up a mass or a radius would be caught.
    """
    out = _run()
    records = [r for chain in out['impacts'].values() for r in chain]
    assert len(records) >= 4, 'reference system should produce several impacts'

    for r in records:
        assert set(_SCHEMA) <= set(r), f'record is missing fields: {set(_SCHEMA) - set(r)}'
        for key in _SCHEMA:
            assert np.isfinite(r[key]), f'{key} is not finite'

        M_t, M_i = r['M_target_before'], r['M_impactor']
        # Perfect merger: the reported merged mass is the plain sum, exactly.
        assert r['M_merged_after'] == pytest.approx(M_t + M_i, rel=1e-12)
        # Discrimination: reporting the model's post-loss target mass instead
        # would, even dry, differ by the whole impactor mass, ~30% here.
        assert abs(r['M_merged_after'] - M_t) > 0.1 * r['M_merged_after']

        # Every extensive quantity is strictly positive.
        for key in ('M_target_before', 'M_impactor', 'M_merged_after', 'v_impact',
                    'v_esc', 'R_target_before', 'R_impactor', 'rho_target',
                    'rho_impactor', 'a_before', 'a_after'):
            assert r[key] > 0.0, f'{key} must be positive'

        # Collision speed cannot fall below the mutual escape speed.
        assert r['v_impact'] >= r['v_esc'] * (1.0 - 1e-9)
        # And that escape speed is the mutual one for this pair; a wrong mass
        # or radius in the record would not reproduce it to this tolerance.
        v_esc_expected = np.sqrt(2.0 * G * (M_t + M_i)
                                 / (r['R_target_before'] + r['R_impactor']))
        assert r['v_esc'] == pytest.approx(v_esc_expected, rel=1e-9)

        # Geometry and eccentricity stay in their physical ranges.
        assert 0.0 <= r['impact_parameter'] <= 1.0
        assert 0.0 <= r['e_after'] < 1.0


@pytest.mark.unit
@pytest.mark.physics_invariant
def test_a_body_grows_monotonically_along_its_impact_chain():
    """Consecutive impacts on one body describe it gaining mass, in order.

    A survivor is the target of every impact it appears in, so its chain
    must advance in time and its target mass must pick up where the last
    merger left off. Run dry, that handover is exact, which is the sharp
    discriminator: it holds only because the merged mass is reported as the
    plain sum, so a naive post-loss report would break the equality.
    """
    out = _run()
    chains = [c for c in out['impacts'].values() if len(c) >= 2]
    assert chains, 'reference system should give at least one multi-impact chain'

    for chain in chains:
        times = [r['time'] for r in chain]
        assert times == sorted(times), 'impacts must be returned in time order'
        assert len(set(times)) == len(times), 'two impacts share a time'
        for earlier, later in zip(chain, chain[1:]):
            # Dry: the next target mass equals the previous merged mass exactly.
            assert later['M_target_before'] == pytest.approx(
                earlier['M_merged_after'], rel=1e-12
            )
            # The body is strictly heavier after absorbing an impactor.
            assert later['M_merged_after'] > earlier['M_merged_after']


@pytest.mark.unit
def test_impacts_are_keyed_only_by_survivors_and_each_is_present():
    """The impact histories belong to survivors and to nobody else.

    A body that is a target early and then dies leaves a partial history
    that is not a survivor's and can even hold an unbound orbit, so the
    returned histories must be exactly the survivors, each present so a
    caller can always look one up, empty if that body never merged.
    """
    out = _run()
    survivor_ids = {s['id'] for s in out['survivors']}
    assert set(out['impacts']) == survivor_ids
    # A body that never merged is still queryable, with an empty history.
    never_hit = [sid for sid in survivor_ids if not out['impacts'][sid]]
    assert never_hit, 'this system leaves at least one body untouched'
    # No survivor's own chain carries an unbound post-merge orbit.
    for chain in out['impacts'].values():
        assert all(r['e_after'] < 1.0 for r in chain)


@pytest.mark.unit
@pytest.mark.physics_invariant
def test_survivors_report_initial_and_final_state_consistently():
    """Every survivor carries a finite, positive initial and final state.

    Survivor selection downstream reads these, so a survivor must expose a
    positive initial and final mass and orbit, and its final mass must be
    at least its initial mass, because a surviving body only ever accretes.
    """
    out = _run()
    assert out['survivors'], 'reference system leaves survivors'
    for s in out['survivors']:
        for key in ('mass_initial', 'a_initial', 'mass_final', 'a_final'):
            assert np.isfinite(s[key]) and s[key] > 0.0, f'{key} must be finite and positive'
        # A survivor never loses net mass over the run (dry dynamics).
        assert s['mass_final'] >= s['mass_initial'] * (1.0 - 1e-9)


@pytest.mark.unit
def test_a_head_on_impact_reports_a_zero_impact_parameter():
    """The impact parameter follows the configured impact angle.

    It is the sine of the impact angle, so a head-on run reports zero and a
    grazing run reports one; a run in between is strictly inside the open
    interval. This is the edge behaviour of the geometry the schema carries.
    """
    head_on = _run(impact_angle=0.0)
    grazing = _run(impact_angle=90.0)
    oblique = _run(impact_angle=20.0)

    def _b_values(out):
        return [r['impact_parameter'] for chain in out['impacts'].values() for r in chain]

    assert all(b == pytest.approx(0.0, abs=1e-12) for b in _b_values(head_on))
    assert all(b == pytest.approx(1.0, rel=1e-12) for b in _b_values(grazing))
    oblique_b = _b_values(oblique)
    assert oblique_b and all(0.0 < b < 1.0 for b in oblique_b)
