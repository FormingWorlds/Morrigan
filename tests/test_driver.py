"""Tests for ``src/morrigan/driver.py``, the in-memory entry point.

``run_system`` is the interface another program drives the model through,
so what it must guarantee is the shape and the physics of the records it
returns: a per-impact schema whose masses close, whose collision speed
sits above the mutual escape speed, whose geometry is in range, and whose
successive impacts along one body describe that body growing. These are
the same invariants the consuming framework re-checks, pinned here so a
schema regression is caught in this repository rather than downstream.
The file also covers the driver's own plumbing: the seed-to-timeline
determinism contract, the settings-file error path, the event-landing
timestep, and the Hill-radius embryo spacing.
"""

from __future__ import annotations

import numpy as np
import pytest

import morrigan
from morrigan.constants import G, M_earth, M_sun, au2m
from morrigan.driver import allocate_a, read_config, time_step
from morrigan.helper_functions import hill_sphere

pytestmark = [pytest.mark.smoke, pytest.mark.timeout(60)]

#a system that produces several mergers, a body hit more than once, and a body
#left untouched, so the chain and empty-history assertions are all exercised
_MASSES = [0.2, 0.9, 0.4, 1.1, 0.3, 0.7, 1.3, 0.5]
# The bulk density every body is given, and the anchor the recorded radii and
# densities are checked against.
_DENSITY = 5500.0
_SCHEMA = (
    'time', 'M_target_before', 'M_impactor', 'M_merged_after', 'v_impact',
    'v_esc', 'impact_parameter', 'R_target_before', 'R_impactor', 'rho_target',
    'rho_impactor', 'a_before', 'a_after', 'e_after', 'id_target', 'id_impactor',
)


def _run(seed=7, impact_angle=20.0):
    """Evolve one reference system."""
    return morrigan.run_system(
        seed=seed,
        masses=[m * M_earth for m in _MASSES],
        eccentricity=0.05,
        inner_edge=0.05 * au2m,
        spacing=10,
        density=_DENSITY,
        impact_angle=impact_angle,
        evolution_time=1.0,
        inner_cutoff=0.005 * au2m,
        stellar_mass=1.0,
    )


@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_every_impact_record_is_physically_self_consistent():
    """Each returned impact obeys the invariants the consumer will enforce.

    The record is the whole interface, so every one must close its mass as
    a perfect merger, keep its collision speed at or above the mutual
    escape speed, and carry a geometry in range. The escape speed is the
    discriminating check and the analytical anchor: the record value must
    reproduce the closed form sqrt(2 G (M_t + M_i) / (R_t + R_i)) rebuilt
    from the record's own masses and radii, so a record that mixed up a
    mass or a radius would be caught.
    """
    out = _run()
    records = [r for chain in out['impacts'].values() for r in chain]
    assert len(records) >= 4, 'reference system should produce several impacts'

    for r in records:
        assert set(_SCHEMA) <= set(r), f'record is missing fields: {set(_SCHEMA) - set(r)}'
        for key in _SCHEMA:
            assert np.isfinite(r[key]), f'{key} is not finite'

        M_t, M_i = r['M_target_before'], r['M_impactor']
        # Perfect merger: the merged mass the dynamics produced equals the
        # plain sum of the two bodies. The record carries the model's own
        # value rather than re-adding these two fields, so this is a real
        # cross-check: any mass sink in the merger would break it.
        assert r['M_merged_after'] == pytest.approx(M_t + M_i, rel=1e-12)
        # Discrimination: the merged body is strictly heavier than the target
        # alone, by the whole impactor, so reporting the target's own mass
        # unchanged would fail.
        assert abs(r['M_merged_after'] - M_t) > 0.1 * r['M_merged_after']

        # Every extensive quantity is strictly positive.
        for key in ('M_target_before', 'M_impactor', 'M_merged_after', 'v_impact',
                    'v_esc', 'R_target_before', 'R_impactor', 'rho_target',
                    'rho_impactor', 'a_before', 'a_after'):
            assert r[key] > 0.0, f'{key} must be positive'

        # Collision speed cannot fall below the mutual escape speed.
        assert r['v_impact'] >= r['v_esc'] * (1.0 - 1e-9)
        # And that escape speed is the mutual one for this pair. Rebuilding it
        # from the record's own radii cannot catch a swap between them, since
        # R_t + R_i is symmetric, so anchor the radii on the configured bulk
        # density instead: each body's radius must follow from its own mass.
        for mass, radius in ((M_t, r['R_target_before']), (M_i, r['R_impactor'])):
            r_expected = (3.0 * mass / (4.0 * np.pi * _DENSITY)) ** (1.0 / 3.0)
            assert radius == pytest.approx(r_expected, rel=1e-9)
        v_esc_expected = np.sqrt(2.0 * G * (M_t + M_i)
                                 / (r['R_target_before'] + r['R_impactor']))
        assert r['v_esc'] == pytest.approx(v_esc_expected, rel=1e-9)

        # The model carries one bulk density, so both recorded densities must
        # equal the configured value exactly. Building either from the other
        # body's radius leaves them unequal by (M_i/M_t)**2, which feeds the
        # erosion law's (rho_i/rho_t)**0.5 term straight into PROTEUS.
        assert r['rho_target'] == pytest.approx(_DENSITY, rel=1e-9)
        assert r['rho_impactor'] == pytest.approx(_DENSITY, rel=1e-9)

        # The target is the body that survives, so it is the heavier of the
        # pair, and the two identifiers name different bodies.
        assert M_t >= M_i
        assert r['id_target'] != r['id_impactor']

        # The orbit moves: a merger rewrites the survivor's semi-major axis,
        # so reporting a_before twice would leave the coupled planet's orbit
        # frozen at every impact.
        assert r['a_after'] != r['a_before']

        # Geometry and eccentricity stay in their physical ranges.
        assert 0.0 <= r['impact_parameter'] <= 1.0
        assert 0.0 <= r['e_after'] < 1.0

    # A merger leaves the survivor on an eccentric orbit, and the value is
    # carried per impact rather than being a constant. The consumer writes
    # this straight into the coupled planet's eccentricity, so reporting a
    # fixed zero would silently circularise its orbit at every impact while
    # every bound above still held.
    eccentricities = [r['e_after'] for r in records]
    assert any(e > 0.0 for e in eccentricities)
    assert len(set(eccentricities)) > 1


@pytest.mark.physics_invariant
def test_a_body_grows_monotonically_along_its_impact_chain():
    """Consecutive impacts on one body describe it gaining mass, in order.

    A survivor is the target of every impact it appears in, so its chain
    must advance in time and its target mass must pick up where the last
    merger left off. The handover is exact, since nothing removes mass
    between impacts: a body only ever grows by absorbing another.
    """
    out = _run()
    chains = [c for c in out['impacts'].values() if len(c) >= 2]
    assert chains, 'reference system should give at least one multi-impact chain'

    for chain in chains:
        times = [r['time'] for r in chain]
        assert times == sorted(times), 'impacts must be returned in time order'
        assert len(set(times)) == len(times), 'two impacts share a time'
        for earlier, later in zip(chain, chain[1:]):
            # The next target mass equals the previous merged mass exactly.
            assert later['M_target_before'] == pytest.approx(
                earlier['M_merged_after'], rel=1e-12
            )
            # The body is strictly heavier after absorbing an impactor.
            assert later['M_merged_after'] > earlier['M_merged_after']


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


def test_a_head_on_impact_reports_a_zero_impact_parameter():
    """The impact parameter follows the configured impact angle.

    It is the sine of the impact angle, so a head-on run reports zero and a
    grazing run reports one; a run in between is strictly inside the open
    interval. This is the edge behaviour of the geometry the schema carries.
    The 30 degree case is the discriminating one: 0 and 90 degrees are the
    two fixed points of x -> x**2, so squaring the geometry passes both, and
    only an intermediate angle with its value pinned resolves sin from sin**2.
    """
    head_on = _run(impact_angle=0.0)
    grazing = _run(impact_angle=90.0)
    oblique = _run(impact_angle=20.0)
    half = _run(impact_angle=30.0)

    def _b_values(out):
        return [r['impact_parameter'] for chain in out['impacts'].values() for r in chain]

    assert all(b == pytest.approx(0.0, abs=1e-12) for b in _b_values(head_on))
    assert all(b == pytest.approx(1.0, rel=1e-12) for b in _b_values(grazing))
    oblique_b = _b_values(oblique)
    assert oblique_b and all(0.0 < b < 1.0 for b in oblique_b)

    # sin(30 deg) = 0.5 exactly; the squared variant gives 0.25, far outside
    # the tolerance, and the cosine variant gives 0.866.
    half_b = _b_values(half)
    assert half_b and all(b == pytest.approx(0.5, rel=1e-12) for b in half_b)
    assert abs(0.25 - half_b[0]) > 1e-3 and abs(0.8660254 - half_b[0]) > 1e-3


@pytest.mark.physics_invariant
def test_the_same_seed_reproduces_the_timeline_exactly():
    """One seed, one timeline: the run is a pure function of its inputs.

    The determinism contract of the Monte Carlo layer: two runs under
    seed 7 must agree bit for bit in every survivor and every impact
    record, because the consuming framework replays histories and any
    drift would desynchronise its bookkeeping. The discriminating
    counterpart is that seed 8 must produce a different outcome, so the
    equality above cannot be satisfied by the seed being ignored.
    """
    first = _run(seed=7)
    second = _run(seed=7)
    assert first == second

    other = _run(seed=8)
    assert first != other


def test_read_config_names_the_missing_file(tmp_path):
    """A missing settings file fails loudly and names the path it tried.

    Error contract: the default settings path only resolves from a
    checkout root, so the exception must carry the absolute path that
    was tried, which is what makes the failure diagnosable from a batch
    log. The round trip through a real file is the companion check: a
    written TOML table comes back as the same parsed mapping.
    """
    with pytest.raises(FileNotFoundError, match='no/such/settings'):
        read_config('/no/such/settings.toml')

    path = tmp_path / 'settings.toml'
    path.write_text('[init_par]\nN = 3\ne = 0.05\n')
    cfg = read_config(str(path))
    assert cfg['init_par']['N'] == 3
    assert cfg['init_par']['e'] == pytest.approx(0.05, rel=1e-12)


@pytest.mark.physics_invariant
def test_time_step_lands_exactly_on_scheduled_events():
    """The adaptive timestep grows away from events and lands on them.

    Far from any event the step is a tenth of the elapsed time plus a
    hundred-year floor, hand-derived to 3.1536e8 s at t = 0 with the
    model's 365-day year. Close to an event the step is clipped to the
    remaining gap plus one second, so the next time strictly crosses
    the event; at the event itself the step degenerates to exactly one
    second, the edge that keeps the loop advancing.
    """
    dt_far = time_step(0.0, 1e20)
    assert dt_far == pytest.approx(3.1536e8, rel=1e-9)

    dt_near = time_step(1e9, 1e9 + 50.0)
    assert dt_near == pytest.approx(51.0, rel=1e-12)
    assert 1e9 + dt_near >= 1e9 + 50.0  # the step crosses the event

    dt_at = time_step(1e9, 1e9)
    assert dt_at == pytest.approx(1.0, rel=1e-12)
    assert dt_at > 0.0  # the loop always advances


@pytest.mark.physics_invariant
def test_allocate_a_spaces_embryos_by_mutual_hill_radii():
    """Initial embryos are laid out by mutual-Hill-radius spacing.

    For two equal-mass embryos from 0.1 au at spacing 10 the closed
    form is a2 = a1 (1 + 10 ((M1+M2)/(3 Ms))^(1/3)), hand-derived to
    1.6890514743e10 m with the model constants. Spacing zero is the
    degenerate edge, collapsing the system onto one orbit, and a longer
    chain must be strictly increasing: the layout can never fold back
    inward.
    """
    masses = np.array([1.0, 1.0]) * M_earth
    a = allocate_a(2, M_sun, masses, 0.1, 10)
    expected = (0.1 + 10 * 0.1 * ((2 * M_earth) / (3 * M_sun)) ** (1 / 3)) * au2m
    assert a[0] == pytest.approx(0.1 * au2m, rel=1e-12)
    assert a[1] == pytest.approx(expected, rel=1e-12)
    # Discrimination: single-mass Hill spacing (no mutual sum) differs.
    single = (0.1 + 10 * 0.1 * (M_earth / (3 * M_sun)) ** (1 / 3)) * au2m
    assert abs(single - a[1]) > 1e-3 * a[1]

    # Edge: zero spacing stacks every embryo on the inner edge.
    a_zero = allocate_a(3, M_sun, np.array([1.0, 1.0, 1.0]) * M_earth, 0.1, 0)
    assert np.all(a_zero == pytest.approx(0.1 * au2m, rel=1e-12))

    # A five-body chain is strictly increasing.
    chain = allocate_a(5, M_sun, np.ones(5) * M_earth, 0.1, 10)
    assert np.all(np.diff(chain) > 0.0)

    # The mutual Hill radius of the first pair is what sets the first gap.
    r_hill = hill_sphere(0.1, 2 * M_earth, M_sun)
    assert (a[1] - a[0]) == pytest.approx(10 * r_hill * au2m, rel=1e-12)


def test_the_file_writing_path_produces_the_documented_tables(tmp_path):
    """A command-line style run writes the three result tables with the
    documented columns and self-consistent contents.

    The file-writing path is the standalone interface, so one small
    system is run into a temporary directory and the outputs are held
    to their contract: every documented merger column present, each
    merger row closing its dry-run mass sum exactly with a positive
    collision speed, survivor masses
    positive, and the full-system table covering every recorded body.
    The error contract is the schema itself: a renamed or dropped
    column fails here before any consumer sees it.
    """
    from astropy.io import ascii as astropy_ascii

    from morrigan.driver import run_once

    config = {
        'run_simulation': {
            't': 0.0, 't_ref': 0.0, 't_event': 0.0, 'flag_event': 1,
            'a_min': 0.005, 'max_time': 1.0, 'random_seed': 7,
            'save_directory': str(tmp_path),
        },
        'init_par': {
            'N': len(_MASSES), 'e': 0.05, 'impact_angle': 20.0,
            'Mp': list(_MASSES),
            'Ms': 1.0, 'rho_p': 5500.0, 'inner_edge': 0.05, 'spacing': 10,
        },
    }
    summary = run_once(0, config, collect=False)
    assert summary['n_survivors'] >= 1

    mergers = astropy_ascii.read(tmp_path / 'data' / 'mergers' / 'mergers_00.csv',
                                 format='fixed_width')
    expected_cols = ['t', 'id_target', 'id_impactor', 'M_target_before',
                     'M_impactor_before', 'M_merged_after', 'v_c', 'a_final_AU']
    assert list(mergers.colnames) == expected_cols
    assert len(mergers) >= 1
    for row in mergers:
        # A merger is perfect: the merged mass is the exact sum.
        assert row['M_merged_after'] == pytest.approx(
            row['M_target_before'] + row['M_impactor_before'], rel=1e-12
        )
        assert row['v_c'] > 0.0 and row['t'] > 0.0

    survivors = astropy_ascii.read(tmp_path / 'data' / 'survivors' / 'survivors_00.csv',
                                   format='fixed_width')
    assert list(survivors.colnames) == ['id', 'Mp', 'a_AU', 'ecc']
    assert len(survivors) == summary['n_survivors']
    assert all(survivors['Mp'] > 0.0) and all(survivors['a_AU'] > 0.0)

    full = astropy_ascii.read(tmp_path / 'data' / 'full_systems' / 'full_system_00.csv',
                              format='fixed_width')
    # Every surviving body appears in the history, and the clock runs forward.
    assert set(survivors['id']) <= set(full['id'])
    assert min(full['t']) == pytest.approx(0.0, abs=0.0)
