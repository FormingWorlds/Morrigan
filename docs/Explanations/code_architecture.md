# Code architecture

Eleven source files under `src/morrigan/`, laid out so that each stage of the physics has one home. This page maps them onto what the [model overview](model.md) describes.

## The call graph

```
   morrigan -c settings.toml
        |
        v
   driver.run_batch                 parallel over systems
        |
        v
   driver.run_once                  ONE system, start to finish
        |
        +-- allocate_a              lay out embryos by mutual Hill spacing
        |
        +-- while t <= max_time:
        |     |
        |     +-- sort_planet       drop the dead, re-sort by orbit
        |     |
        |     +-- secular_solution  eccentricities between events
        |     |
        |     +-- crossing_pair     which pair crosses next, and when
        |     |     |
        |     |     +-- interaction_timescales   viscous and collision times
        |     |
        |     +-- orbit_cross_K25   resolve the event
        |           |
        |           +-- merge_embryo         if they collide
        |           +-- (scatter in place)   if they do not
        |           +-- (ejection branch)    if one is excited past e = 1
        |
        +-- write tables / return records
```

`run_system` is the in-process entry point beside `run_batch`. It drives the same `run_once` core but returns the records in memory instead of writing files, which is what PROTEUS consumes.

## What each file owns

| File | Owns | Physics or plumbing |
| --- | --- | --- |
| `driver.py` | Orchestration: settings, layout, the event loop, the record schema, survivor selection, output | physics |
| `crossing_pair.py` | Which adjacent pair goes unstable next and when | physics |
| `interaction_timescales.py` | Viscous stirring and collision timescales, the stability criterion | physics |
| `secular_solution.py` | Laplace-Lagrange eigenmodes between events | physics |
| `orbit_cross_K25.py` | Resolving a crossing into a merger, a scattering or an ejection | physics |
| `merge_embryo.py` | Merger bookkeeping: masses, orbits, eccentricities | physics |
| `helper_functions.py` | Kepler period, Hill radius, escape eccentricity, Rayleigh draws | physics |
| `sort_planet.py` | Dead-planet removal and re-sorting | plumbing |
| `constants.py` | Physical constants and unit conversions | plumbing |
| `__init__.py` | Public API: `run_system` | plumbing |
| `_version.py` | Generated from git tags | plumbing |

The split matters for testing: every physics file must carry a test that pins it against a published value, an analytical limit or a cross-implementation check, recorded under [Validation anchors](../Validation/index.md). Plumbing files are exempt from that requirement but not from the rest of the test contract.

## State and where it lives

There is no state object. A system is a set of parallel numpy arrays, one entry per body, threaded through the loop as arguments:

`a` (semi-major axes), `masses`, `ecc`, `Rp` (radii), `live_status`, `interact`, `densities`, `planet_id`.

They are mutated in place by the event functions, which is why `sort_planet` can drop the dead bodies from all of them in one shared permutation, and why a test that reads an input array after calling an event function is reading post-mutation state.

`planet_id` is the identity that survives sorting: array positions change as bodies die and the system is re-sorted, so it is the id, not the index, that appears in the output records.

## Randomness

All Monte Carlo draws go through numpy's global random state, seeded once per system at the top of `run_once`. The seed is mixed from the settings-file seed and the system's index, so a batch's systems are independent of each other and two batches with neighbouring seeds do not overlap.

`run_system` additionally saves and restores the caller's random state around its own draws, so running a system in-process does not perturb an embedding program's sequence. That matters in a coupled run, where several modules share the same global state.

## Units

SI internally: kg, m, s. The settings file is read in Earth masses, AU, degrees and Gyr, converted once in `driver.py`. The returned records are SI apart from times, which are in years using the model's own 365-day year. A physical constant appearing as a literal in a function body rather than coming from `constants.py` is a defect.
