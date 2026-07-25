# Validation: crossing-pair scheduler

Source: `src/morrigan/crossing_pair.py`. Tests: `tests/test_crossing_pair.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_pair_selection_follows_the_perihelion_gap_geometry` | On the unstable triplet (1.0, 1.06, 1.13) au, e = 0.05 on the outer body flips the selection from the inner to the outer pair | Kimura et al. (2025), eq. 4 |

The flip discriminates the perihelion-gap geometry from plain semi-major-axis differences, which order the gaps the other way in this triplet. The two-body branches are pinned alongside: a stable pair defers to the 1e20 sentinel, an unstable pair schedules at exactly 1.5 t, and a packed five-planet system must schedule finite forward events that respond to packing.
