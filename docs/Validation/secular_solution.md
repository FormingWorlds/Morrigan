# Validation: secular solution

Source: `src/morrigan/secular_solution.py`. Tests: `tests/test_secular_solution.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_recovers_the_initial_eccentricities_at_epoch` | The returned modes reconstruct the input eccentricities at the reference epoch to 1e-9, down to e = 1e-8 | The Laplace-Lagrange secular solution (Murray & Dermott 1999, ch. 7) |

Epoch recovery is the defining property of the eigen-decomposition, so a single assertion discriminates errors in the interaction matrix, the eigen-solve, the constant fit, and the amplitude scaling at once. The overlap guard branch (zero eigenfrequencies for radially overlapping neighbours), seed determinism, and the amplitude envelope bound are asserted alongside.
