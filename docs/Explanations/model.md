# Model overview

Morrigan implements the semi-analytical Monte Carlo model of [Kimura et al. (2025)](https://iopscience.iop.org/article/10.3847/1538-4357/ade992/meta) for the giant-impact phase of terrestrial planet formation. Instead of integrating orbits directly, the model predicts *when* a system of embryos becomes unstable, resolves each instability into a collision, a scattering, or an ejection with analytic prescriptions calibrated against N-body simulations, and steps the system forward between events with the linear secular solution. Equation numbers below refer to Kimura et al. (2025).

## System setup

A system starts as $N$ embryos on nearly circular orbits, spaced by a fixed number of mutual Hill radii,

$$ r_\mathrm{H} = \left( \frac{M_i + M_{i+1}}{3 M_\ast} \right)^{1/3} \frac{a_i + a_{i+1}}{2}, $$

from an inner edge outward. Each embryo carries a mass, a bulk density (setting its radius), an eccentricity, and an atmospheric mass fraction used by the internal loss bookkeeping. All state is SI internally; the settings file is read once in Earth masses, au, degrees, and Gyr.

## Secular evolution between events

Between events the eccentricities evolve under the classical Laplace-Lagrange secular solution (e.g. Murray & Dermott 1999, ch. 7), computed with Laplace coefficients from the `pylaplace` package. The eigen-decomposition of the interaction matrix yields the secular frequencies, and the integration constants are fitted to the state at the last event, so the reconstructed eccentricities recover the epoch state exactly. Radially overlapping neighbours are excluded from the coupling, where the linear expansion is invalid.

## Timing the next instability

The time to the next orbit crossing is evaluated for every adjacent triplet from the resonance-overlap crossing time of Petit et al. (2020), scaled up from the three-planet case with a resonance-density factor that counts how many neighbours are packed within the critical spacing (eq. 5). A two-planet pair is first screened with the stability criterion of eq. 28: a pair whose eccentricities and separation put it on the stable side never crosses. The interacting pair within a triplet is the one with the smaller physical gap, the outer body's perihelion minus the inner body's aphelion (eq. 4). The scheduled event time adds the crossing time and the shorter of two interaction timescales: the viscous stirring time (eqs. 8-10, 23) and the collision time (eq. 11).

## Resolving an event

When a crossing falls due, the pair collides with probability (eqs. 12, 14)

$$ p_\mathrm{col} = 1 - e^{-\lambda}, \qquad \lambda = \left( \frac{2 e_{ij}}{e_\mathrm{esc}} \right)^2 \left( 1 + \frac{e_{ij}^2}{e_\mathrm{esc}^2} \right) \frac{1}{\ln \Lambda}, $$

where $e_{ij}$ is the pair's relative eccentricity and $e_\mathrm{esc}$ the mutual escape eccentricity, the surface escape speed of the pair over the local Kepler speed. A uniform draw against $p_\mathrm{col}$ selects the branch:

- **Collision.** The eccentricities are re-drawn from a truncated Rayleigh distribution until the epicycles geometrically overlap, the contact speed is $v_c = \sqrt{v_\infty^2 + v_\mathrm{esc}^2}$, and the pair merges (see below).
- **Scattering.** Both orbits are shifted apart by the excited epicycle amplitudes, weighted with the opposite body's mass fraction (eqs. 18-20), which conserves the mass-weighted sum of semi-major axes exactly.
- **Ejection.** If scattering excites an eccentricity to or past unity, the more excited body escapes; the survivor is placed on the tighter orbit that conserves the pair's orbital energy, and the escaping body is removed.

A body whose perihelion falls inside the inner cutoff is removed as having fallen into the star.

## Merging

A merger (eqs. 15-17) sums the masses, places the merged body at

$$ a_\mathrm{new} = \frac{M_1 + M_2}{M_1/a_1 + M_2/a_2}, $$

the orbital-energy-weighted orbit, and composes the eccentricity as a mass-weighted vector sum with a pericentre alignment drawn uniformly outside the geometrically forbidden range. The heavier body survives; the lighter one is removed. The combined atmosphere of the pair is then reduced by the [Kegerreis et al. (2020) loss fraction](mass_loss.md), clamped to $[0, 1]$, with the rock mass conserved exactly.

## Randomness and reproducibility

All Monte Carlo draws go through numpy's global random state, seeded once per system from the settings-file seed plus the system index, so a given settings file reproduces its systems exactly. `run_system` additionally saves and restores the caller's random state, so an embedding program's own draws are unaffected.

## References

- Kimura, T., et al. (2025), ApJ, 989, 109 (the model; "Paper I")
- Petit, A. C., Pichierri, G., Davies, M. B., & Johansen, A. (2020), A&A, 641, A176 (three-planet crossing time)
- Kegerreis, J. A., et al. (2020), ApJL, 901, L31 (impact atmosphere loss)
- Murray, C. D., & Dermott, S. F. (1999), Solar System Dynamics (Laplace-Lagrange secular theory)
