# Model overview

Morrigan implements the semi-analytical Monte Carlo model of Kimura et al. (2025)[^cite-kimura2025] for the giant-impact phase of terrestrial planet formation. Instead of integrating orbits directly, the model predicts *when* a system of embryos becomes unstable, resolves each instability into a collision, a scattering, or an ejection with analytic prescriptions calibrated against N-body simulations, and steps the system forward between events with the linear secular solution. Equation numbers below refer to that paper, which does not treat ejection; the ejection branch follows the companion application paper[^cite-kimura2025b].

## System setup

A system starts as $N$ embryos on nearly circular orbits, spaced by a fixed number of mutual Hill radii,

$$ r_\mathrm{H} = \left( \frac{M_i + M_{i+1}}{3 M_\ast} \right)^{1/3} \frac{a_i + a_{i+1}}{2}, $$

laid out from an inner edge outward, each orbit placed a fixed multiple of
$r_\mathrm{H}$ beyond the one inside it. The radius is evaluated at the mean of
the pair's orbits, the standard mutual Hill radius, so the outer orbit appears
on both sides of the spacing condition and the layout follows from solving it:
$a_{i+1} = a_i (1 + sC/2)/(1 - sC/2)$ with $C = ((M_i + M_{i+1})/3M_\ast)^{1/3}$.
Every gap then measures exactly the configured number of mutual Hill radii.

Each embryo carries a mass, a bulk density (setting its radius), and an
eccentricity. Bodies are bare: the model tracks no atmosphere, so a mass is
always a bulk mass and a radius always follows from it at the configured
density. All state is SI internally; the settings file is read once in Earth
masses, au, degrees, and Gyr.

## Secular evolution between events

Between events the eccentricities evolve under the classical Laplace-Lagrange secular solution (e.g. Murray & Dermott 1999, ch. 7[^cite-murray1999]), computed with Laplace coefficients from the `pylaplace` package. The eigen-decomposition of the interaction matrix yields the secular frequencies, and the integration constants are fitted to the state at the last event, so the reconstructed eccentricities recover the epoch state exactly. Radially overlapping neighbours are excluded from the coupling, where the linear expansion is invalid.

## Timing the next instability

The time to the next orbit crossing is evaluated for every adjacent triplet from the resonance-overlap crossing time of Petit et al. (2020)[^cite-petit2020], scaled up from the three-planet case with a resonance-density factor that counts how many neighbours are packed within the critical spacing (eq. 5). A two-planet pair is first screened with the stability criterion of eq. 28: a pair whose eccentricities and separation put it on the stable side never crosses. The interacting pair within a triplet is the one with the smaller physical gap, the outer body's perihelion minus the inner body's aphelion (eq. 4). The scheduled event time adds the crossing time and the shorter of two interaction timescales: the viscous stirring time (eqs. 8-10, 23) and the collision time (eq. 11).

## Resolving an event

When a crossing falls due, the pair collides with probability (eqs. 12, 14)

$$ p_\mathrm{col} = 1 - e^{-\lambda}, \qquad \lambda = \left( \frac{2 e_{ij}}{e_\mathrm{esc}} \right)^2 \left( 1 + \frac{e_{ij}^2}{e_\mathrm{esc}^2} \right) \frac{1}{\ln \Lambda}, $$

where $e_{ij}$ is the pair's relative eccentricity and $e_\mathrm{esc}$ the mutual escape eccentricity, the surface escape speed of the pair over the local Kepler speed. A uniform draw against $p_\mathrm{col}$ selects the branch:

- **Collision.** The eccentricities are re-drawn from a truncated Rayleigh distribution until the epicycles geometrically overlap, the contact speed is $v_c = \sqrt{v_\infty^2 + v_\mathrm{esc}^2}$, and the pair merges (see below).
- **Scattering.** Both orbits are shifted apart by the excited epicycle amplitudes, weighted with the opposite body's mass fraction (eqs. 18-20), which conserves the mass-weighted sum of semi-major axes exactly.
- **Ejection.** If scattering excites an eccentricity to or past unity, the more excited body escapes[^cite-kimura2025b]. The survivor absorbs the escaper's binding energy, so its orbit tightens to conserve the pair's $M/a$ sum, and it must still pass through the place the encounter happened: that place is now its apocentre, which fixes its eccentricity at $a_\mathrm{old}/a_\mathrm{new} - 1$. For a comparable-mass pair that eccentricity reaches unity, meaning no bound orbit satisfies both conditions, and both bodies leave.

A body whose perihelion falls inside the inner cutoff is removed as having fallen into the star. Scattering can also leave a body on a hyperbolic orbit, which removes it as unbound; in practice this is the channel that fires most often, ahead of both the eccentricity-driven ejection and the infall cutoff.

## Merging

A merger (eqs. 15-17) sums the masses, places the merged body at

$$ a_\mathrm{new} = \frac{M_1 + M_2}{M_1/a_1 + M_2/a_2}, $$

the orbital-energy-weighted orbit, and composes the eccentricity as a mass-weighted vector sum with a pericentre alignment drawn uniformly outside the geometrically forbidden range. The heavier body survives; the lighter one is removed.

Merging is perfect: the surviving body's mass is the plain sum of the two, with no fragmentation, no debris, and nothing shed. A giant impact does erode a real atmosphere, and in a coupled run that erosion is applied downstream by PROTEUS through ZEPHYRUS, from the impact geometry each record carries, following the scaling law of Kegerreis et al. (2020)[^cite-kegerreis2020]. Keeping it out of the dynamical model means one body of physics owns the atmosphere rather than two, and it is what lets the reported merged mass be an exact sum that a consumer can start its own accounting from.

## Randomness and reproducibility

All Monte Carlo draws go through numpy's global random state, seeded once per system from the settings-file seed plus the system index, so a given settings file reproduces its systems exactly. `run_system` additionally saves and restores the caller's random state, so an embedding program's own draws are unaffected.

## References

[^cite-kimura2025]: Kimura, T., Hoshino, H., Kokubo, E., Matsumoto, Y. & Ikoma, M., *[Semi-analytical Model for the Dynamical Evolution of Planetary Systems via Giant Impacts](https://doi.org/10.3847/1538-4357/ade992)*, The Astrophysical Journal, 989, 109, 2025.

[^cite-kimura2025b]: Kimura, T., Kokubo, E., Matsumoto, Y., Mordasini, C. & Ikoma, M., *[Semi-analytical model for the dynamical evolution of planetary system II: Application to systems formed by a planet formation model](https://arxiv.org/abs/2507.03906)*, arXiv:2507.03906, 2025.

[^cite-petit2020]: Petit, A.C., Pichierri, G., Davies, M.B. & Johansen, A., *[The path to instability in compact multi-planetary systems](https://doi.org/10.1051/0004-6361/202038764)*, Astronomy & Astrophysics, 641, A176, 2020.

[^cite-kegerreis2020]: Kegerreis, J.A., Eke, V.R., Catling, D.C., Massey, R.J., Teodoro, L.F.A. & Zahnle, K.J., *[Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density](https://doi.org/10.3847/2041-8213/abb5fb)*, The Astrophysical Journal Letters, 901, L31, 2020.

[^cite-murray1999]: Murray, C.D. & Dermott, S.F., *Solar System Dynamics*, Cambridge University Press, 1999.
