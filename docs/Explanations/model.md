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
- **Ejection.** If scattering excites an eccentricity to or past unity, the more excited body escapes[^cite-kimura2025b]. The survivor absorbs the escaper's binding energy, so its orbit tightens to conserve the pair's $M/a$ sum, and it must still pass through the place the encounter happened: that place is now its apocentre, which fixes its eccentricity at $a_\mathrm{old}/a_\mathrm{new} - 1$. That closure has a closed form, $e_\mathrm{rem} = (M_\mathrm{ejc}/M_\mathrm{rem})(a_\mathrm{rem,0}/a_\mathrm{ejc,0})$, and it reaches unity when the survivor starts outside the escaper by more than their mass ratio. No bound orbit then satisfies both conditions and both bodies leave. Equipartition always gives the lighter body the larger kick, so the survivor is the heavier one and this outcome needs it to be the outer body as well; a pair with the heavier body on the inside never reaches it, however close their masses. For an exactly equal-mass pair the two excitations are identical and the paper's rule does not say which body goes; the outer one is ejected, which leaves the survivor inside it and therefore bound. Note also what the closure discards: the scattering step has already given the survivor an eccentricity of its own, frequently above unity, and that value is replaced rather than carried forward, so the surviving orbit records the orbit ratio and not how hard the body was hit.

A body whose perihelion falls inside the inner cutoff is removed as having fallen into the star. Scattering can also leave a body on a hyperbolic orbit, which removes it as unbound; in practice this is the channel that fires most often, ahead of both the eccentricity-driven ejection and the infall cutoff.

## Merging

A merger (eqs. 15-17) sums the masses, places the merged body at

$$ a_\mathrm{new} = \frac{M_1 + M_2}{M_1/a_1 + M_2/a_2}, $$

the orbital-energy-weighted orbit, and composes the eccentricity as a mass-weighted vector sum with a pericentre alignment drawn uniformly outside the geometrically forbidden range. The heavier body survives; the lighter one is removed.

Merging is perfect: the surviving body's mass is the plain sum of the two, with no fragmentation, no debris, and nothing shed. A giant impact does erode a real atmosphere, and in a coupled run that erosion is applied downstream by PROTEUS through ZEPHYRUS, from the impact geometry each record carries, following the scaling law of Kegerreis et al. (2020)[^cite-kegerreis2020]. Keeping it out of the dynamical model means one body of physics owns the atmosphere rather than two, and it is what lets the reported merged mass be an exact sum that a consumer can start its own accounting from.

## Mantle melt of a merger

Every merger also reports the bulk melt fraction of the post-impact mantle from the scaling law of Nakajima et al. (2021)[^cite-nakajima2021], a fit to SPH giant impacts between bodies with a 70 % dunite mantle and an iron core. The energy the collision releases is the kinetic energy at contact plus the potential energy of the two bodies minus that of the merged body (their eqs. 2 and 3),

$$ \mathrm{KE}_0 = \tfrac12 \frac{M_t M_i}{M_t + M_i} v_\mathrm{imp}^2, \qquad \Delta \mathrm{PE} = -\tfrac35 \frac{G M_t^2}{R_t} - \tfrac35 \frac{G M_i^2}{R_i} - \frac{G M_t M_i}{R_t + R_i} + \tfrac35 \frac{G (M_t + M_i)^2}{R'}, $$

with Morrigan's own masses, radii and contact speed, and $R'$ the merged radius at the target's density. Legendre fits in $\cos\theta$, with $\sin\theta$ the impact parameter, give the share of that energy that becomes internal energy, $\Delta\mathrm{IE}$ (eq. 4), the mantle's share $h$ of it (eq. 6) and the post-impact mantle mass (eq. 7). The melt fraction (eq. 14) is the mantle's energy gain per unit of the fitted mantle mass, $0.70\,k(\theta)\,M$ with $k$ from 0.924 at grazing to 1.0015 head-on, over $E_M = 5.2 \times 10^6$ J/kg, clipped at 1:

$$ f_\mathrm{melt} = \frac{h\, \Delta \mathrm{IE}}{M_\mathrm{mantle} E_M}. $$

The coefficients are those of the authors' reference code, MeltScalingLaw at commit 2eca30f2, which adds 45° impacts and differs from the paper's Tables S.5 and S.6; the merging coefficients come from the runs at $v_\mathrm{imp} = v_\mathrm{esc}$ only, while the campaign reaches $2 v_\mathrm{esc}$. The reference's sixth-order Legendre term is twice the standard $P_6$, and the port evaluates that basis, the one its coefficients come with. Each merger also reports `f_melt_raw`, the fit value before the clip, which keeps the melt information when $f_\mathrm{melt}$ saturates at 1. The 1-sigma, `f_melt_sigma`, belongs to `f_melt_raw`; it comes from the fit errors at the simulated angle nearest to $\theta$ (0, 30, 45, 60 or 90°, ties to the lower one); the reference code accepts only those five angles, so this rule is the port's. Morrigan uses one impact angle for every merger of a run (`impact_angle` in the settings, 0 to 90°), so a run shows the melt's dependence on mass ratio and speed, not on the angle.

Every merger uses the merging coefficients, since Morrigan merges every collision. `above_v_cr` marks a merger faster than the critical speed of Genda et al. (2012)[^cite-genda2012], where the reference code would switch to its hit-and-run fit. `melt_in_range` is false outside the fitted range, 1 to 53 Mars masses in total, an impactor-to-total mass ratio of 0.03 to 0.5, and 1 to 2 times the mutual escape speed; the value there is the fit extrapolated, and a run with any such merger warns once, naming its seed and system. Python's default warning filter shows an identical warning text once per process, so two runs that differ only in settings the text does not name warn once together; the per-merger fields `melt_in_range` and `above_v_cr` carry the full information. A merger below $v_\mathrm{cr}$ and inside that range still lies beyond the $v_\mathrm{imp} = v_\mathrm{esc}$ runs of the merging fit when it is faster than $v_\mathrm{esc}$, and carries no flag for it. Mergers between Earth-mass embryos at 1.5 to 2 times the escape speed can release enough energy to melt the whole mantle, and $f_\mathrm{melt}$ is then 1.

## Randomness and reproducibility

All Monte Carlo draws go through numpy's global random state, seeded once per system from the settings-file seed plus the system index, so a given settings file reproduces its systems exactly. `run_system` additionally saves and restores the caller's random state, so an embedding program's own draws are unaffected.

## References

[^cite-kimura2025]: Kimura, T., Hoshino, H., Kokubo, E., Matsumoto, Y. & Ikoma, M., *[Semi-analytical Model for the Dynamical Evolution of Planetary Systems via Giant Impacts](https://doi.org/10.3847/1538-4357/ade992)*, The Astrophysical Journal, 989, 109, 2025.

[^cite-kimura2025b]: Kimura, T., Kokubo, E., Matsumoto, Y., Mordasini, C. & Ikoma, M., *[Semi-analytical model for the dynamical evolution of planetary system II: Application to systems formed by a planet formation model](https://arxiv.org/abs/2507.03906)*, arXiv:2507.03906, 2025.

[^cite-petit2020]: Petit, A.C., Pichierri, G., Davies, M.B. & Johansen, A., *[The path to instability in compact multi-planetary systems](https://doi.org/10.1051/0004-6361/202038764)*, Astronomy & Astrophysics, 641, A176, 2020.

[^cite-kegerreis2020]: Kegerreis, J.A., Eke, V.R., Catling, D.C., Massey, R.J., Teodoro, L.F.A. & Zahnle, K.J., *[Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density](https://doi.org/10.3847/2041-8213/abb5fb)*, The Astrophysical Journal Letters, 901, L31, 2020.

[^cite-murray1999]: Murray, C.D. & Dermott, S.F., *Solar System Dynamics*, Cambridge University Press, 1999.

[^cite-nakajima2021]: Nakajima, M., Golabek, G.J., Wünnemann, K., Rubie, D.C., Burger, C., Melosh, H.J., Jacobson, S.A., Manske, L. & Hull, S.D., *[Scaling laws for the geometry of an impact-induced magma ocean](https://doi.org/10.1016/j.epsl.2021.116983)*, Earth and Planetary Science Letters, 568, 116983, 2021.

[^cite-genda2012]: Genda, H., Kokubo, E. & Ida, S., *Merging criteria for giant impacts of protoplanets*, The Astrophysical Journal, 744, 137, 2012.
