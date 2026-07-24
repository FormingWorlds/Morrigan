# Impact atmosphere loss

Every merger strips a fraction of the colliding pair's combined atmosphere following the scaling law of [Kegerreis et al. (2020, ApJL, 901, L31)](https://iopscience.iop.org/article/10.3847/2041-8213/abb5fb), fitted to smoothed-particle hydrodynamics simulations of giant impacts onto terrestrial planets:

$$ X = 0.64 \left[ \left( \frac{v_c}{v_\mathrm{esc}} \right)^2 \left( \frac{M_i}{M_\mathrm{tot}} \right)^{1/2} \left( \frac{\rho_i}{\rho_t} \right)^{1/2} f(b) \right]^{0.65}, $$

where $v_c$ is the speed at first contact, $v_\mathrm{esc} = \sqrt{2 G (M_t + M_i) / (R_t + R_i)}$ the mutual escape speed, $M_\mathrm{tot} = M_t + M_i$, and $b = \sin \beta$ the impact parameter of an impact at angle $\beta$.

## The geometry factor

Morrigan evaluates the interacting-volume form of the geometry factor,

$$ f(b) = \frac{1}{4} \, \frac{(R_t + R_i)^3}{R_t^3 + R_i^3} \, (1 - b)^2 (1 + 2b), $$

the fraction of the combined volume inside the spherical caps cut by the impact chord. At equal bulk densities, which is how Morrigan's population dynamics uses it, this equals the density-weighted cap-mass fraction of the paper's appendix B. The factor is an exact algebraic zero at fully grazing geometry, $b = 1$.

## Conventions and clamping

- The contact speed convention is $v_c = \sqrt{v_\infty^2 + v_\mathrm{esc}^2}$, so $v_c \geq v_\mathrm{esc}$ always; the law's speed ratio never falls below one.
- The raw power law exceeds one for energetic head-on impacts. `mass_loss` returns the raw value; the merger bookkeeping in `merge_embryo` clamps it into $[0, 1]$ at the call site before applying it to the combined atmosphere.
- The loss applies to the combined atmosphere of both bodies: the impactor's atmosphere is added to the target's before the fraction is taken, and the surviving atmosphere is stored as a fraction of the merged mass.

## Fitted domain

The published fit covers target masses of roughly 0.3 to 3 Earth masses, impactors above 0.05 Earth masses, bulk densities within a factor of about two of Earth's, contact speeds of one to three times the mutual escape speed, and thin atmospheres of order one percent of the planet mass. The reported scatter of the fit is at the ten to twenty percent level, tightest for head-on geometries. Morrigan applies the law across whatever its dynamics produces, so impacts outside this envelope carry the extrapolation uncertainty of the fit.

## Relation to the ZEPHYRUS implementation

The coupled PROTEUS framework applies its own atmosphere loss through `zephyrus.collision.mass_loss`, which implements the appendix-B cap geometry with independent densities for the two bodies. Morrigan's internal copy serves the population statistics of the standalone model. At equal bulk densities the two forms are algebraically identical; the test suite cross-checks them against each other, with the b-dependence agreeing to machine precision and the absolute values to the small residual left by the two packages' gravitational constants.
