# Limitations

Morrigan is a statistical model of the giant-impact phase, not an N-body integrator. Its approximations are inherited from Kimura et al. (2025) and from the choices of this implementation; the ones a user should hold in mind are listed here.

## Dynamical approximations

- **Semi-analytic timing.** Instability times come from the Petit et al. (2020) crossing-time fit and calibrated interaction timescales, not from direct integration. Individual system histories are statistical realisations; only ensemble statistics are meaningful for comparison with N-body results.
- **Perfect merging.** Collisions always merge the pair. There is no fragmentation, no hit-and-run channel, and no debris; the only non-merger mass channel is the atmosphere stripped by the erosion law.
- **Single impact angle.** One configured impact angle applies to every collision of a run; the model does not draw a distribution of impact geometries.
- **Linear secular theory.** Eccentricities between events follow the Laplace-Lagrange solution, which is linear in eccentricity and excludes inclination; the model is two-dimensional, and radially overlapping neighbours are simply decoupled from the secular solution.
- **No dissipation.** There is no gas drag, no tidal damping, and no dynamical friction from a planetesimal population; excitation is removed only by mergers, ejections, and the star-infall cutoff.

## Atmosphere bookkeeping

- The internal atmosphere fractions serve the population statistics of the standalone model. The [erosion law](mass_loss.md) is applied with the two bodies at their bulk densities and its fitted-domain caveats apply; the coupled PROTEUS framework does its own atmosphere accounting through ZEPHYRUS and treats Morrigan's merged masses as perfect-merger sums.
- Bodies are bulk objects: radii do not distinguish a rocky core from an envelope, and the density entering the erosion law is the bulk density.

## Numerical conventions

- The gravitational constant is carried at three significant figures (6.67e-11), the au as 1.5e11 m, and the year as 365 days; results agree with CODATA-constant implementations only to the corresponding relative level, which matters when cross-checking against other codes.
- All Monte Carlo draws share numpy's global random state, seeded per system. Library code must never reseed mid-run; a hidden reseed would silently correlate ensemble members.
- Ejections remove the body excited past unity eccentricity; its unbound orbit is not followed, and the survivor's orbit conserves the pair's orbital energy at the moment of the event.
