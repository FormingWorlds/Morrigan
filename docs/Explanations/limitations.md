# Limitations

Morrigan is a statistical model of the giant-impact phase, not an N-body integrator. Its approximations are inherited from Kimura et al. (2025)[^cite-kimura2025] and from the choices of this implementation; the ones a user should hold in mind are listed here.

## Dynamical approximations

- **Semi-analytic timing.** Instability times come from the Petit et al. (2020)[^cite-petit2020] crossing-time fit and calibrated interaction timescales, not from direct integration. Individual system histories are statistical realisations; only ensemble statistics are meaningful for comparison with N-body results.
- **Perfect merging.** Collisions always merge the pair. There is no fragmentation, no hit-and-run channel, and no debris, so a body only ever grows: the merged mass is the exact sum of the two, and mass leaves a system only when a body is ejected or falls into the star.
- **Single impact angle.** One configured impact angle applies to every collision of a run; the model does not draw a distribution of impact geometries.
- **Linear secular theory.** Eccentricities between events follow the Laplace-Lagrange solution, which is linear in eccentricity and excludes inclination; the model is two-dimensional, and radially overlapping neighbours are simply decoupled from the secular solution.
- **No dissipation.** There is no gas drag, no tidal damping, and no dynamical friction from a planetesimal population; excitation is removed only by mergers, ejections, and the star-infall cutoff.

## Atmospheres

- **The model carries no atmosphere.** Bodies are bare: a mass is a bulk mass, a radius follows from it at the configured bulk density, and no envelope is tracked or shed. A giant impact erodes a real atmosphere, and that erosion is applied downstream by the consumer rather than here; in a coupled PROTEUS run it comes from `zephyrus.collision.mass_loss`, evaluated from the impact geometry each record carries. The scaling law itself is Kegerreis et al. (2020)[^cite-kegerreis2020].
- Bodies are bulk objects in the dynamics too: radii do not distinguish a rocky core from an envelope, so a body with a heavy envelope is represented by its bulk mass and radius alone.

## Numerical conventions

- The gravitational constant is carried at three significant figures (6.67e-11), the au as 1.5e11 m (0.27 % from the defined value), and the year as 365 days; results agree with CODATA-constant implementations only to the corresponding relative level, which matters when cross-checking against other codes.
- All Monte Carlo draws share numpy's global random state, seeded per system. Library code must never reseed mid-run; a hidden reseed would silently correlate ensemble members.
- Ejections remove the body excited past unity eccentricity; its unbound orbit is not followed, and the survivor's orbit conserves the pair's orbital energy at the moment of the event.

## References

[^cite-kimura2025]: Kimura, T., Hoshino, H., Kokubo, E., Matsumoto, Y. & Ikoma, M., *[Semi-analytical Model for the Dynamical Evolution of Planetary Systems via Giant Impacts](https://doi.org/10.3847/1538-4357/ade992)*, The Astrophysical Journal, 989, 109, 2025.

[^cite-petit2020]: Petit, A.C., Pichierri, G., Davies, M.B. & Johansen, A., *[The path to instability in compact multi-planetary systems](https://doi.org/10.1051/0004-6361/202038764)*, Astronomy & Astrophysics, 641, A176, 2020.

[^cite-kegerreis2020]: Kegerreis, J.A., Eke, V.R., Catling, D.C., Massey, R.J., Teodoro, L.F.A. & Zahnle, K.J., *[Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density](https://doi.org/10.3847/2041-8213/abb5fb)*, The Astrophysical Journal Letters, 901, L31, 2020.
