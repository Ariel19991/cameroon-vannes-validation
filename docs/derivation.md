# Physical derivation of the corrected method

Let the observation point be at radius (r_a=R+h_a), and a terrain cell be at
angular distance \(\psi\) with height (h_p). The original Vannes program
integrates absolute topography above the geoid and applies a condensation term
on the geoid. That residual observable is not, in general, the positive
classical terrain correction relative to the horizontal surface through the
station.

The corrected method uses only the signed residual column between station and
terrain:

\[
r_1=R+\min(h_a,h_p),\qquad r_2=R+\max(h_a,h_p),\qquad
s=\operatorname{sign}(h_p-h_a).
\]

The radial contribution per unit solid angle is evaluated by Gauss–Legendre
quadrature:

\[
I_t=s\int_{r_1}^{r_2}
\frac{r'^2(r'\cos\psi-r_a)}
{(r_a^2+r'^2-2r_ar'\cos\psi)^{3/2}}\,dr'.
\]

The same signed residual volume per solid angle,

\[
V_\Omega=s\frac{r_2^3-r_1^3}{3},
\]

is condensed on the sphere through the station, not on the geoid:

\[
I_c=V_\Omega\,
\frac{r_a(\cos\psi-1)}
{[2r_a^2(1-\cos\psi)]^{3/2}}.
\]

For far cells the correction is

\[
TC_{far}=G\rho\sum_j (I_{t,j}-I_{c,j})\Delta\Omega_j.
\]

Near the observation point, midpoint angular quadrature is not adequate and
the self-cell can be singular. The production estimator therefore uses the
exact Nagy rectangular-prism integral inside a declared radius:

\[
TC=TC_{near}^{Nagy}+TC_{far}^{CRCS}.
\]

This is a hybrid estimator, not an empirical fit. The near-field radius must be
reported and its sensitivity tested. Equality with Nagy in that region is by
construction and cannot be counted as independent evidence.

