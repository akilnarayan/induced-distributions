# induced-distributions
Computation of induced orthogonal polynomial distributions

Evaluating distribution functions
---

The following scripts evaluate distribution functions
- `idist_jacobi.m` : Induced distribution function for Jacobi weights 
- `idist_freud.m`  : Induced distribution function for Freud weights
- `idist_hfreud.m` : Induced distribution function for half-line Freud weights

Evaluating inverse distribution functions
---

The following scripts evaluate inverses of induced distribution functions
- `idistinv_jacobi.m`  : Inverse induced distribution function for Jacobi weights
- `idistinv_freud.m`   : Inverse induced distribution function for Freud weights
- `idistinv_hfreud.m`  : Inverse induced distribution function for half-line Freud weights
- `fidistinv_jacobi.m` : Fast evaluation of `idistinv_jacobi.m`. Requires expensive one-time setup 
- `fidistinv_freud.m`  : Fast evaluation of `idistinv_freud.m`. Requires expensive one-time setup 
- `fidistinv_hfreud.m` : Fast evaluation of `idistinv_hfreud.m`. Requires expensive one-time setup 

Random sampling from mixtures of tensorial induced distributions
---

The script `idist_mixture_sampling.m` samples from an additive mixture of tensorial univariate induced distributions. The method only supports tensorial constructions where each of the univariate marginals are identical.

Demos
---

The following demo files illustrate usage of the methods above.
- `demo_equilibrium_measure.m` : Empirical evidence for multivariate equilibrium measure conjecture. Showcases `idist_mixture_sampling.m`
- `demo_fast_sampling_jacobi.m` : Usage of `idistinv_jacobi.m` versus the faster `fidistinv_jacobi.m`
- `demo_fast_sampling_hfreud.m` : Usage of `idistinv_hfreud.m` versus the faster `fidistinv_hfreud.m`
- `demo_hfreud_error.m` : Error plot for various quadrature sizes for half-line Freud weights. Showcases `idist_hfreud.m`
- `demo_jacobi_error.m` : Error plot for various quadrature sizes for Jacobi weights. Showcases `idist_jacobi.m`
