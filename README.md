# FoodWebDynamics.jl

Food-web structure and dynamics in Julia:

- **Structure:** niche-model webs (Williams & Martinez 2000), and observed webs from the
  GATEWAy database, which can be aggregated into trophic species.
- **Dynamics:** the allometric bioenergetic model, calibrated to reproduce Domínguez-García,
  Dakos & Kéfi (2019, *PNAS* 116:25714), with an analytic Jacobian.
- **Assembly:** community assembly as in their SI §8–9.
- **Perturbations:** species removal and mortality presses, each followed by a return to
  steady state.

```julia
using FoodWebDynamics
using Random: Xoshiro

rng = Xoshiro(1)
c = nothing
while c === nothing
    global c, status = assemble_community(rng, 30)   # draw, integrate, polish, filter
end
richness(c), biomass(c)
J = jacobian(c)                                      # analytic Jacobian at the equilibrium
```

A `Community` implements the community interface of
[ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl). Its simulation-based
metrics therefore run on it directly: mortality tolerance, press and removal responses, and
robustness to extinctions.
[EcologicalResilience.jl](https://github.com/ecorecipes/EcologicalResilience.jl) builds the
paper's datasets from it through a package extension.

`reaction_network(model)` gives the model as a Catalyst reaction network, after `using Catalyst`.

## Installation

Not yet registered. Its dependencies ResilienceCore and StructuredPopulationCore are not
registered either:

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/ProjectionModels.jl")   # StructuredPopulationCore
Pkg.add(url = "https://github.com/ecorecipes/ResilienceCore.jl")
Pkg.add(url = "https://github.com/ecorecipes/FoodWebDynamics.jl")
```

## Allometric trophic networks

The ATN family is implemented in one general form (`ATNModel`). It covers:
- the scaled model of Delmas et al. (2017), `scaled_atn`;
- the temperature-dependent unscaled model of Binzer et al. (2016), `unscaled_atn`;
- the model with explicit nutrient pools of Schneider et al. (2016), `nutrient_atn`.

All three have interference competition and per-consumer Hill exponents. The
allometric L-matrix web generator is `lmatrix_web`. The implementation was written from the
published equations, not from other packages' code. It is checked against
[ATNr](https://cran.r-project.org/package=ATNr) 1.1.2: on the same web and masses, the vector
fields agree to 10⁻¹⁰ and the trajectories to 10⁻⁶, and the default parametrisations of the
scaled and unscaled models are ATNr's. `test/reference/atnr/export_reference.R` regenerates the
reference.

## Status

Planned: adaptive foraging and pollutant bioaccumulation.
