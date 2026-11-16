# API Reference

## Food webs

```@docs
FoodWeb
niche_model
richness
connectance
isbasal
trophic_levels
```

## The bioenergetic model

```@docs
BioenergeticParameters
AllometricExponents
draw_exponents
FoodWebModel
kill!
set_mortality!
rhs!
jacobian!
jacobian
reaction_network
```

## Assembly and perturbation

```@docs
SteadyStateOptions
settle!
polish!
consumers
isolated_plants
Community
biomass
assemble_community
PERTURBATION_STEADY
```

## Observed webs

```@docs
ObservedWeb
read_gateway
trophic_species
relative_masses
```

## Allometric trophic networks

The scaled, unscaled and nutrient models of the ATN family in one general form. Their vector
fields and trajectories agree with ATNr 1.1.2; see `test/test_atnr.jl`.

```@docs
ATNModel
PlantGrowth
LogisticGrowth
NutrientGrowth
atn_rhs!
steady_state!
steady_community
ATNCommunity
scaled_atn
unscaled_atn
nutrient_atn
arrhenius
link_strength
lmatrix_web
to_ode_problem
```
