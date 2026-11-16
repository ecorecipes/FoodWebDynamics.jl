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
