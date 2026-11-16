"""
    FoodWebDynamics

Food-web structure and dynamics. It covers:

- **Structure:** niche-model webs (Williams & Martinez 2000) and observed GATEWAy webs, which
  can be aggregated into trophic species.
- **Dynamics:** the allometric bioenergetic model, with the calibration that reproduces
  Domínguez-García, Dakos & Kéfi (2019).
- **Assembly:** community assembly as in their SI §8–9.
- **Perturbations:** removing species, pressing mortality and returning to steady state.

A [`Community`](@ref) implements the community interface of ResilienceCore (`state`,
`present_species`, `with_mortality`, `without_species`, `jacobian`), so ResilienceCore's
simulation-based metrics and EcologicalResilience.jl's datasets work on it directly.
"""
module FoodWebDynamics

using CSV: CSV
using DataFrames: DataFrame, groupby
using LinearAlgebra: I, eigvals, issuccess, lu
using OrdinaryDiffEqTsit5: Tsit5
using Random: Random, AbstractRNG
using SciMLBase: SciMLBase, ODEFunction, ODEProblem, DiscreteCallback
using Statistics: mean
using ResilienceCore: ResilienceCore
import ResilienceCore: richness, jacobian, state, present_species, with_mortality, without_species

export FoodWeb, niche_model, richness, connectance, isbasal, trophic_levels
export BioenergeticParameters, AllometricExponents, draw_exponents, FoodWebModel, kill!, set_mortality!, rhs!,
       jacobian!, jacobian
export SteadyStateOptions, settle!, polish!, consumers, isolated_plants, Community, biomass, assemble_community
export state, present_species, with_mortality, without_species, PERTURBATION_STEADY
export ObservedWeb, read_gateway, trophic_species, relative_masses
export reaction_network

include("niche.jl")
include("bioenergetic.jl")
include("assembly.jl")
include("community.jl")
include("gateway.jl")

"""
    reaction_network(model::FoodWebModel, B = nothing; name = :foodweb)

The living part of `model` as a Catalyst `ReactionSystem` (requires `using Catalyst`). There is
one species `Bᵢ` per living population, and the reactions are:

- `∅ → Bᵢ` for plant growth;
- `Bᵢ → ∅` for metabolic and natural losses, with a symbolic mortality multiplier `μᵢ`;
- `Bⱼ → eⱼ B_c` for consumer `c` feeding on `j`, with the full functional response as its rate,
  so the assimilation efficiency is a stoichiometry.

`q` and `K` are symbolic too. Returns `(network, species, mortality, u0, parameters)`, where
`u0` maps species to `B` (if given) and `parameters` maps parameters to the model's values.
Its ODEs are those of [`rhs!`](@ref), and its symbolic Jacobian equals [`jacobian!`](@ref).
"""
function reaction_network end

end
