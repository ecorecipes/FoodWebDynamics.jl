using Test
using FoodWebDynamics
using ResilienceCore: ResilienceCore, PerturbationOptions, global_mortality_tolerance, local_mortality_tolerances,
                      removal_responses, extinction_tolerance, linearise, ContinuousJacobian, isstable
using LinearAlgebra: Diagonal, I, eigvals, norm
using Random: Xoshiro

@testset "FoodWebDynamics" begin
    include("test_niche.jl")
    include("test_bioenergetic.jl")
    include("test_assembly.jl")
    include("test_atn.jl")
    include("test_atnr.jl")
end
