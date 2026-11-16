@testset "community assembly and the community interface" begin
    rng = Xoshiro(2024)
    c, status = nothing, :none
    while c === nothing
        c, status = assemble_community(rng, 12)
    end
    @test status == :ok
    J = jacobian(c)
    @test J == jacobian(c.B, c.model)
    @test maximum(real, eigvals(J)) < 0
    @test all(biomass(c) .> c.model.threshold[c.model.alive])
    @test isempty(isolated_plants(c.model))
    @test richness(c) ≥ floor(Int, 0.9 * c.initial_richness)   # Dataset S1: N ≥ floor(0.9 iniN)
    @test state(c) === c.B && present_species(c) == findall(c.model.alive)

    # ResilienceCore's simulation experiments run on a Community
    opts = PerturbationOptions(sequences = 3)
    TM = global_mortality_tolerance(c; options = opts)
    @test isapprox(round((TM - 1.001) / 0.1), (TM - 1.001) / 0.1; atol = 1e-9)
    loc = local_mortality_tolerances(c; options = opts)
    @test length(loc) == richness(c) && all(>(0), loc)
    rem = removal_responses(c; options = opts)
    @test all(r -> r.survivors ≤ richness(c) - 1, rem)
    @test 1 ≤ extinction_tolerance(rng, c; options = opts) ≤ richness(c)

    # perturbations return a new community and leave the original untouched
    B0, alive0 = copy(c.B), copy(c.model.alive)
    p = with_mortality(c, fill(1.1, length(c.B)))
    r = without_species(c, first(present_species(c)))
    @test p isa Community && r isa Community && richness(r) ≤ richness(c) - 1
    @test c.B == B0 && c.model.alive == alive0
    @test linearise(c) isa ContinuousJacobian
end

@testset "Trophic species of an observed web" begin
    # plant 1; herbivores 2 and 3 both eat only 1; 4 eats both; detritus 5 eaten by 2 and 3
    A = falses(5, 5)
    A[2, 1] = A[3, 1] = A[4, 2] = A[4, 3] = A[2, 5] = A[3, 5] = true
    w = ObservedWeb("toy", "lakes", FoodWeb(A), ["p", "h1", "h2", "c", "d"], [1.0, 4.0, 16.0, 100.0, NaN],
                    BitVector([1, 1, 1, 1, 0]), ["photo", "invertebrate", "invertebrate", "ectotherm", "detritus"])
    t = trophic_species(w)
    @test richness(t.web) == 4
    @test "h1 | h2" in t.taxa
    @test t.mass[findfirst(==("h1 | h2"), t.taxa)] ≈ 8.0
    @test count(t.web.A) == 3                    # H→p, H→d, c→H
    @test trophic_species(t).taxa == t.taxa      # idempotent
end
