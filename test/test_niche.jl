@testset "niche model" begin
    rng = Xoshiro(11)
    webs = [niche_model(rng, 40, 0.15) for _ in 1:400]
    @test isapprox(sum(connectance, webs) / length(webs), 0.15; atol = 0.01)
    @test all(w -> any(isbasal(w)), webs)
    # species are in niche order, so species 1 (smallest niche value, feeding range 0) is basal
    @test all(w -> isbasal(w)[1], webs)
    @test_throws ArgumentError niche_model(rng, 10, 0.6)

    # a chain 3 → 2 → 1 plus an omnivore 4 eating 1 and 3
    A = falses(4, 4)
    A[2, 1] = A[3, 2] = A[4, 1] = A[4, 3] = true
    L = trophic_levels(FoodWeb(A))
    @test L ≈ [1, 2, 3, 3]
    # a consumer that eats only itself has no path to a basal species
    B = falses(2, 2)
    B[2, 2] = true
    @test trophic_levels(FoodWeb(B)) === nothing
end
