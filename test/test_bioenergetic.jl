@testset "bioenergetic model" begin
    rng = Xoshiro(5)
    p = BioenergeticParameters()
    web = niche_model(rng, 25, 0.15)
    while trophic_levels(web) === nothing
        web = niche_model(rng, 25, 0.15)
    end
    model = FoodWebModel(web, p, draw_exponents(rng, p))
    @test all(model.mass[model.basal] .== 1)
    @test all(model.mass[.!model.basal] .> 1)

    # analytic Jacobian against central differences
    B = 0.1 .+ rand(rng, 25)
    J = zeros(25, 25)
    jacobian!(J, B, model)
    Jfd = similar(J)
    f1, f2 = zeros(25), zeros(25)
    for l in 1:25
        Bp, Bm = copy(B), copy(B)
        Bp[l] += 1e-6
        Bm[l] -= 1e-6
        rhs!(f1, Bp, model, 0.0)
        rhs!(f2, Bm, model, 0.0)
        Jfd[:, l] = (f1 - f2) / 2e-6
    end
    @test maximum(abs, J - Jfd) < 1e-7

    # a lone plant settles at K(1 − d) by default, and at K(1 − x − d) if plants pay metabolism
    A = falses(1, 1)
    lone = FoodWebModel(FoodWeb(A), p, draw_exponents(rng, p))
    Beq, ok = settle!(lone, [0.3])
    @test ok
    @test Beq[1] ≈ p.K * (1 - p.x_plant * p.d0)
    pm = BioenergeticParameters(plant_metabolism = true)
    lone = FoodWebModel(FoodWeb(A), pm, draw_exponents(rng, pm))
    Beq, _ = settle!(lone, [0.3])
    @test Beq[1] ≈ pm.K * (1 - pm.x_plant * (1 + pm.d0))

    # removing a consumer's only resource starves it
    A = falses(2, 2)
    A[2, 1] = true
    pair = FoodWebModel(FoodWeb(A), p, draw_exponents(rng, p))
    Bk = [0.5, 0.2]
    kill!(pair, Bk, 1)
    @test pair.w[2] == 1          # the diet weight is kept (renormalise_diet = false) ...
    Bs, _ = settle!(pair, Bk)
    @test richness(pair) == 0     # ... but nothing is left to eat

    # a consumer of two plants: its weight stays 1/2 when one dies, unless diets are renormalised
    A = falses(3, 3)
    A[3, 1] = A[3, 2] = true
    for (rn, w) in ((false, 0.5), (true, 1.0))
        q = BioenergeticParameters(renormalise_diet = rn)
        trio = FoodWebModel(FoodWeb(A), q, draw_exponents(rng, q))
        @test trio.w[3] == 0.5
        kill!(trio, [0.5, 0.5, 0.2], 1)
        @test trio.w[3] == w
    end
end
