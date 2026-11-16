function fd_jacobian(m, u; ε = 1e-6)
    n = length(u)
    J = zeros(n, n)
    for l in 1:n
        h = ε * max(abs(u[l]), 1.0)
        up, um = copy(u), copy(u)
        up[l] += h
        um[l] -= h
        J[:, l] = (atn_rhs!(zeros(n), up, m, 0.0) - atn_rhs!(zeros(n), um, m, 0.0)) / 2h
    end
    return J
end

@testset "ATN: the general form reproduces the calibrated model" begin
    rng = Xoshiro(7)
    c = nothing
    while c === nothing
        c, _ = assemble_community(rng, 20)
    end
    atn = ATNModel(c.model)
    u = c.B .* (1 .+ 0.1 .* rand(rng, length(c.B)))
    du1, du2 = zeros(length(u)), zeros(length(u))
    rhs!(du1, u, c.model, 0.0)
    atn_rhs!(du2, u, atn, 0.0)
    @test du2 ≈ du1 rtol = 1e-10 atol = 1e-14
    idx = findall(c.model.alive)
    @test jacobian(atn, c.B)[idx, idx] ≈ jacobian(c.B, c.model) rtol = 1e-8 atol = 1e-12
end

@testset "ATN: steady states with known solutions" begin
    web = FoodWeb(falses(1, 1))
    # a lone plant in the scaled model (paying metabolism x = 0.314) settles at K (1 − x/r)
    c, ok = steady_community(scaled_atn(web; mass = [1.0]), [0.5])
    @test ok
    @test state(c)[1] ≈ 10.0 * (1 - 0.314) rtol = 1e-2
    # a lone plant on one nutrient: G* = x/r, N* = K x/(r − x), B* = D (S − N*)/(v x)
    r, x, K, v, S, D = 1.0, 0.2, 0.15, 1.0, 10.0, 0.25
    g = NutrientGrowth([r], fill(K, 1, 1), fill(v, 1, 1), [S], D)
    m = ATNModel(web, [1.0], g; metabolism = x, rate = 1.0, numerator = zeros(1, 1), denominator = zeros(1, 1),
                 half_saturation = 1.0, hill = 1.2, efficiency = 0.45)
    cn, okn = steady_community(m, [1.0, 5.0]; tmax = 1e6)
    Nstar = K * x / (r - x)
    @test okn
    @test cn.u[2] ≈ Nstar rtol = 1e-2
    @test cn.u[1] ≈ D * (S - Nstar) / (v * x) rtol = 1e-2
    @test isstable(linearise(cn)) && size(jacobian(cn)) == (2, 2)   # plant and nutrient
end

@testset "ATN: variants, temperature, interference and the L-matrix" begin
    @test arrhenius(0.5, 293.15) == 1.0
    @test arrhenius(-0.84, 303.15) > 1          # faster when warmer, in this sign convention
    mass = [1.0, 1.0, 100.0, 1e4]
    web, L = lmatrix_web(mass; nbasal = 2)
    @test L[3, 1] ≈ 1 && L[4, 3] ≈ 1            # at the optimal mass ratio
    @test web.A[3, 1] && web.A[3, 2] && web.A[4, 3] && !web.A[4, 1] && !any(web.A[1:2, :])
    for (model, nn) in ((scaled_atn(web; mass), 0), (unscaled_atn(web; mass), 0),
                        (nutrient_atn(Xoshiro(1), web; mass), 2))
        u0 = [0.5, 0.8, 0.3, 0.2, fill(5.0, nn)...]
        du = atn_rhs!(zeros(length(u0)), u0, model, 0.0)
        @test all(isfinite, du)
        @test jacobian(model, u0) ≈ fd_jacobian(model, u0) rtol = 1e-5 atol = 1e-12
    end
    # interference lowers what the top predator gains
    u = [0.5, 0.8, 0.3, 0.2]
    plain, crowded = scaled_atn(web; mass), scaled_atn(web; mass, interference = 1.0)
    @test atn_rhs!(zeros(4), u, crowded, 0.0)[4] < atn_rhs!(zeros(4), u, plain, 0.0)[4]

    # the ResilienceCore experiments run on an ATN community
    c, ok = steady_community(scaled_atn(web; mass), [0.5, 0.8, 0.3, 0.2]; tmax = 1e5)
    @test ok && richness(c) ≥ 1
    rem = removal_responses(c)
    @test length(rem) == richness(c) && all(r -> r.survivors ≤ richness(c) - 1, rem)
end
