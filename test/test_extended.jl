@testset "adaptive foraging and pollutants" begin
    mass = [1.0, 2.0, 100.0, 1e4]
    web, _ = lmatrix_web(mass; nbasal = 2)          # consumer 3 eats plants 1 and 2; 4 eats 3
    S = 4
    B = [0.6, 0.9, 0.4, 0.2]

    # fixed diets reproduce the base model, for logistic and nutrient growth
    for m in (scaled_atn(web; mass), nutrient_atn(Xoshiro(3), web; mass))
        x = ExtendedATN(m)
        u = initial_state(x, B)
        @test extended_rhs!(zeros(length(u)), u, x, 0.0) ≈ atn_rhs!(zeros(length(u)), u, m, 0.0)
    end

    # adaptive foraging: equal efforts give the base model; efforts stay on the simplex and move
    # towards the more profitable prey
    m = scaled_atn(web; mass)
    x = ExtendedATN(m; adaptation = 0.5)
    u = initial_state(x, B)
    du = extended_rhs!(zeros(length(u)), u, x, 0.0)
    @test du[1:S] ≈ atn_rhs!(zeros(S), B, m, 0.0)
    _, _, ia, _ = FoodWebDynamics.offsets(x)
    k3 = findall(l -> first(l) == 3, x.links)        # consumer 3's two links
    @test sum(du[ia .+ k3]) ≈ 0 atol = 1e-14
    richer = argmax([B[last(x.links[k])] for k in k3])
    @test du[ia + k3[richer]] > 0

    # pollutant bookkeeping: with full retention, no depuration, no decay and no mortality, the
    # total body burden changes only by uptake from the environment
    p = Pollutant(; uptake = fill(0.2, S), depuration = zeros(S), sensitivity = zeros(S), assimilation = 1.0,
                  decay = 0.0, load = 0.0, peak = 0.0, width = 1.0)
    xp = ExtendedATN(m; pollutant = p)
    up = initial_state(xp, B)
    _, _, _, ie = FoodWebDynamics.offsets(xp)
    up[ie + 1] = 1.0                                   # environmental concentration
    up[(ie + 2):end] .= 0.1 .* B                       # body burdens
    dup = extended_rhs!(zeros(length(up)), up, xp, 0.0)
    @test sum(dup[(ie + 2):end]) ≈ sum(0.2 * 1.0 .* B)

    # only plants take it up from the environment: consumers acquire it through their prey, and
    # toxicity lowers biomass
    function exposed(γ)
        q = Pollutant(; uptake = [0.2, 0.2, 0.0, 0.0], depuration = fill(0.005, S), sensitivity = fill(γ, S),
                      decay = 0.0, load = 0.0, peak = 0.0, width = 1.0)
        y = ExtendedATN(m; pollutant = q)
        u0 = initial_state(y, B)
        u0[ie + 1] = 1.0
        sol, _ = simulate(y, u0, (0.0, 200.0); saveat = [200.0])
        return sol.u[end]
    end
    clean, toxic = exposed(0.0), exposed(0.5)
    @test clean[ie + 1 + 3] > 0 && clean[ie + 1 + 4] > 0
    @test sum(toxic[1:S]) < sum(clean[1:S])
end
