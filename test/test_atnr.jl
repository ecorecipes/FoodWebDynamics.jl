# Against ATNr 1.1.2 (CRAN): the same web, masses and states give the same derivatives and
# trajectories. The reference files are written by reference/atnr/export_reference.R. ATNr
# stores consumer–resource arrays as species × consumers (rows resources); FoodWebDynamics uses
# [consumer, resource].
using DelimitedFiles: readdlm

const ATNR = joinpath(@__DIR__, "reference", "atnr")
rdm(p...) = readdlm(joinpath(ATNR, p...), ',', Float64)

function atnr_trajectory(m, u0, times)
    prob = to_ode_problem(m, u0, (first(times), last(times)))
    sol = FoodWebDynamics.SciMLBase.solve(prob, FoodWebDynamics.Tsit5(); saveat = times, reltol = 1e-11,
                                          abstol = 1e-13)
    return reduce(hcat, sol.u)'
end

@testset "ATN models agree with ATNr" begin
    BM = vec(rdm("BM.csv"))
    fw = rdm("fw.csv")
    S, nb = length(BM), 3
    web = FoodWeb(permutedims(fw .> 0))

    # the allometric L-matrix (ATNr's create_Lmatrix, γ = 2, threshold 0.01)
    lweb, L = lmatrix_web(BM; nbasal = nb, γ = 2, threshold = 0.01)
    @test lweb.A == web.A
    @test L .* lweb.A ≈ permutedims(rdm("L.csv")) rtol = 1e-12

    # scaled and unscaled: our constructors, with ATNr's defaults, against ATNr's own models
    for (name, model) in (("scaled", scaled_atn(web; mass = BM)), ("unscaled", unscaled_atn(web; mass = BM)))
        ode = rdm(name, "ode.csv")
        @test atn_rhs!(zeros(S), ode[:, 1], model, 0.0) ≈ ode[:, 2] rtol = 1e-10
        traj = rdm(name, "traj.csv")
        @test atnr_trajectory(model, ode[:, 1], traj[:, 1]) ≈ traj[:, 2:end] rtol = 1e-6
    end

    # nutrients: random parameters, so the model is built from ATNr's draws
    p(n) = rdm("nutrients", "$n.csv")
    b, h, w = p("b"), p("h"), p("w")                       # species × consumers
    num, den = zeros(S, S), zeros(S, S)
    for c in (nb + 1):S, r in 1:S
        fw[r, c] > 0 || continue
        num[c, r] = w[r, c - nb] * b[r, c - nb]
        den[c, r] = w[r, c - nb] * h[r, c - nb] * b[r, c - nb]
    end
    nn = 2
    K, V = ones(nn, S), zeros(nn, S)
    K[:, 1:nb], V[:, 1:nb] = p("K"), p("V")
    growth = NutrientGrowth([vec(p("r")); zeros(S - nb)], K, V, vec(p("S")), only(p("D")))
    model = ATNModel(web, BM, growth; metabolism = vec(p("X")), rate = 1 ./ BM, numerator = num, denominator = den,
                     half_saturation = 1.0, hill = [ones(nb); vec(p("q"))], interference = [zeros(nb); vec(p("c"))],
                     efficiency = vec(p("e")))
    ode = p("ode")
    ours(y) = [y[(nn + 1):end]; y[1:nn]]                  # ATNr: nutrients first; ours: species first
    theirs(y) = [y[(S + 1):end]; y[1:S]]
    @test theirs(atn_rhs!(zeros(S + nn), ours(ode[:, 1]), model, 0.0)) ≈ ode[:, 2] rtol = 1e-10
    traj = p("traj")
    sim = atnr_trajectory(model, ours(ode[:, 1]), traj[:, 1])
    @test reduce(vcat, (theirs(sim[k, :])' for k in axes(sim, 1))) ≈ traj[:, 2:end] rtol = 1e-6
end
