# Adaptation-speed sweep for the pollution vignette (01_pollution.qmd).
#
# Every combination of web, strategy, adaptation rate G and pulse load is run once; the
# results are written to adaptation_sweep.csv, which the vignette reads. The full sweep takes
# about an hour with 4 threads, so it is not rerun at render time.
#
#   cd vignettes && julia -t 4 --gcthreads=1 --project=. 01_pollution/adaptation_sweep.jl
#
# Web draws, settling and the pulse follow `persistence` in the vignette exactly.
using FoodWebDynamics, Random, Printf, LinearAlgebra
BLAS.set_num_threads(1)

function draw_web(rng; S = 30, nbasal = 8)
    while true
        mass = sort(10 .^ (4 .* rand(rng, S)))
        web, _ = lmatrix_web(mass; nbasal = nbasal)
        all(c -> any(web.A[c, :]), (nbasal + 1):S) && return web, mass
    end
end

# Settle one web with adaptation rate G, then apply each load. Returns (initial, persisting...).
function sweep_web(seed, G, avoid, loads)
    rng = Xoshiro(seed)
    web, mass = draw_web(rng)
    base = scaled_atn(web; mass)
    x0 = ExtendedATN(base; adaptation = G)
    pre, settled = simulate(x0, initial_state(x0, 0.5 .+ rand(rng, 30)), (0.0, 2000.0); saveat = [2000.0])
    n0 = count(settled.base.alive)
    after_counts = map(loads) do L
        x = ExtendedATN(settled.base; adaptation = G, avoidance = avoid,
                        pollutant = allometric_pollutant(settled.base; load = L))
        u = initial_state(x, pre.u[end][1:30])
        G > 0 && (u[31:(30 + length(x.links))] .= pre.u[end][31:(30 + length(x.links))])
        _, after = simulate(x, u, (0.0, 3000.0); saveat = [3000.0])
        count(after.base.alive)
    end
    return n0, after_counts
end

const LOADS = [0.0, 0.002, 0.005, 0.01, 0.02, 0.05]
const RATES = [0.1, 0.3, 1.0, 3.0, 10.0]
const WEBS = 1:6

configs = vcat([("fixed", 0.0, false)],
               [("adaptive", G, false) for G in RATES],
               [("avoiding", G, true) for G in RATES])
jobs = [(c, w) for c in configs for w in WEBS]
out = Vector{Any}(undef, length(jobs))
t0 = time()
Threads.@threads :dynamic for i in eachindex(jobs)
    ((s, G, avoid), w) = jobs[i]
    t = time()
    out[i] = (s, G, w, sweep_web(w, G, avoid, LOADS)...)
    @printf(stderr, "%-9s G=%5.1f web %d  %6.0fs\n", s, G, w, time() - t)
end
path = joinpath(@__DIR__, "adaptation_sweep.csv")
open(path, "w") do io
    println(io, "web,strategy,G,load,initial,persisting")
    for (s, G, w, n0, counts) in out, (L, n) in zip(LOADS, counts)
        println(io, join((w, s, G, L, n0, n), ","))
    end
end
@printf(stderr, "wrote %s in %.0fs\n", path, time() - t0)
