# Adaptive foraging and pollutant bioaccumulation as an optional layer over any ATN model.
#
# Adaptive foraging (Kondoh 2003; Valdovinos et al. 2010): each consumer c splits its foraging
# effort α_cr over its prey (Σ_r α_cr = 1), and moves effort towards the prey that pays most,
#     dα_cr/dt = G α_cr (g_cr − Σ_k α_ck g_ck),
# with g_cr the consumer's per-capita gain per unit effort on r. Efforts replace the fixed diet
# weights of the base model: the link coefficients are scaled by α_cr n_c (n_c resources), which
# is 1 at equal efforts, so fixed equal efforts reproduce the base model exactly.
#
# Pollutant bioaccumulation (toxicokinetics in the spirit of Garay-Narváez et al. 2013): an
# environmental concentration E receives a release pulse and decays at rate ψ. Species take the
# pollutant up from the environment in proportion to their biomass (ω_i E B_i), consumers take
# it up with their prey (a fraction f of what is ingested, at the prey's concentration
# C_r = P_r/B_r), and the body burden P_i is lost by depuration (ρ_i) and with dying biomass.
# The pollutant is toxic: it adds mortality γ_i C_i, i.e. a loss γ_i P_i of biomass. The
# environment is not depleted by uptake.

"""
    Pollutant(; uptake, depuration, sensitivity, assimilation = 1, decay, load, peak, width)

Parameters of a pollutant (vectors over species where named so):

- `uptake` (ω): uptake from the environment per unit biomass and concentration;
- `depuration` (ρ): elimination rate of the body burden;
- `sensitivity` (γ): extra mortality per unit internal concentration;
- `assimilation` (f): fraction of the pollutant in ingested prey that is retained;
- `decay` (ψ): loss rate of the environmental concentration;
- `load`, `peak`, `width`: a Gaussian release pulse of total amount `load` centred at `peak`.

See [`allometric_pollutant`](@ref) for mass-scaled defaults.
"""
Base.@kwdef struct Pollutant
    uptake::Vector{Float64}
    depuration::Vector{Float64}
    sensitivity::Vector{Float64}
    assimilation::Float64 = 1.0
    decay::Float64
    load::Float64
    peak::Float64
    width::Float64
end

"""
    release(p::Pollutant, t)

The release rate of the pollutant into the environment at time `t`.
"""
release(p::Pollutant, t) = p.load * exp(-(t - p.peak)^2 / (2p.width^2)) / (p.width * sqrt(2π))

"""
    allometric_pollutant(model; ω0 = 0.2, ρ0 = 0.005, γ = 1, assimilation = 1, decay = 1e-3,
                         load = 1, peak = 500, width = 50)

A [`Pollutant`](@ref) whose uptake and depuration rates scale with body mass as `m^−0.25` relative
to the lightest species (`ω = ω0 (m/m₁)^−0.25`, `ρ = ρ0 (m/m₁)^−0.25`), with the same sensitivity
`γ` for every species. The defaults are illustrative, chosen so that a pulse is absorbed over a
few hundred time units of a scaled ATN; they are not fitted to any chemical.
"""
function allometric_pollutant(model::ATNModel; ω0::Real = 0.2, ρ0::Real = 0.005, γ::Real = 1.0,
                              assimilation::Real = 1.0, decay::Real = 1e-3, load::Real = 1.0, peak::Real = 500.0,
                              width::Real = 50.0)
    s = (model.mass ./ minimum(model.mass)) .^ -0.25
    return Pollutant(; uptake = ω0 .* s, depuration = ρ0 .* s, sensitivity = fill(Float64(γ), length(s)),
                     assimilation = assimilation, decay = decay, load = load, peak = peak, width = width)
end

"""
    ExtendedATN(model::ATNModel; adaptation = 0, pollutant = nothing, avoidance = false)

An ATN model with adaptive foraging (speed `adaptation`, the `G` of the replicator equation; 0
keeps diets fixed) and/or a [`Pollutant`](@ref). Foraging effort follows energy gain per unit
effort, `g = e_r × gain`. With `avoidance = true` it also counts the toxic cost of the pollutant
in the prey: `g = gain × (e_r − f C_r γ_c/ρ_c)`. Here `f C_r` is the pollutant ingested per unit
food, and `γ_c/ρ_c` the mortality a unit of internal pollutant causes over its residence time
`1/ρ_c`. This is a modelling choice made here, not taken from a publication. The state is `[B; N; α; E; P]`:

- `B`: species biomasses;
- `N`: nutrients, for [`NutrientGrowth`](@ref);
- `α`: foraging efforts, one per link in the order of `links`, present when `adaptation > 0`;
- `E`, `P`: the environmental concentration and the body burdens, present with a pollutant.

Build initial states with [`initial_state`](@ref).
"""
struct ExtendedATN{G}
    base::ATNModel{G}
    links::Vector{Tuple{Int,Int}}
    adaptation::Float64
    pollutant::Union{Nothing,Pollutant}
    avoidance::Bool
end
function ExtendedATN(model::ATNModel; adaptation::Real = 0.0, pollutant::Union{Nothing,Pollutant} = nothing,
                     avoidance::Bool = false)
    links = [(c, r) for c in eachindex(model.mass) for r in model.resources[c]]
    return ExtendedATN(model, links, Float64(adaptation), pollutant, avoidance)
end

adaptive(x::ExtendedATN) = x.adaptation > 0
function offsets(x::ExtendedATN)
    S, nn = length(x.base.mass), nnutrients(x.base.growth)
    ia = S + nn
    ie = ia + (adaptive(x) ? length(x.links) : 0)
    return (S, nn, ia, ie)
end

"""
    initial_state(x::ExtendedATN, B; N = nothing)

The state `[B; N; α; E; P]` with equal foraging efforts, no pollutant, and nutrients `N` (for a
nutrient model; defaults to the supply concentrations).
"""
function initial_state(x::ExtendedATN, B::AbstractVector; N = nothing)
    m = x.base
    u = Float64[B...]
    nnutrients(m.growth) > 0 && append!(u, something(N, m.growth.supply))
    adaptive(x) && append!(u, [1 / length(m.resources[c]) for (c, _) in x.links])
    x.pollutant === nothing || append!(u, zeros(1 + length(m.mass)))
    return u
end

"""
    extended_rhs!(du, u, x::ExtendedATN, t)

The vector field of an [`ExtendedATN`](@ref).
"""
function extended_rhs!(du, u, x::ExtendedATN, t)
    m = x.base
    S, nn, ia, ie = offsets(x)
    z = zero(eltype(u))
    fill!(du, zero(eltype(du)))
    pol = x.pollutant
    conc(i) = pol === nothing || u[i] ≤ 0 ? z : max(u[ie + 1 + i], z) / u[i]
    l = 0
    for c in 1:S
        res = m.resources[c]
        n = length(res)
        if !(m.alive[c] && n > 0)
            l += n
            continue
        end
        Bq(r) = m.alive[r] ? max(u[r], z)^m.hill[c] : z
        φ(k) = adaptive(x) ? max(u[ia + l + k], z) * n : one(z)
        den = m.half_saturation[c] + m.interference[c] * max(u[c], z)
        for (k, r) in enumerate(res)
            den += φ(k) * m.denominator[c, r] * Bq(r)
        end
        avg = z
        for (k, r) in enumerate(res)
            perunit = m.rate[c] * n * m.numerator[c, r] * Bq(r) / den      # gain per unit effort and biomass
            flux = φ(k) / n * perunit * u[c]
            du[c] += m.efficiency[r] * flux
            du[r] -= flux
            if pol !== nothing
                du[ie + 1 + c] += pol.assimilation * flux * conc(r)
                du[ie + 1 + r] -= flux * conc(r)
            end
            if adaptive(x)
                g = m.efficiency[r] * perunit
                if x.avoidance && pol !== nothing
                    g -= perunit * pol.assimilation * conc(r) * pol.sensitivity[c] / max(pol.depuration[c], eps())
                end
                du[ia + l + k] = g                                      # completed below
                avg += max(u[ia + l + k], z) * g
            end
        end
        if adaptive(x)
            for k in 1:n
                du[ia + l + k] = x.adaptation * max(u[ia + l + k], z) * (du[ia + l + k] - avg)
            end
        end
        l += n
    end
    plant_growth!(du, u, m, m.growth)
    for i in 1:S
        mort = m.mortality[i] * m.base_mortality[i]
        if m.alive[i]
            du[i] -= (m.metabolism[i] + mort) * u[i]
            if pol !== nothing
                du[i] -= pol.sensitivity[i] * max(u[ie + 1 + i], z)
                du[ie + 1 + i] += pol.uptake[i] * max(u[ie + 1], z) * u[i] -
                                  (pol.depuration[i] + mort) * u[ie + 1 + i]
            end
        else
            du[i] = zero(eltype(du))
            pol === nothing || (du[ie + 1 + i] = zero(eltype(du)))
        end
    end
    pol === nothing || (du[ie + 1] = release(pol, t) - pol.decay * u[ie + 1])
    return du
end

to_ode_problem(x::ExtendedATN, u0::AbstractVector, tspan) =
    ODEProblem(ODEFunction((du, u, p, t) -> extended_rhs!(du, u, p, t)), Vector{Float64}(u0), tspan, x)

"""
    simulate(x::ExtendedATN, u0, tspan; saveat = [], reltol = 1e-8, abstol = 1e-10)

Integrate an extended ATN with extinctions: a species whose biomass falls below its threshold is
removed (with its body burden). Works on a copy of the model and returns `(solution, model)`;
`model.base.alive` gives the survivors.
"""
function simulate(x::ExtendedATN, u0::AbstractVector, tspan; saveat = Float64[], reltol::Real = 1e-8,
                  abstol::Real = 1e-10)
    y = ExtendedATN(copy(x.base), x.links, x.adaptation, x.pollutant, x.avoidance)
    m = y.base
    S, _, _, ie = offsets(y)
    dying(u) = any(i -> m.alive[i] && u[i] < m.threshold[i], 1:S)
    function extinguish!(integ)
        for i in 1:S
            if m.alive[i] && integ.u[i] < m.threshold[i]
                m.alive[i] = false
                integ.u[i] = 0.0
                y.pollutant === nothing || (integ.u[ie + 1 + i] = 0.0)
            end
        end
    end
    stops = y.pollutant === nothing ? Float64[] : [y.pollutant.peak]
    sol = SciMLBase.solve(to_ode_problem(y, u0, tspan), Tsit5();
                          callback = DiscreteCallback((u, t, i) -> dying(u), extinguish!), saveat = saveat,
                          reltol = reltol, abstol = abstol, tstops = stops)
    return sol, y
end
