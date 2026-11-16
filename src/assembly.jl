# Community assembly (SI §9): draw a niche-model web, integrate the bioenergetic model to a
# steady state with extinctions below 10⁻⁶ mᵢ, discard webs that lose more than 10% of their
# species, then remove plants left without consumers.

"""
    SteadyStateOptions(; window = 100.0, tolerance = 5e-4, tmax = 1e6, abstol = 1e-10,
                         reltol = 1e-8, polish = true)

How [`settle!`](@ref) decides that a community is at steady state. The integration proceeds in
windows of length `window`; it stops when every living species' log biomass changed by less
than `tolerance` over the last window (SI §9 uses 0.0005 and leaves the window unstated), or
fails at `tmax` (T_SAFE = 10⁶ in SI §9). With `polish = true` the state is then refined to the
equilibrium by Newton's method. Also with `polish = true`, once the change over a window is
below `early_polish * tolerance`, a Newton step is tried, and accepted if it lands on a locally
stable equilibrium within `early_polish_distance` (in log biomass) of the current state; this
avoids integrating along slow modes for very long times. `alg` is any SciML ODE solver;
implicit solvers are given the model's analytic Jacobian. `maxiters` bounds the number of
solver steps of one integration. A run that has not converged stops once no species has gone
extinct for `stall` time units (default: never); it is then reported as not converged, with
its current state (a point on its trajectory, so that further perturbations can continue
from it).
"""
Base.@kwdef struct SteadyStateOptions{A}
    alg::A = Tsit5()
    window::Float64 = 100.0
    tolerance::Float64 = 5e-4
    tmax::Float64 = 1e6
    abstol::Float64 = 1e-10
    reltol::Float64 = 1e-8
    polish::Bool = true
    early_polish::Float64 = 20.0
    early_polish_distance::Float64 = 0.05
    maxiters::Int = 10^6
    stall::Float64 = Inf
end

function extinction_callback()
    condition(u, t, integ) = any(i -> integ.p.alive[i] && u[i] < integ.p.threshold[i], eachindex(u))
    function affect!(integ)
        model = integ.p
        dead = [i for i in eachindex(integ.u) if model.alive[i] && integ.u[i] < model.threshold[i]]
        kill!(model, integ.u, dead)
        return SciMLBase.derivative_discontinuity!(integ, true)
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

"""
    settle!(model::FoodWebModel, B0; options = SteadyStateOptions())

Integrate `model` from biomasses `B0` until steady state, killing species that fall below their
extinction threshold (which updates `model.alive`). Returns `(B, converged)`, where `B` is the
full biomass vector (zero for dead species) and `converged` is `false` if `options.tmax` was
reached first, the integration failed, or the Newton polish did not land on a locally stable
equilibrium with every survivor above its threshold.
"""
function settle!(model::FoodWebModel, B0::AbstractVector; options::SteadyStateOptions = SteadyStateOptions())
    B = collect(Float64, B0)
    for i in eachindex(B)
        model.alive[i] && B[i] < model.threshold[i] && kill!(model, B, i)
    end
    any(model.alive) || return (B, true)
    f = ODEFunction{true}(rhs!; jac = jacobian!)
    prob = ODEProblem(f, B, (0.0, options.tmax), model)
    integ = SciMLBase.init(prob, options.alg; callback = extinction_callback(), abstol = options.abstol,
                           reltol = options.reltol, save_everystep = false, save_start = false,
                           save_end = false, maxiters = options.maxiters)
    prev = copy(integ.u)
    prev_alive = copy(model.alive)
    converged = false
    last_change = integ.t
    while integ.t < options.tmax
        SciMLBase.step!(integ, min(options.window, options.tmax - integ.t), true)
        SciMLBase.successful_retcode(integ.sol.retcode) || integ.sol.retcode == SciMLBase.ReturnCode.Default ||
            break
        u = integ.u
        if model.alive == prev_alive
            Δ = maximum((abs(log(u[i]) - log(prev[i])) for i in eachindex(u) if model.alive[i]); init = 0.0)
            if Δ < options.tolerance
                converged = true
                break
            end
            # Slow approach to a stable equilibrium: jump there by Newton once the state is
            # nearly steady and the equilibrium is close (avoids crawling along slow modes).
            if options.polish && Δ < options.early_polish * options.tolerance
                Bn = copy(u)
                Bn[.!model.alive] .= 0.0
                if polish!(model, Bn) &&
                   maximum((abs(log(Bn[i]) - log(u[i])) for i in eachindex(u) if model.alive[i]); init = 0.0) <
                   options.early_polish_distance
                    return (Bn, true)
                end
            end
        end
        if model.alive != prev_alive
            last_change = integ.t
        elseif integ.t - last_change ≥ options.stall
            break
        end
        copyto!(prev, u)
        copyto!(prev_alive, model.alive)
    end
    B = copy(integ.u)
    B[.!model.alive] .= 0.0
    if converged && options.polish
        converged = polish!(model, B)
    end
    return (B, converged)
end

"""
    polish!(model::FoodWebModel, B; maxiter = 50, tol = 1e-11)

Refine `B` in place to an equilibrium of the living species by Newton's method. Returns `true`
if it converged to a point where every living species is above its extinction threshold and
the Jacobian is locally stable (every eigenvalue has negative real part).
"""
function polish!(model::FoodWebModel, B::AbstractVector; maxiter::Integer = 50, tol::Real = 1e-11)
    idx = findall(model.alive)
    isempty(idx) && return true
    S = length(B)
    f = zeros(S)
    J = zeros(S, S)
    x = B[idx]
    for _ in 1:maxiter
        rhs!(f, B, model, 0.0)
        r = f[idx]
        maximum(abs, r) < tol && break
        jacobian!(J, B, model)
        Jr = J[idx, idx]
        F = lu(Jr; check = false)
        issuccess(F) || return false
        x .-= F \ r
        B[idx] .= x
    end
    rhs!(f, B, model, 0.0)
    maximum(abs, f[idx]) < 1e3 * tol || return false
    all(i -> B[i] > model.threshold[i], idx) || return false
    return maximum(real, eigvals(jacobian(B, model))) < 0
end

"""
    consumers(model::FoodWebModel, j)

The living consumers of species `j`.
"""
consumers(model::FoodWebModel, j::Integer) =
    [c for c in eachindex(model.alive) if model.alive[c] && j in model.resources[c]]

"""
    isolated_plants(model::FoodWebModel)

Living basal species that no living species eats.
"""
isolated_plants(model::FoodWebModel) =
    [j for j in eachindex(model.alive) if model.alive[j] && model.basal[j] && isempty(consumers(model, j))]

"""
    Community

A dynamically assembled food web at a locally stable equilibrium: the `model` (whose `alive`
marks the survivors), the equilibrium biomasses `B` (zero for dead species), the initial
richness `initial_richness`, and the richness `dynamic_richness` before isolated plants were
removed.
"""
struct Community
    model::FoodWebModel
    B::Vector{Float64}
    initial_richness::Int
    dynamic_richness::Int
end

"""
    richness(c::Community)

The number of species in the assembled community.
"""
richness(c::Community) = richness(c.model)

"""
    biomass(c::Community)

Equilibrium biomasses of the living species, in the order of `findall(c.model.alive)`.
"""
biomass(c::Community) = c.B[c.model.alive]

"""
    assemble_community([rng], S; params = BioenergeticParameters(), C = 0.15,
                       survival = 0.9, initial_biomass = (0.05, 1.0),
                       options = SteadyStateOptions())

One attempt at the assembly procedure of SI §9: draw a niche web of `S` species with
connectance `C` and the community's allometric exponents, start from biomasses uniform on
`initial_biomass`, integrate to steady state and remove plants that no living species eats.
Returns `(community, status)`. `community` is `nothing` and `status` says why when the attempt
is discarded:

- `:trophic_levels`: undefined trophic levels (SI: "the network is discarded");
- `:not_converged`: no steady state by T_SAFE, or no locally stable equilibrium;
- `:extinctions`: fewer than `floor(survival * S)` species are left once isolated plants are
  removed;
- `:empty`: every survivor was an isolated plant, so nothing is left once they are removed
  (possible in small webs without consumers).
- `:neutral`: the equilibrium has an eigenvalue 0 (to within 10⁻⁹), as when two dynamically
  identical consumers survive. It is not asymptotically stable, and `polish!`'s strict `< 0`
  test admits it about half the time because the zero is computed as ±10⁻¹⁷. The published
  Dataset S1 has none (its smallest Rinf is 1.6×10⁻⁴).

The 90% rule is applied as SI §9 orders it, after isolated plants are removed, to the integer
`floor(survival * S)`. Every one of the 19,816 communities of Dataset S1 satisfies
`N ≥ floor(0.9 iniN)` and 5,679 sit exactly at it. The real-valued `N ≥ 0.9 iniN` would reject
those. Applying the rule before the removal would let about 4% of communities fall below the
floor; none of the published ones does.
"""
function assemble_community(rng::AbstractRNG, S::Integer; params::BioenergeticParameters = BioenergeticParameters(),
                            C::Real = 0.15, survival::Real = 0.9, initial_biomass = (0.05, 1.0),
                            options::SteadyStateOptions = SteadyStateOptions())
    web = niche_model(rng, S, C)
    ex = draw_exponents(rng, params)
    trophic_levels(web) === nothing && return (nothing, :trophic_levels)
    model = FoodWebModel(web, params, ex)
    lo, hi = initial_biomass
    B0 = lo .+ (hi - lo) .* rand(rng, S)
    B, converged = settle!(model, B0; options = options)
    converged || return (nothing, :not_converged)
    nd = richness(model)
    iso = isolated_plants(model)
    isempty(iso) || kill!(model, B, iso)
    richness(model) > 0 || return (nothing, :empty)
    richness(model) ≥ floor(Int, survival * S) || return (nothing, :extinctions)
    -maximum(real, eigvals(jacobian(B, model))) > 1e-9 || return (nothing, :neutral)   # return rate −Re λ₁
    return (Community(model, B, S, nd), :ok)
end
assemble_community(S::Integer; kwargs...) = assemble_community(Random.default_rng(), S; kwargs...)
