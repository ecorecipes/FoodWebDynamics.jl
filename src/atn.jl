# The allometric trophic network (ATN) family in one general form.
#
# For consumer c feeding on resource r, the biomass flux is
#
#     Φ_cr = rate_c · B_c · num_cr · B_r^h / (half_c + interference_c · B_c + Σ_k den_ck · B_k^h),
#
# the consumer gains Σ_r e_r Φ_cr and the resource loses Σ_c Φ_cr. Plants grow logistically
# with competition, or on explicit nutrient pools, and every species loses
# (metabolism + mortality multiplier × base mortality) × B. Each published variant is a choice of
# these coefficients:
#
# - unscaled ATN (Binzer et al. 2016): num = a, den = a h, half = 1, rate = 1/m, with
#   temperature-dependent (Arrhenius) rates;
# - scaled ATN (Delmas et al. 2017, as in BioEnergeticFoodWebs): num = den = w, half = B₀^h,
#   rate = x y;
# - nutrient ATN (Schneider et al. 2016): num = w b, den = w h b with allometric link strengths
#   b ∝ L, nutrient-limited plant growth;
# - the calibrated model of Domínguez-García et al. (2019): num = w a, den = w a h, half = 1,
#   rate = 1/m (the same as `FoodWebModel`; see `ATNModel(::FoodWebModel)`).
#
# Defaults follow the description of these models in the documentation of the ATNr package
# (Gauzens et al.); where it leaves a value open, the choice is stated in the docstring.

"""
    PlantGrowth

How plants grow in an [`ATNModel`](@ref): [`LogisticGrowth`](@ref) or [`NutrientGrowth`](@ref).
"""
abstract type PlantGrowth end

"""
    LogisticGrowth(r, K, α)

Logistic plant growth `rᵢ (1 − Σⱼ αᵢⱼ Bⱼ / Kᵢ) Bᵢ` over the living basal species: `α = I` gives
each plant its own resource, `α` of ones a shared resource pool.
"""
struct LogisticGrowth <: PlantGrowth
    r::Vector{Float64}
    K::Vector{Float64}
    α::Matrix{Float64}
end

"""
    NutrientGrowth(r, K, v, supply, turnover)

Nutrient-limited plant growth (Schneider et al. 2016): `rᵢ Gᵢ Bᵢ` with
`Gᵢ = minₗ Nₗ / (Kₗᵢ + Nₗ)`, and nutrient pools `dNₗ/dt = D (Sₗ − Nₗ) − Σᵢ vₗᵢ rᵢ Gᵢ Bᵢ` (`D` the
`turnover`, `S` the `supply`). `K` and `v` are nutrients × species.
"""
struct NutrientGrowth <: PlantGrowth
    r::Vector{Float64}
    K::Matrix{Float64}
    v::Matrix{Float64}
    supply::Vector{Float64}
    turnover::Float64
end

nnutrients(::LogisticGrowth) = 0
nnutrients(g::NutrientGrowth) = length(g.supply)

"""
    ATNModel(web, mass, growth; metabolism, rate, numerator, denominator, half_saturation,
             hill, efficiency, interference = 0, base_mortality = 0, threshold = 1e-6)

An allometric trophic network in the general form above. The state is the vector of species
biomasses, followed by the nutrient concentrations for [`NutrientGrowth`](@ref). Usually built
with [`scaled_atn`](@ref), [`unscaled_atn`](@ref), [`nutrient_atn`](@ref) or from a
[`FoodWebModel`](@ref). Mortality multipliers (`mortality`, all one at baseline) scale
`base_mortality` for press experiments.
"""
mutable struct ATNModel{G<:PlantGrowth}
    web::FoodWeb
    mass::Vector{Float64}
    basal::BitVector
    growth::G
    metabolism::Vector{Float64}
    base_mortality::Vector{Float64}
    mortality::Vector{Float64}
    rate::Vector{Float64}
    numerator::Matrix{Float64}
    denominator::Matrix{Float64}
    half_saturation::Vector{Float64}
    interference::Vector{Float64}
    hill::Vector{Float64}
    efficiency::Vector{Float64}
    alive::BitVector
    threshold::Vector{Float64}
    resources::Vector{Vector{Int}}
end

function ATNModel(web::FoodWeb, mass::AbstractVector, growth::PlantGrowth; metabolism, rate, numerator, denominator,
                  half_saturation, hill, efficiency, interference = 0.0, base_mortality = 0.0,
                  threshold = 1e-6)
    S = richness(web)
    asvec(x) = x isa Real ? fill(Float64(x), S) : Vector{Float64}(x)
    return ATNModel(web, Vector{Float64}(mass), isbasal(web), growth, asvec(metabolism), asvec(base_mortality),
                    ones(S), asvec(rate), Matrix{Float64}(numerator), Matrix{Float64}(denominator),
                    asvec(half_saturation), asvec(interference), asvec(hill), asvec(efficiency), trues(S),
                    asvec(threshold), [findall(web.A[i, :]) for i in 1:S])
end

Base.copy(m::ATNModel) = ATNModel(m.web, m.mass, m.basal, m.growth, m.metabolism, m.base_mortality, copy(m.mortality),
                                  m.rate, m.numerator, m.denominator, m.half_saturation, m.interference, m.hill,
                                  m.efficiency, copy(m.alive), m.threshold, m.resources)

richness(m::ATNModel) = count(m.alive)
nstate(m::ATNModel) = length(m.mass) + nnutrients(m.growth)

"""
    atn_rhs!(du, u, m::ATNModel, t)

The vector field of the ATN model (species biomasses, then nutrients). Negative states are
treated as zero in the rates; dead species do not change.
"""
function atn_rhs!(du, u, m::ATNModel, t)
    S = length(m.mass)
    z = zero(eltype(u))
    fill!(du, zero(eltype(du)))
    for c in 1:S
        (m.alive[c] && !isempty(m.resources[c])) || continue
        Bq(r) = m.alive[r] ? max(u[r], z)^m.hill[c] : z      # Hill exponent of the consumer
        den = m.half_saturation[c] + m.interference[c] * max(u[c], z)
        for r in m.resources[c]
            den += m.denominator[c, r] * Bq(r)
        end
        scale = m.rate[c] * u[c] / den
        for r in m.resources[c]
            flux = scale * m.numerator[c, r] * Bq(r)
            du[c] += m.efficiency[r] * flux
            du[r] -= flux
        end
    end
    plant_growth!(du, u, m, m.growth)
    for i in 1:S
        du[i] = m.alive[i] ? du[i] - (m.metabolism[i] + m.mortality[i] * m.base_mortality[i]) * u[i] : zero(eltype(du))
    end
    return du
end

function plant_growth!(du, u, m::ATNModel, g::LogisticGrowth)
    plants = [i for i in eachindex(m.mass) if m.basal[i] && m.alive[i]]
    for i in plants
        s = sum(g.α[i, j] * u[j] for j in plants)
        du[i] += g.r[i] * (1 - s / g.K[i]) * u[i]
    end
end
function plant_growth!(du, u, m::ATNModel, g::NutrientGrowth)
    S, z = length(m.mass), zero(eltype(u))
    N = [max(u[S + l], z) for l in eachindex(g.supply)]
    for l in eachindex(g.supply)
        du[S + l] = g.turnover * (g.supply[l] - u[S + l])
    end
    for i in eachindex(m.mass)
        (m.basal[i] && m.alive[i]) || continue
        G = minimum(N[l] / (g.K[l, i] + N[l]) for l in eachindex(g.supply))
        uptake = g.r[i] * G * u[i]
        du[i] += uptake
        for l in eachindex(g.supply)
            du[S + l] -= g.v[l, i] * uptake
        end
    end
end

"""
    jacobian(m::ATNModel, u)

The Jacobian of [`atn_rhs!`](@ref) at state `u`, by forward-mode automatic differentiation.
"""
jacobian(m::ATNModel, u::AbstractVector) =
    ForwardDiff.jacobian((du, x) -> atn_rhs!(du, x, m, 0.0), zeros(length(u)), Vector{Float64}(u))

"""
    to_ode_problem(m::ATNModel, u0, tspan)

The ATN model as a SciML `ODEProblem` (the shared lowering verb of the ecorecipes packages).
"""
to_ode_problem(m::ATNModel, u0::AbstractVector, tspan) =
    ODEProblem(ODEFunction((du, u, p, t) -> atn_rhs!(du, u, p, t)), Vector{Float64}(u0), tspan, m)

"""
    steady_state!(m::ATNModel, u0; tmax = 1e6, window = 100.0, tolerance = 5e-4,
                  abstol = 1e-10, reltol = 1e-8)

Integrate from `u0` in windows of length `window` until every living species' log biomass
changes by less than `tolerance` over a window (the criterion of Domínguez-García et al. 2019),
or until `tmax`. A species is extinct when its biomass falls below its `threshold`; it is
removed from `m` (`m.alive`) and set to zero. Returns `(u, converged)`.
"""
function steady_state!(m::ATNModel, u0::AbstractVector; tmax::Real = 1e6, window::Real = 100.0,
                       tolerance::Real = 5e-4, abstol::Real = 1e-10, reltol::Real = 1e-8)
    S = length(m.mass)
    dying(u) = any(i -> m.alive[i] && u[i] < m.threshold[i], 1:S)
    function extinguish!(integ)
        for i in 1:S
            if m.alive[i] && integ.u[i] < m.threshold[i]
                m.alive[i] = false
                integ.u[i] = 0.0
            end
        end
    end
    cb = DiscreteCallback((u, t, integ) -> dying(u), extinguish!)
    u = Vector{Float64}(u0)
    for i in 1:S
        m.alive[i] || (u[i] = 0.0)
    end
    integ = SciMLBase.init(to_ode_problem(m, u, (0.0, Float64(tmax))), Tsit5(); callback = cb, abstol = abstol,
                           reltol = reltol, save_everystep = false)
    prev = copy(integ.u)
    while integ.t < tmax
        SciMLBase.step!(integ, min(window, tmax - integ.t), true)
        SciMLBase.successful_retcode(integ.sol.retcode) || integ.sol.retcode == SciMLBase.ReturnCode.Default ||
            return (copy(integ.u), false)
        cur = integ.u
        live = [i for i in 1:S if m.alive[i]]
        isempty(live) && return (copy(cur), true)
        all(abs(log(cur[i]) - log(prev[i])) < tolerance for i in live if prev[i] > 0) && return (copy(cur), true)
        prev = copy(cur)
    end
    return (copy(integ.u), false)
end

"""
    ATNCommunity(model, u)

An ATN model at a (steady) state `u`. It implements the ResilienceCore community interface, so
the simulation experiments run on it. Its [`jacobian`](@ref) covers the living species and the
nutrient pools, and `linearise` includes the nutrients in the equilibrium.
"""
struct ATNCommunity{G}
    model::ATNModel{G}
    u::Vector{Float64}
end

"""
    steady_community(m::ATNModel, u0; kwargs...)

Integrate to steady state ([`steady_state!`](@ref)) and return `(ATNCommunity, converged)`.
"""
function steady_community(m::ATNModel, u0::AbstractVector; kwargs...)
    model = copy(m)
    u, ok = steady_state!(model, u0; kwargs...)
    return ATNCommunity(model, u), ok
end

state(c::ATNCommunity) = c.u[1:length(c.model.mass)]
present_species(c::ATNCommunity) = findall(c.model.alive)
richness(c::ATNCommunity) = count(c.model.alive)
atn_index(c::ATNCommunity) = [present_species(c); length(c.model.mass) .+ (1:nnutrients(c.model.growth))]
jacobian(c::ATNCommunity) = (idx = atn_index(c); jacobian(c.model, c.u)[idx, idx])
ResilienceCore.linearise(c::ATNCommunity) = ResilienceCore.ContinuousJacobian(jacobian(c), c.u[atn_index(c)])

function with_mortality(c::ATNCommunity, factors::AbstractVector; newton::Bool = true, steady = nothing)
    model = copy(c.model)
    copyto!(model.mortality, factors)
    u, _ = steady_state!(model, c.u; something(steady, (;))...)
    return ATNCommunity(model, u)
end

function without_species(c::ATNCommunity, i; steady = nothing)
    model = copy(c.model)
    model.alive[i] = false
    u = copy(c.u)
    u[i] = 0.0
    u, _ = steady_state!(model, u; something(steady, (;))...)
    return ATNCommunity(model, u)
end

"""
    ATNModel(m::FoodWebModel)

The calibrated bioenergetic model of Domínguez-García et al. (2019) in the general ATN form, with
the diet weights and survivors of `m` at the time of conversion. Its vector field equals
[`rhs!`](@ref) of `m`, which the tests check.
"""
function ATNModel(m::FoodWebModel)
    S = length(m.mass)
    num = zeros(S, S)
    den = zeros(S, S)
    for c in 1:S, r in m.resources[c]
        num[c, r] = m.w[c] * m.attack[c, r]
        den[c, r] = m.w[c] * m.attack[c, r] * m.handling[c, r]
    end
    growth = LogisticGrowth(copy(m.growth), fill(m.params.K, S), Matrix{Float64}(I, S, S))
    model = ATNModel(m.web, m.mass, growth; metabolism = m.metabolism, base_mortality = m.base_mortality,
                     rate = 1 ./ m.mass, numerator = num, denominator = den, half_saturation = 1.0, hill = 1 + m.params.q,
                     efficiency = m.efficiency, threshold = m.threshold)
    copyto!(model.mortality, m.mortality)
    copyto!(model.alive, m.alive)
    return model
end

const BOLTZMANN = 8.617333262e-5   # eV K⁻¹

"""
    arrhenius(E, T; T0 = 293.15)

The temperature factor `exp(E (T₀ − T) / (k T T₀))` of the unscaled ATN models (`k` the Boltzmann
constant in eV K⁻¹), equal to 1 at `T₀`.
"""
arrhenius(E::Real, T::Real; T0::Real = 293.15) = exp(E * (T0 - T) / (BOLTZMANN * T * T0))

basal_masses(web, mass) = (b = isbasal(web); any(b) || throw(ArgumentError("the web has no basal species")); b)

"""
    scaled_atn(web; mass, K = 10, x_ratio = 0.314, y = 8, B0 = 0.5, hill = 1.2, interference = 0.8,
               e_plant = 0.545, e_animal = 0.906, competition = :independent, plant_metabolism = true,
               threshold = 1e-6)

The scaled ATN (Delmas et al. 2017; BioEnergeticFoodWebs). Rates are expressed relative to the
growth rate of the smallest plant: `rᵢ = (mᵢ/m₁)^−0.25` for plants and `xᵢ = x_ratio (mᵢ/m₁)^−0.25`,
with the maximum consumption `y` relative to metabolism. The functional response is
`F = w B^h / (B₀^h + c B_c + Σ w B^h)` with diet weights `w = 1/(number of resources)`.
`competition = :shared` puts all plants on one resource pool (`α` of ones).

The defaults are those ATNr 1.1.2 actually initialises (`initialise_default_Scaled`):
x/r = 0.314, y = 8, K = 10, B₀ = 0.5, Hill exponent 1.2, interference 0.8, efficiencies 0.545
(plant food) and 0.906 (animal food), and plants paying metabolism. The test suite checks the
vector field against ATNr's.
"""
function scaled_atn(web::FoodWeb; mass::AbstractVector, K::Real = 10.0, x_ratio::Real = 0.314, y::Real = 8.0,
                    B0::Real = 0.5, hill::Real = 1.2, interference::Real = 0.8, e_plant::Real = 0.545,
                    e_animal::Real = 0.906, competition::Symbol = :independent, plant_metabolism::Bool = true,
                    threshold::Real = 1e-6)
    S = richness(web)
    basal = basal_masses(web, mass)
    m1 = minimum(mass[basal])
    r = [basal[i] ? (mass[i] / m1)^-0.25 : 0.0 for i in 1:S]
    x = [basal[i] && !plant_metabolism ? 0.0 : x_ratio * (mass[i] / m1)^-0.25 for i in 1:S]
    num = zeros(S, S)
    for c in 1:S
        res = findall(web.A[c, :])
        isempty(res) || (num[c, res] .= 1 / length(res))
    end
    α = competition === :shared ? ones(S, S) : Matrix{Float64}(I, S, S)
    return ATNModel(web, mass, LogisticGrowth(r, fill(Float64(K), S), α); metabolism = x, rate = x .* y,
                    numerator = num, denominator = num, half_saturation = B0^hill, hill = hill,
                    interference = interference, efficiency = [basal[j] ? e_plant : e_animal for j in 1:S],
                    threshold = threshold)
end

"""
    unscaled_atn(web; mass, T = 293.15, k0 = 40, hill = 1.2, interference = 0.8,
                 e_plant = 0.545, e_animal = 0.906, competition = :independent, threshold = 1e-6)

The unscaled, temperature-dependent ATN (Binzer et al. 2016), with rates per second and
Arrhenius factors ([`arrhenius`](@ref)). For consumer `c` and resource `r`:

- attack rate `a = e^−13.1 m_r^0.25 m_c^−0.8 A(−0.38)`;
- handling time `h = e^9.66 m_r^−0.45 m_c^0.47 A(0.26)`;
- carrying capacity `K = k0 m^0.28 A(0.71)`;
- plant growth `r = e^−15.68 m^−0.25 A(−0.84)`;
- metabolism `X = e^−16.54 m^−0.31 A(−0.69)` for every species.

The functional response is `F = a B^h / (1 + c B_c + Σ a h B^h) / m_c`. The values are those ATNr 1.1.2 initialises
(`initialise_default_Unscaled`), including `k0 = 40`, a Hill exponent of 1.2 and efficiencies
0.545 and 0.906. The test suite checks the vector field and trajectories against ATNr's.
"""
function unscaled_atn(web::FoodWeb; mass::AbstractVector, T::Real = 293.15, k0::Real = 40.0, hill::Real = 1.2,
                      interference::Real = 0.8, e_plant::Real = 0.545, e_animal::Real = 0.906,
                      competition::Symbol = :independent, threshold::Real = 1e-6)
    S = richness(web)
    basal = basal_masses(web, mass)
    num, den = zeros(S, S), zeros(S, S)
    for c in 1:S, r in findall(web.A[c, :])
        a = exp(-13.1) * mass[r]^0.25 * mass[c]^-0.8 * arrhenius(-0.38, T)
        h = exp(9.66) * mass[r]^-0.45 * mass[c]^0.47 * arrhenius(0.26, T)
        num[c, r], den[c, r] = a, a * h
    end
    r = [basal[i] ? exp(-15.68) * mass[i]^-0.25 * arrhenius(-0.84, T) : 0.0 for i in 1:S]
    K = [k0 * mass[i]^0.28 * arrhenius(0.71, T) for i in 1:S]
    X = [exp(-16.54) * mass[i]^-0.31 * arrhenius(-0.69, T) for i in 1:S]
    α = competition === :shared ? ones(S, S) : Matrix{Float64}(I, S, S)
    return ATNModel(web, mass, LogisticGrowth(r, K, α); metabolism = X, rate = 1 ./ mass, numerator = num,
                    denominator = den, half_saturation = 1.0, hill = hill, interference = interference,
                    efficiency = [basal[j] ? e_plant : e_animal for j in 1:S], threshold = threshold)
end

"""
    link_strength(m_consumer, m_resource; Ropt = 100, γ = 2)

The allometric link strength `L = ((m_c / (m_r R_opt)) e^{1 − m_c/(m_r R_opt)})^γ`: one at the
optimal consumer–resource mass ratio `R_opt`, falling off on either side (Schneider et al. 2016).
The defaults are those of ATNr's `create_Lmatrix`, against which the tests check it.
"""
function link_strength(mc::Real, mr::Real; Ropt::Real = 100.0, γ::Real = 2.0)
    x = mc / (mr * Ropt)
    return (x * exp(1 - x))^γ
end

"""
    lmatrix_web(mass; nbasal, Ropt = 100, γ = 2, threshold = 0.01)

An allometric food web: species are sorted by body mass, the `nbasal` lightest are plants, and
every other species `c` eats every species `r` whose [`link_strength`](@ref) `L(m_c, m_r)` exceeds
`threshold` (the L-matrix of Schneider et al. 2016). Returns `(web, L)` with `L[c, r]`.
"""
function lmatrix_web(mass::AbstractVector; nbasal::Integer, Ropt::Real = 100.0, γ::Real = 2.0, threshold::Real = 0.01)
    issorted(mass) || throw(ArgumentError("masses must be sorted in increasing order"))
    S = length(mass)
    L = [c > nbasal ? link_strength(mass[c], mass[r]; Ropt = Ropt, γ = γ) : 0.0 for c in 1:S, r in 1:S]
    return FoodWeb(L .> threshold), L
end

"""
    nutrient_atn(rng, web; mass, nutrients = 2, T = 293.15, Ropt = 100, hill = 1.2, turnover = 0.25,
                 e_plant = 0.545, e_animal = 0.906, threshold = 1e-6)

The ATN with explicit nutrient dynamics (Schneider et al. 2016), with random parameters drawn as
in ATNr's default parametrisation:

- **Capture coefficient:** `b = 50 m_r^{b₁} m_c^{b₂} A(−0.38) L`, with `b₁ ~ N(0.15, 0.03)`,
  `b₂ ~ N(0.47, 0.04)` and the link strength `L` ([`link_strength`](@ref)).
- **Handling time:** `h = 0.4 m_r^{h₁} m_c^{h₂} A(0.26)`, with `h₁ ~ N(−0.66, 0.02)` and
  `h₂ ~ N(−0.48, 0.03)`.
- **Interference:** `c ~ N(0.8, 0.2)` per consumer, truncated at 0.
- **Plant growth:** `r = m^−0.25 A(−0.84)`.
- **Metabolism:** `X = x₀ m^−0.25 A(−0.69)`, with `x₀ = 0.138` for plants and 0.314 for animals.
- **Nutrients:** half-saturation `Kₗᵢ ~ U(0.1, 0.2)`, plant content `vₗᵢ ~ U(1, 2)` normalised
  over nutrients, supply `Sₗ ~ N(10, 2)`, turnover `D = 0.25`.

The functional response is `F = w b B^h / (1 + c B_c + Σ w h b B^h) / m_c`. The allometric
exponents are drawn once per community here. ATNr's description does not say whether they are
drawn per community or per link. `R_opt = 100` and the L-matrix's γ = 2 are ATNr's defaults.
"""
function nutrient_atn(rng::AbstractRNG, web::FoodWeb; mass::AbstractVector, nutrients::Integer = 2, T::Real = 293.15,
                      Ropt::Real = 100.0, hill::Real = 1.2, turnover::Real = 0.25, e_plant::Real = 0.545,
                      e_animal::Real = 0.906, threshold::Real = 1e-6)
    S = richness(web)
    basal = basal_masses(web, mass)
    b1, b2 = 0.15 + 0.03randn(rng), 0.47 + 0.04randn(rng)
    h1, h2 = -0.66 + 0.02randn(rng), -0.48 + 0.03randn(rng)
    num, den = zeros(S, S), zeros(S, S)
    for c in 1:S
        res = findall(web.A[c, :])
        isempty(res) && continue
        w = 1 / length(res)
        for r in res
            b = 50 * mass[r]^b1 * mass[c]^b2 * arrhenius(-0.38, T) * link_strength(mass[c], mass[r]; Ropt = Ropt)
            h = 0.4 * mass[r]^h1 * mass[c]^h2 * arrhenius(0.26, T)
            num[c, r], den[c, r] = w * b, w * h * b
        end
    end
    c = [basal[i] ? 0.0 : max(0.8 + 0.2randn(rng), 0.0) for i in 1:S]
    K = 0.1 .+ 0.1 .* rand(rng, nutrients, S)
    v = 1 .+ rand(rng, nutrients, S)
    v ./= sum(v; dims = 1)
    supply = 10 .+ 2 .* randn(rng, nutrients)
    r = [basal[i] ? mass[i]^-0.25 * arrhenius(-0.84, T) : 0.0 for i in 1:S]
    X = [(basal[i] ? 0.138 : 0.314) * mass[i]^-0.25 * arrhenius(-0.69, T) for i in 1:S]
    return ATNModel(web, mass, NutrientGrowth(r, K, v, supply, turnover); metabolism = X, rate = 1 ./ mass,
                    numerator = num, denominator = den, half_saturation = 1.0, hill = hill, interference = c,
                    efficiency = [basal[j] ? e_plant : e_animal for j in 1:S], threshold = threshold)
end
