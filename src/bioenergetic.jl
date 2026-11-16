# The allometric bioenergetic food-web model of Domínguez-García, Dakos & Kéfi (2019), SI §8,
# after Brose, Williams & Martinez (2006) and Kéfi et al. (2016):
#
#     dBᵢ/dt = rᵢ(1 − Bᵢ/K)Bᵢ + Bᵢ Σⱼ eⱼ Fᵢⱼ − Σₖ Bₖ Fₖᵢ − (xᵢ + dᵢ)Bᵢ,
#     Fᵢⱼ    = wᵢ aᵢⱼ Bⱼ^{1+q} / (mᵢ (1 + wᵢ Σₖ aᵢₖ hᵢₖ Bₖ^{1+q})),
#
# with masses mᵢ = Z^{Lᵢ−1} from prey-averaged trophic levels Lᵢ. The details the paper leaves
# open, and how its dataset settles them, are in lab/wiki/model/simulator-reconstruction.md.

"""
    BioenergeticParameters(; kwargs...)

Constants of the bioenergetic model of Domínguez-García et al. (2019, SI §8). Keywords and
defaults:

- `q = 0.3`: Hill coefficient (the functional response has exponent `1 + q`);
- `Z = 1.5`: predator–prey body-mass ratio, `mᵢ = Z^(Lᵢ − 1)`;
- `K = 1.0`: plant carrying capacity; plant growth rate is 1;
- `e_plant = 0.45`, `e_animal = 0.85`: assimilation efficiency on plant and animal prey;
- `x_plant = 0.138`, `x_animal = 0.314`: metabolic constants, `xᵢ = x mᵢ^metabolic_exponent`;
- `metabolic_exponent = -0.25`: printed as `+0.25` in the SI;
- `d0 = 0.1`: natural mortality as a fraction of metabolism, `dᵢ = d0 xᵢ`;
- `a0 = 50.0`, `a0_plant = 10.0`: capture coefficients for animal and plant prey;
- `h0 = 0.4`: handling-time constant (not given in the paper);
- `alpha_consumer = 0.7`, `alpha_resource = 0.4`, `eta_consumer = -0.6`, `eta_resource = 0.5`:
  means of the per-community allometric exponents of capture (α) and handling (η);
- `exponent_cv = 0.1`: their standard deviation as a fraction of the absolute mean;
- `extinction_threshold = 1e-6`: a species is extinct below `extinction_threshold * mᵢ`;
- `plant_metabolism = false`: whether plants pay the metabolic loss `xᵢBᵢ`. SI §8 lists a
  plant metabolic rate, but in Brose et al. (2006) producers have none, and only `false`
  reproduces Dataset S1's total biomass per plant (median 0.98). Natural mortality
  `d0 x mᵢ^metabolic_exponent` applies to plants either way.
- `renormalise_diet = false`: whether a consumer's diet weight `wᵢ` is reset to
  1/(number of living resources) when a resource dies. By default it stays at 1/(number of
  resources in the initial web). SI §9 only deletes the dead species' links. Keeping `wᵢ`
  reproduces Dataset S1's tolerance to extinctions, cascading extinctions and mortality
  tolerance at large richness (median TE 0.332 against 0.334 published at N 61–99;
  renormalising gives 0.393).

The defaults of `metabolic_exponent`, `h0` and `plant_metabolism` were chosen by matching
summaries of Dataset S1 (lab/wiki/model/calibration.md).
"""
Base.@kwdef struct BioenergeticParameters
    q::Float64 = 0.3
    Z::Float64 = 1.5
    K::Float64 = 1.0
    e_plant::Float64 = 0.45
    e_animal::Float64 = 0.85
    x_plant::Float64 = 0.138
    x_animal::Float64 = 0.314
    metabolic_exponent::Float64 = -0.25
    d0::Float64 = 0.1
    a0::Float64 = 50.0
    a0_plant::Float64 = 10.0
    h0::Float64 = 0.4
    alpha_consumer::Float64 = 0.7
    alpha_resource::Float64 = 0.4
    eta_consumer::Float64 = -0.6
    eta_resource::Float64 = 0.5
    exponent_cv::Float64 = 0.1
    extinction_threshold::Float64 = 1e-6
    plant_metabolism::Bool = false
    renormalise_diet::Bool = false
end

"""
    AllometricExponents(alpha_consumer, alpha_resource, eta_consumer, eta_resource)

The exponents of one community: capture `aᵢⱼ ∝ mᵢ^alpha_consumer mⱼ^alpha_resource` and
handling `hᵢⱼ ∝ mᵢ^eta_consumer mⱼ^eta_resource`. The dataset of Domínguez-García et al.
(2019) records them per community as `a_i`, `a_j`, `h_i` and `h_j`.
"""
struct AllometricExponents
    alpha_consumer::Float64
    alpha_resource::Float64
    eta_consumer::Float64
    eta_resource::Float64
end

"""
    draw_exponents([rng], p::BioenergeticParameters)

Draw one community's [`AllometricExponents`](@ref), each normal with mean from `p` and standard
deviation `p.exponent_cv` times the absolute mean.
"""
function draw_exponents(rng::AbstractRNG, p::BioenergeticParameters)
    draw(μ) = μ + p.exponent_cv * abs(μ) * randn(rng)
    return AllometricExponents(draw(p.alpha_consumer), draw(p.alpha_resource), draw(p.eta_consumer),
                               draw(p.eta_resource))
end
draw_exponents(p::BioenergeticParameters) = draw_exponents(Random.default_rng(), p)

"""
    FoodWebModel(web::FoodWeb, p::BioenergeticParameters, ex::AllometricExponents; mass = nothing)

The bioenergetic model on `web`, with rates from `p` and `ex`. Throws an `ArgumentError` if the
trophic levels of `web` are undefined (see [`trophic_levels`](@ref)).

Body masses are `Z^(L − 1)` from prey-averaged trophic levels `L` (SI §9), unless `mass` gives
them, e.g. measured masses for an observed web. Measured masses must be relative to the
producers' mass, because producers have mass 1 in the model; producers' entries of `mass`
are ignored.

The model is mutable state: `alive` marks surviving species, `mortality` holds a multiplier of
each species' natural mortality (all ones at baseline) and `w` the diet weights, `1/(number of
resources)` (living resources only if `params.renormalise_diet`). Change them with [`kill!`](@ref) and [`set_mortality!`](@ref).
Use [`rhs!`](@ref) and [`jacobian!`](@ref) for the vector field and its Jacobian.
"""
mutable struct FoodWebModel
    web::FoodWeb
    params::BioenergeticParameters
    exponents::AllometricExponents
    trophic_level::Vector{Float64}
    mass::Vector{Float64}
    basal::BitVector
    growth::Vector{Float64}
    metabolism::Vector{Float64}
    base_mortality::Vector{Float64}
    mortality::Vector{Float64}
    efficiency::Vector{Float64}
    attack::Matrix{Float64}
    handling::Matrix{Float64}
    resources::Vector{Vector{Int}}
    alive::BitVector
    w::Vector{Float64}
    threshold::Vector{Float64}
    Bp::Vector{Float64}
    dBp::Vector{Float64}
    gbuf::Vector{Float64}
    link_attack::Vector{Vector{Float64}}
    link_ah::Vector{Vector{Float64}}
end

function FoodWebModel(web::FoodWeb, p::BioenergeticParameters, ex::AllometricExponents;
                      mass::Union{Nothing,AbstractVector{<:Real}} = nothing)
    S = richness(web)
    L = trophic_levels(web)
    L === nothing && throw(ArgumentError("trophic levels of this web are undefined"))
    basal = isbasal(web)
    if mass === nothing
        m = [basal[i] ? 1.0 : p.Z^(L[i] - 1) for i in 1:S]
    else
        length(mass) == S || throw(ArgumentError("mass has $(length(mass)) entries for $S species"))
        all(i -> basal[i] || (isfinite(mass[i]) && mass[i] > 0), 1:S) ||
            throw(ArgumentError("consumer masses must be positive and finite"))
        m = [basal[i] ? 1.0 : Float64(mass[i]) for i in 1:S]
    end
    xs = [basal[i] ? p.x_plant : p.x_animal for i in 1:S]
    xfull = xs .* m .^ p.metabolic_exponent
    d = p.d0 .* xfull
    x = [basal[i] && !p.plant_metabolism ? 0.0 : xfull[i] for i in 1:S]
    e = [basal[j] ? p.e_plant : p.e_animal for j in 1:S]
    a = zeros(S, S)
    h = zeros(S, S)
    for i in 1:S, j in 1:S
        web.A[i, j] || continue
        a[i, j] = basal[j] ? p.a0_plant * m[i]^ex.alpha_consumer :
                  p.a0 * m[i]^ex.alpha_consumer * m[j]^ex.alpha_resource
        h[i, j] = p.h0 * m[i]^ex.eta_consumer * m[j]^ex.eta_resource
    end
    res = [findall(web.A[i, :]) for i in 1:S]
    la = [a[i, res[i]] for i in 1:S]
    lah = [a[i, res[i]] .* h[i, res[i]] for i in 1:S]
    model = FoodWebModel(web, p, ex, L, m, basal, Float64.(basal), x, d, ones(S), e, a, h, res,
                         trues(S), zeros(S), p.extinction_threshold .* m, zeros(S), zeros(S), zeros(S), la, lah)
    update_diets!(model)
    return model
end

"""
    richness(model::FoodWebModel)

The number of species still alive in `model`.
"""
richness(model::FoodWebModel) = count(model.alive)

function update_diets!(model::FoodWebModel)
    for i in eachindex(model.w)
        n = model.params.renormalise_diet ? count(j -> model.alive[j], model.resources[i]) :
            length(model.resources[i])
        model.w[i] = n == 0 ? 0.0 : 1 / n
    end
    return model
end

"""
    kill!(model::FoodWebModel, B, i)

Remove species `i`: set `B[i] = 0`, mark it dead and remove its links. Its consumers lose that
share of their foraging, or split it over their remaining resources if
`params.renormalise_diet`. `i` may be a collection.
"""
function kill!(model::FoodWebModel, B::AbstractVector, i)
    for k in i
        B[k] = 0.0
        model.alive[k] = false
    end
    update_diets!(model)
    return model
end

"""
    set_mortality!(model::FoodWebModel, factor)

Set every species' natural mortality to `factor` times its baseline (a scalar applies to all
species, a vector species by species).
"""
set_mortality!(model::FoodWebModel, factor::Real) = (fill!(model.mortality, factor); model)
set_mortality!(model::FoodWebModel, factor::AbstractVector) = (copyto!(model.mortality, factor); model)

"""
    rhs!(du, B, model::FoodWebModel, t)

The vector field of the bioenergetic model in the in-place SciML form. Negative biomasses are
treated as zero in the functional response.
"""
function rhs!(du, B, model::FoodWebModel, t)
    p1 = 1 + model.params.q
    Bp = model.Bp
    @inbounds for j in eachindex(B)
        Bp[j] = model.alive[j] ? max(B[j], 0.0)^p1 : 0.0
    end
    K = model.params.K
    @inbounds for i in eachindex(B)
        loss = model.metabolism[i] + model.mortality[i] * model.base_mortality[i]
        du[i] = model.alive[i] ? (model.growth[i] * (1 - B[i] / K) - loss) * B[i] : 0.0
    end
    @inbounds for c in eachindex(B)
        (model.alive[c] && model.w[c] > 0) || continue
        wc = model.w[c]
        res, la, lah = model.resources[c], model.link_attack[c], model.link_ah[c]
        D = 0.0
        for k in eachindex(res)
            D += lah[k] * Bp[res[k]]
        end
        scale = wc / (model.mass[c] * (1 + wc * D))
        Bc = B[c]
        gain = 0.0
        for k in eachindex(res)
            j = res[k]
            F = scale * la[k] * Bp[j]
            gain += model.efficiency[j] * F
            du[j] -= Bc * F
        end
        du[c] += Bc * gain
    end
    return du
end

"""
    jacobian!(J, B, model::FoodWebModel, t = 0.0)

The Jacobian of [`rhs!`](@ref) at `B`, computed analytically. Rows and columns of dead species
are zero.
"""
function jacobian!(J, B, model::FoodWebModel, t = 0.0)
    q = model.params.q
    p1 = 1 + q
    K = model.params.K
    fill!(J, 0.0)
    Bp, dBp, g = model.Bp, model.dBp, model.gbuf
    @inbounds for j in eachindex(B)
        if model.alive[j]
            b = max(B[j], 0.0)
            bq = b^q
            Bp[j] = bq * b
            dBp[j] = p1 * bq
        else
            Bp[j] = dBp[j] = 0.0
        end
    end
    @inbounds for i in eachindex(B)
        model.alive[i] || continue
        loss = model.metabolism[i] + model.mortality[i] * model.base_mortality[i]
        J[i, i] = model.growth[i] * (1 - 2B[i] / K) - loss
    end
    # For consumer c with D = 1 + w Σₖ a h Bₖ^{1+q} and gⱼ = w a_cj Bⱼ^{1+q}, F_cj = gⱼ/(m D) and
    # ∂F_cj/∂B_l = δⱼₗ w a_cl dB_l/(m D) − gⱼ v_l with v_l = w a_cl h_cl dB_l/(m D²),
    # so each consumer adds a diagonal part and a rank-one block on its resources.
    @inbounds for c in eachindex(B)
        (model.alive[c] && model.w[c] > 0) || continue
        wc = model.w[c]
        res, la, lah = model.resources[c], model.link_attack[c], model.link_ah[c]
        D = 1.0
        for k in eachindex(res)
            D += wc * lah[k] * Bp[res[k]]
        end
        inv = 1 / (model.mass[c] * D)
        Bc = B[c]
        G = 0.0
        for k in eachindex(res)
            j = res[k]
            g[j] = wc * la[k] * Bp[j]
            F = g[j] * inv
            J[c, c] += model.efficiency[j] * F
            J[j, c] -= F
            G += model.efficiency[j] * g[j]
        end
        for k in eachindex(res)
            l = res[k]
            direct = wc * la[k] * dBp[l] * inv
            v = wc * lah[k] * dBp[l] * inv / D
            J[c, l] += Bc * (model.efficiency[l] * direct - G * v)
            J[l, l] -= Bc * direct
            Bcv = Bc * v
            for j in res
                J[j, l] += Bcv * g[j]
            end
        end
    end
    return J
end

"""
    jacobian(B, model::FoodWebModel)

The Jacobian of the living species at `B`, as a dense matrix on `findall(model.alive)`.
"""
function jacobian(B, model::FoodWebModel)
    S = length(B)
    J = zeros(S, S)
    jacobian!(J, B, model)
    idx = findall(model.alive)
    return J[idx, idx]
end
