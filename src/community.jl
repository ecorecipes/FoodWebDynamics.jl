# Perturbing a community and returning it to steady state; the ResilienceCore community
# interface for `Community`.

function Base.copy(m::FoodWebModel)
    return FoodWebModel(m.web, m.params, m.exponents, m.trophic_level, m.mass, m.basal, m.growth, m.metabolism,
                        m.base_mortality, copy(m.mortality), m.efficiency, m.attack, m.handling, m.resources,
                        copy(m.alive), copy(m.w), m.threshold, zeros(length(m.Bp)), zeros(length(m.Bp)), zeros(length(m.Bp)),
                        m.link_attack, m.link_ah)
end

"""
    PERTURBATION_STEADY

The steady-state settings used after a perturbation. Integration is limited to `tmax = 10⁴`
(assembly uses T_SAFE = 10⁶). It stops once no species has died for 1000 time units: such a
community is oscillating, and its state at that time is used. SI §10 does not say how
non-stationary perturbed states were treated.
"""
const PERTURBATION_STEADY = SteadyStateOptions(tmax = 1e4, maxiters = 2 * 10^5, stall = 1e3)

# Move `model` (already modified) to its new steady state from `B`. A Newton step from the old
# equilibrium is accepted when it lands on a stable point with all survivors above threshold
# (`newton = true`); otherwise the model is integrated, which also handles extinctions.
function equilibrate!(model::FoodWebModel, B::Vector{Float64}, steady::SteadyStateOptions; newton::Bool = true)
    if newton
        Bn = copy(B)
        alive0 = copy(model.alive)
        if polish!(model, Bn)
            copyto!(B, Bn)
            return true
        end
        copyto!(model.alive, alive0)
    end
    Bs, ok = settle!(model, B; options = steady)
    copyto!(B, Bs)
    return ok
end

state(c::Community) = c.B
present_species(c::Community) = findall(c.model.alive)
jacobian(c::Community) = jacobian(c.B, c.model)

function with_mortality(c::Community, factors::AbstractVector; newton::Bool = true, steady = nothing)
    model = copy(c.model)
    B = copy(c.B)
    set_mortality!(model, factors)
    equilibrate!(model, B, something(steady, PERTURBATION_STEADY); newton = newton)
    return Community(model, B, c.initial_richness, c.dynamic_richness)
end

function without_species(c::Community, i; steady = nothing)
    model = copy(c.model)
    B = copy(c.B)
    kill!(model, B, i)
    equilibrate!(model, B, something(steady, PERTURBATION_STEADY); newton = false)
    return Community(model, B, c.initial_richness, c.dynamic_richness)
end
