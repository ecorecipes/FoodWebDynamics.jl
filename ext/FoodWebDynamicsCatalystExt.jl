module FoodWebDynamicsCatalystExt

using FoodWebDynamics
using FoodWebDynamics: FoodWebModel
using Catalyst
import FoodWebDynamics: reaction_network

# The bioenergetic model of a community as a chemical reaction network: one species per living
# population, and reactions for plant growth, metabolic and natural losses, and feeding.
# Feeding of consumer c on resource j is the reaction Bⱼ → eⱼ B_c, so the assimilation
# efficiency is a (non-integer) stoichiometry and the functional response its rate.
function reaction_network(model::FoodWebModel, B = nothing; name::Symbol = :foodweb)
    t = Catalyst.default_t()
    idx = findall(model.alive)
    pos = Dict(i => k for (k, i) in enumerate(idx))
    X = [only(@species $(Symbol(:B, i))(t)) for i in idx]
    μ = [only(@parameters $(Symbol(:μ, i))) for i in idx]
    @parameters q K
    sp(i) = X[pos[i]]
    rxs = Reaction[]
    for i in idx
        x = sp(i)
        model.basal[i] && push!(rxs, Reaction(model.growth[i] * (1 - x / K) * x, nothing, [x]; only_use_rate = true))
        loss = (model.metabolism[i] + μ[pos[i]] * model.base_mortality[i]) * x
        push!(rxs, Reaction(loss, [x], nothing; only_use_rate = true))
    end
    for c in idx
        res = [j for j in model.resources[c] if model.alive[j]]
        isempty(res) && continue
        w = 1 / length(res)
        D = 1 + sum(w * model.attack[c, j] * model.handling[c, j] * sp(j)^(1 + q) for j in res)
        for j in res
            F = w * model.attack[c, j] * sp(j)^(1 + q) / (model.mass[c] * D)
            push!(rxs, Reaction(sp(c) * F, [sp(j)], [sp(c)], [1], [model.efficiency[j]]; only_use_rate = true))
        end
    end
    rn = complete(ReactionSystem(rxs, t, X, [μ; q; K]; name = name, combinatoric_ratelaws = false))
    pmap = Dict([μ[k] => model.mortality[i] for (k, i) in enumerate(idx)]..., q => model.params.q,
                K => model.params.K)
    u0 = B === nothing ? nothing : Dict(X[k] => B[i] for (k, i) in enumerate(idx))
    return (network = rn, species = X, mortality = μ, u0 = u0, parameters = pmap)
end

end
