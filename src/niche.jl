# The niche model of Williams & Martinez (2000), "Simple rules yield complex food webs",
# Nature 404:180–183, doi:10.1038/35004572, and prey-averaged trophic levels.
#
# Each of S species has a niche value nᵢ ~ U(0, 1), a feeding range rᵢ = xᵢ nᵢ with
# xᵢ ~ Beta(1, β), β = 1/(2C) − 1 (so that the expected connectance is C), and a range centre
# cᵢ ~ U(rᵢ/2, nᵢ). Species i eats every species whose niche value lies in
# [cᵢ − rᵢ/2, cᵢ + rᵢ/2]. The species with the smallest niche value gets rᵢ = 0, so every
# web has at least one basal species.

"""
    FoodWeb(A::AbstractMatrix{Bool})

A directed feeding network on `S = size(A, 1)` species, with `A[i, j] == true` when species `i`
eats species `j` (rows are consumers, columns are resources). Basal species (producers) are
those with no resources. Construct one directly or with [`niche_model`](@ref).
"""
struct FoodWeb
    A::BitMatrix
    function FoodWeb(A::AbstractMatrix{Bool})
        S = size(A, 1)
        size(A, 2) == S || throw(ArgumentError("the adjacency matrix must be square, got $(size(A))"))
        return new(BitMatrix(A))
    end
end

"""
    richness(web::FoodWeb)

The number of species in `web`.
"""
richness(web::FoodWeb) = size(web.A, 1)

"""
    connectance(web::FoodWeb)

Directed connectance L/S², with L the number of feeding links (cannibalistic loops included).
"""
connectance(web::FoodWeb) = count(web.A) / richness(web)^2

"""
    isbasal(web::FoodWeb)

A `BitVector` that is `true` for species with no resources (primary producers).
"""
isbasal(web::FoodWeb) = BitVector(vec(.!any(web.A; dims = 2)))

"""
    niche_model([rng], S, C)

Draw a food web of `S` species from the niche model of Williams & Martinez (2000) with expected
directed connectance `C` (`0 < C < 1/2`). The species with the smallest niche value has a
feeding range of zero, so the web has at least one basal species. Species are indexed in
increasing order of niche value. The model does not depend on this order, but the published
median metrics do (see `published_pulse_covariance` in EcologicalResilience.jl). No webs are rejected here;
[`trophic_levels`](@ref) returns `nothing` for webs that a community-assembly step should
reject.
"""
function niche_model(rng::AbstractRNG, S::Integer, C::Real)
    S ≥ 1 || throw(ArgumentError("S must be positive, got $S"))
    0 < C < 0.5 || throw(ArgumentError("C must lie in (0, 1/2), got $C"))
    β = 1 / (2C) - 1
    n = rand(rng, S)
    # Beta(1, β) by inversion: F(x) = 1 − (1 − x)^β
    r = n .* (1 .- (1 .- rand(rng, S)) .^ (1 / β))
    r[argmin(n)] = 0.0
    c = [rᵢ / 2 + rand(rng) * (nᵢ - rᵢ / 2) for (nᵢ, rᵢ) in zip(n, r)]
    o = sortperm(n)
    n, r, c = n[o], r[o], c[o]
    A = falses(S, S)
    for i in 1:S, j in 1:S
        A[i, j] = r[i] > 0 && c[i] - r[i] / 2 ≤ n[j] ≤ c[i] + r[i] / 2
    end
    return FoodWeb(A)
end
niche_model(S::Integer, C::Real) = niche_model(Random.default_rng(), S, C)

"""
    trophic_levels(web::FoodWeb; alive = trues(richness(web)))

Prey-averaged trophic levels among the `alive` species: 1 for a basal species, and one plus the
mean trophic level of its living resources for a consumer. Returns `nothing` when they are
undefined, i.e. when some living consumer has no path to a basal species (for example a closed
loop of consumers). Dead species get `NaN`.
"""
function trophic_levels(web::FoodWeb; alive::AbstractVector{Bool} = trues(richness(web)))
    S = richness(web)
    idx = findall(alive)
    A = web.A[idx, idx]
    k = vec(sum(A; dims = 2))
    n = length(idx)
    # (I − W) L = 1 with W the row-normalised diet matrix; basal rows are identity rows.
    M = Matrix{Float64}(I, n, n)
    for i in 1:n
        k[i] == 0 && continue
        for j in 1:n
            A[i, j] && (M[i, j] -= 1 / k[i])
        end
    end
    F = lu(M; check = false)
    issuccess(F) || return nothing
    L = F \ ones(n)
    all(isfinite, L) && all(≥(1 - 1e-9), L) || return nothing
    out = fill(NaN, S)
    out[idx] = L
    return out
end
