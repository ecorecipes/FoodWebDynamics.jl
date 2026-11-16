# Observed food webs from GATEWAy v1.0 (Brose et al. 2019, Nature Ecology & Evolution;
# doi:10.25829/IDIV.283-3-756, CC BY 4.0): 290 webs, one row per feeding link, with consumer and
# resource body masses. See lab/data/gateway/README.md for provenance and data problems.

const GATEWAY_NONLIVING = ("detritus", "dead organic material", "other")

"""
    ObservedWeb

An observed food web: `name`, `ecosystem`, the feeding network `web` (rows eat columns),
`taxa`, mean body masses `mass` in grams (`NaN` where unrecorded), `living` (false for
detritus and other non-living resources) and `metabolic_type` per node.
"""
struct ObservedWeb
    name::String
    ecosystem::String
    web::FoodWeb
    taxa::Vector{String}
    mass::Vector{Float64}
    living::BitVector
    metabolic_type::Vector{String}
end

"""
    read_gateway(path; detritus = :basal)

The food webs of GATEWAy v1.0 from its link table (`283_2_FoodWebDataBase_2018_12_10.csv`).
Names are stripped of surrounding spaces and exact duplicate rows dropped. A node's mass is the
geometric mean of the masses recorded for it (masses are per link and life stage; −999 means
missing). Non-living resources (metabolic type detritus, dead organic material or other) are
kept as basal nodes (`detritus = :basal`) or removed with their links (`:drop`). Returns a
vector of [`ObservedWeb`](@ref).
"""
function read_gateway(path::AbstractString; detritus::Symbol = :basal)
    detritus in (:basal, :drop) || throw(ArgumentError("detritus must be :basal or :drop, got $detritus"))
    cols = ["foodweb.name", "ecosystem.type", "con.taxonomy", "res.taxonomy", "con.mass.mean.g.", "res.mass.mean.g.",
            "con.metabolic.type", "res.metabolic.type"]
    df = CSV.read(path, DataFrame; select = cols, stringtype = String, missingstring = ["NA", ""])
    clean(x) = x === missing ? "" : strip(String(x))
    mass(x) = x === missing || x == -999 || !(x > 0) ? NaN : Float64(x)
    out = ObservedWeb[]
    for g in groupby(df, "foodweb.name")
        rows = unique(DataFrame((c => clean.(g[!, c]) for c in cols[[1, 2, 3, 4, 7, 8]])...,
                                "con.mass" => mass.(g[!, "con.mass.mean.g."]), "res.mass" => mass.(g[!, "res.mass.mean.g."])))
        logm = Dict{String,Vector{Float64}}()
        mtype = Dict{String,String}()
        for r in eachrow(rows)
            for (who, m, t) in ((r["con.taxonomy"], r["con.mass"], r["con.metabolic.type"]),
                                (r["res.taxonomy"], r["res.mass"], r["res.metabolic.type"]))
                isempty(who) && continue
                v = get!(logm, who, Float64[])
                isnan(m) || push!(v, log(m))
                haskey(mtype, who) && mtype[who] != "" || (mtype[who] = lowercase(t))
            end
        end
        nonliving(t) = mtype[t] in GATEWAY_NONLIVING
        taxa = sort(collect(keys(logm)))
        detritus === :drop && filter!(t -> !nonliving(t), taxa)
        pos = Dict(t => k for (k, t) in enumerate(taxa))
        n = length(taxa)
        A = falses(n, n)
        for r in eachrow(rows)
            c, s = r["con.taxonomy"], r["res.taxonomy"]
            haskey(pos, c) && haskey(pos, s) && (A[pos[c], pos[s]] = true)
        end
        for t in taxa
            nonliving(t) && (A[pos[t], :] .= false)     # non-living nodes eat nothing
        end
        m = [isempty(logm[t]) ? NaN : exp(mean(logm[t])) for t in taxa]
        push!(out, ObservedWeb(g[1, "foodweb.name"] |> clean, rows[1, "ecosystem.type"], FoodWeb(A), taxa, m,
                               BitVector([!nonliving(t) for t in taxa]), [mtype[t] for t in taxa]))
    end
    return out
end

"""
    trophic_species(w::ObservedWeb)

Aggregate taxa with identical prey and identical predators (and the same living status) into
trophic species, the unit with which Williams & Martinez (2000) compared the niche model with
observed webs. Repeat until no two nodes share both sets (merging can make further nodes
identical). A trophic species' name joins its taxa with `" | "`, its mass is the geometric mean
of its taxa's recorded masses, and its metabolic type is that of its first taxon.

In the bioenergetic model, taxa with the same links and trophic level get identical parameters.
If two of them survive, the equilibrium is neutral (the Jacobian has an eigenvalue 0), so
observed webs must be aggregated before computing the stability metrics.
"""
function trophic_species(w::ObservedWeb)
    groups = [[i] for i in eachindex(w.taxa)]
    A = w.web.A
    while true
        sig = [(A[i, :], A[:, i], w.living[g[1]]) for (i, g) in zip(axes(A, 1), groups)]
        keys_ = unique(sig)
        length(keys_) == length(sig) && break
        k = Dict(s => j for (j, s) in enumerate(keys_))
        new = [Int[] for _ in keys_]
        for (i, s) in enumerate(sig)
            append!(new[k[s]], groups[i])
        end
        cls = [k[s] for s in sig]
        B = falses(length(keys_), length(keys_))
        for i in axes(A, 1), j in axes(A, 2)
            A[i, j] && (B[cls[i], cls[j]] = true)
        end
        A, groups = B, new
    end
    gmean(v) = (k = filter(!isnan, v); isempty(k) ? NaN : exp(mean(log.(k))))
    return ObservedWeb(w.name, w.ecosystem, FoodWeb(A), [join(w.taxa[g], " | ") for g in groups],
                       [gmean(w.mass[g]) for g in groups], BitVector([w.living[g[1]] for g in groups]),
                       [w.metabolic_type[g[1]] for g in groups])
end

"""
    relative_masses(w::ObservedWeb; Z = 1.5)

Consumer masses relative to the producers, for `FoodWebModel(...; mass)`: divided by the
geometric mean mass of the basal living nodes with a recorded mass, or, if none has one, by the
smallest consumer mass divided by `Z` (so the lightest consumer sits one niche-model trophic
step above the producers). Returns `nothing` if some consumer has no mass.
"""
function relative_masses(w::ObservedWeb; Z::Real = 1.5)
    basal = isbasal(w.web)
    cons = .!basal
    any(isnan, w.mass[cons]) && return nothing
    known = [w.mass[i] for i in eachindex(w.mass) if basal[i] && w.living[i] && !isnan(w.mass[i])]
    ref = isempty(known) ? minimum(w.mass[cons]) / Z : exp(mean(log.(known)))
    return [basal[i] ? 1.0 : w.mass[i] / ref for i in eachindex(w.mass)]
end
