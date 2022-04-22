using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

# Bernoulli sample of vector `m` (all same probability)
function sample_matches(m, r::AbstractFloat, Δt)
    p = cdf(Exponential(), r * Δt)
    randsubseq(m, p)
end

# sample with varying probabilities
function sample_matches(m::AbstractVector{T}, r::AbstractVector{R}, Δt) where {T, R <: AbstractFloat}
    p = cdf(Exponential(), r * Δt)
    runif = rand(length(m))
    samp_idx = Int64[]
    for i in 1:length(m)
        if runif[i] < p[i]
            push!(samp_idx, i)
        end 
    end

    if length(samp_idx) > 0
        return m[samp_idx]
    else 
        return T[]
    end
end


# update matches; return nothing if invalid
function postcompose_partial(kg::ACSetTransformation, kh::ACSetTransformation, m::ACSetTransformation)
  d = Dict()
  for (k,vs) in pairs(components(m))
    vs_ = Int[]
    for v in collect(vs)
      kv = findfirst(==(v), collect(kg[k]))
      if isnothing(kv)
        return nothing
      else
        push!(vs_, kh[k](kv))
      end
    end
    d[k] = vs_
  end
  return ACSetTransformation(dom(m), codom(kh); d...)
end

# stuff to store stuff
struct Rule
    L::ACSetTransformation
    R::ACSetTransformation
end

mutable struct MatchedRule
    rule::Rule
    match::Union{Nothing, ACSetTransformation}
end