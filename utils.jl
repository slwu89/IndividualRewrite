using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

# use a bit array to identify which ones will fire and which will not.
struct queued_updates
  matches::Vector{ACSetTransformation{T}} where {T}
  valid::BitVector
  L::ACSetTransformation
  R::ACSetTransformation

  function queued_updates(matches::Vector{ACSetTransformation{T}}, L::ACSetTransformation, R::ACSetTransformation) where {T}
      # no events are scheduled to fire upon creation
      return new(matches, falses(length(matches)), L, R)
  end
end

# Bernoulli sample of vector `m` (all same probability)
# function sample_matches(m::AbstractVector{T}, r::AbstractFloat, Δt) where {T}
#     p = cdf(Exponential(), r * Δt)
#     randsubseq(m, p)
# end

function sample_matches(m::queued_updates, r::AbstractFloat, Δt)
  p = cdf(Exponential(), r * Δt)
  fire = randsubseq(1:length(m.valid), p)
  m.valid[fire] .= true # these will fire
end

# sample with varying probabilities
# function sample_matches(m::AbstractVector{T}, r::AbstractVector{R}, Δt) where {T, R <: AbstractFloat}
#     p = cdf(Exponential(), r * Δt)
#     runif = rand(length(m))
#     samp_idx = Int64[]
#     for i in 1:length(m)
#         if runif[i] < p[i]
#             push!(samp_idx, i)
#         end 
#     end

#     if length(samp_idx) > 0
#         return m[samp_idx]
#     else 
#         return T[]
#     end
# end

function sample_matches(m::queued_updates, r::AbstractVector{R}, Δt) where {T, R <: AbstractFloat}
    p = cdf(Exponential(), r * Δt)
    runif = rand(length(m.valid))
    fire = Int64[]
    for i in 1:length(m.valid)
        if runif[i] < p[i]
            push!(fire, i)
        end 
    end

    m.valid[fire] .= true
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

# # stuff to store stuff
# struct Rule
#     L::ACSetTransformation
#     R::ACSetTransformation
# end

# mutable struct MatchedRule
#     rule::Rule
#     match::Union{Nothing, ACSetTransformation}
# end

# a very ugly function
function fire_events(state::ACSet, events::Vector{queued_updates})
  newstate = state
  # for each event type
  for i in 1:length(events)
    if sum(events[i].valid) == 0
      continue
    end
    # for each match
    for j in 1:length(events[i].matches)
        if events[i].valid[j] === false
            continue
        else
            # apply rewrite
            _, kg, _, kh = rewrite_match_maps(events[i].L, events[i].R, events[i].matches[j])
            newstate = codom(kh)
            # update remaining matches post rewrite
            if j < length(events[i].matches)                    
                for k in j+1:length(events[i].matches)
                # if the event wasn't going to happen anyway, continue
                if events[i].valid[k] === false
                    continue
                else 
                    m = postcompose_partial(kg, kh, events[i].matches[k])
                    if isnothing(m)
                        events[i].valid[k] = false
                    else
                        events[i].matches[k] = m
                    end    
                end                    
                end
            end
            # update remaining matches for other event types post rewrite
            if i < length(events)
                for k in i+1:length(events)
                    for l in 1:length(events[k].matches)
                    # if the event wasn't going to happen anyway, continue
                    if events[k].valid[l] === false
                        continue
                    else
                        m = postcompose_partial(kg, kh, events[k].matches[l])
                        if isnothing(m)
                            events[k].valid[l] = false
                        else
                            events[k].matches[l] = m
                        end       
                    end                                             
                    end
                end
            end
        end
    end
  end
  return newstate
end
