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

# Bernoulli sample of matches in `m` (all same probability)
function sample_matches(m::queued_updates, r::AbstractFloat, Δt)
  p = cdf(Exponential(), r * Δt)
  fire = randsubseq(1:length(m.valid), p)
  m.valid[fire] .= true # these will fire
end

# Bernoulli sample of matches in `m` (unique probabilities)
function sample_matches(m::queued_updates, r::AbstractVector{R}, Δt) where {R <: AbstractFloat}
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

# update matches; return nothing if invalid (is there a way to do this in a type-stable way?)

"""
Convert a match L->G to a match L->H using a partial morphism G->H, if possible.
       L ========= L
     m ↓           ↓ m'
       G <-- K --> H
"""
function postcompose_partial(kg::ACSetTransformation, kh::ACSetTransformation, m::ACSetTransformation)
  m′ = Dict()
  # morphism (m) is all the functions between Obs in L -> G
  for (m_ob, m_fn) in pairs(components(m)) # m_ob: the ob, m_fn: the function for ob(L)->ob(G)
    m′_fn = Int[] # the new function in m′ between L -> H
    for x in collect(m_fn) # the actual function values of m_fn
      kx = findfirst(==(x), collect(kg[m_ob])) # find elem in K getting sent to the same thing as x
      if isnothing(kx) # it could be that thing was deleted, in which case there is no match
        return nothing
      else
        push!(m′_fn, kh[m_ob](kx)) # x now gets sent here in H 
      end
    end
    m′[m_ob] = m′_fn # put the newly constructed fn for that ob into m′
  end
  return ACSetTransformation(dom(m), codom(kh); m′...)
end

# fire all events which are scheduled to fire this time step
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
