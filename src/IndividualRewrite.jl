module IndividualRewrite
export Rule, sim
using Catlab.CategoricalAlgebra, Catlab.Present, Catlab.Theories
using Random: shuffle

# Data structures
#################

mutable struct Rule
  L::ACSetTransformation
  R::ACSetTransformation
  delay::Int
  frequency::Int
end

mutable struct RuleCount
  r::Rule
  count::Int
end

RuleCount(r::Rule) = RuleCount(r, 0)

mutable struct MatchCount
  name::Symbol
  rule::Rule
  match::ACSetTransformation
  count::Int
end

MatchCount(name::Symbol, r::Rule, m::ACSetTransformation) =
  MatchCount(name, r, m, r.delay)

# Simulation
############
"""
Run simulation for `n` timesteps.
"""
function sim(init_state::StructACSet, rules::Dict{Symbol,Rule}, n::Int;
             verbose=false)
  rule_counts = Dict(k=>RuleCount(v) for (k,v) in collect(rules))
  matches, states = [], [init_state]
  for step_i in 1:n
    if verbose
      println("\nStep $step_i ($(length(matches)) rewrites queued)")
    end
    # Look for matches if frequency count has reached zero
    for (r_name,r) in collect(rule_counts)
      if r.count > 0
        r.count -= 1
      else
        r.count = r.r.frequency
        rmatches = shuffle(homomorphisms(codom(r.r.L), states[end]))
        if !isempty(rmatches)
          if verbose println("\tAdding match for rewrite $r_name") end
          push!(matches, MatchCount(r_name, r.r, first(rmatches)))
        end
      end
    end
    # Execute matches if their delay count has reached zero
    for (i,m) in enumerate(matches)
      if isnothing(m)
        continue
      elseif m.count > 0
        m.count -= 1
      else
        matches[i] = nothing
        _, kg,_,kh = rewrite_match_maps(m.rule.L, m.rule.R, m.match)
        push!(states, codom(kh))

        if verbose
          println("\tnew state after firing $(m.name)")
          print("\t"); show(stdout, "text/plain", states[end])
        end
        for (j, m_) in enumerate(matches )
          matches[j] = postcompose_partial(kg, kh, m_)
        end
      end
    end
    matches = filter(x->!isnothing(x), matches) # remove executed/invalid match
  end
  return states
end


"""
Convert a match L->G to a match L->H using a partial morphism G->H, if possible.
       L ========= L
     m ↓           ↓ m'
       G <-- K --> H
"""
function postcompose_partial(
    kg::ACSetTransformation,kh::ACSetTransformation, mc::MatchCount)
  println("")
  m = mc.match
  d = Dict()
  for (k,vs) in pairs(components(m))
    vs_ = Int[]
    for v in collect(vs)
      kv = findfirst(==(v), collect(kg[k]))
      if isnothing(kv)
        mc = nothing
        return nothing
      else
        push!(vs_, kh[k](kv))
      end
    end
    d[k] = vs_
  end
  MatchCount(mc.name,mc.rule,ACSetTransformation(dom(m), codom(kh); d...), mc.count)
end

postcompose_partial(::ACSetTransformation,::ACSetTransformation,::Nothing) =
  nothing

end # module