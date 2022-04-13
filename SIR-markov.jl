using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

using Plots

# Bernoulli sample of vector `m`
function sample_matches(m, r, Δt)
    p = cdf(Exponential(), r * Δt)
    randsubseq(m, p)
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

# schema
@present ThSIR(FreeSchema) begin
    (S,I,R)::Ob
end

@acset_type SIR(ThSIR)

# infection rules
I = @acset SIR begin I=1 end 
I2 = @acset SIR begin I=2 end 
SI = @acset SIR begin S=1; I=1 end
L_infect =  ACSetTransformation(I, SI; I=[1]) # fn from dom I to codom SI, I is the injective function mapping stuff in dom to codom
R_infect = ACSetTransformation(I, I2; I=[1])

# recovery rules
R = @acset SIR begin R=1 end
L_recovery = ACSetTransformation(SIR(), I)
R_recovery = ACSetTransformation(SIR(), R)

# structs
rule_inf = Rule(L_infect, R_infect)
rule_rec = Rule(L_recovery, R_recovery)

N = 500
I0 = 10
S0 = N - I0

Δt = 0.1
tmax = 100
steps = Int(tmax/Δt)
γ = 1/10 # recovery rate
R0 = 2.5
β = R0 * γ # R0 for corresponding ODEs

state = @acset SIR begin S=S0; I=I0; R=N-S0-I0 end

# simulation loop
out = fill(-1, (steps, 3))

for t = 1:steps
    # possible events
    infections = homomorphisms(SI, state)
    recoveries = homomorphisms(I, state)
    # sample occurances
    infections = sample_matches(infections, β/N, Δt)
    recoveries = sample_matches(recoveries, γ, Δt)
    # apply events
    events = MatchedRule[[MatchedRule(rule_inf, m) for m in infections]..., [MatchedRule(rule_rec, m) for m in recoveries]...]
    while length(events) > 0
        ev = pop!(events)
        # apply event
        _, kg, _, kh = rewrite_match_maps(ev.rule.L, ev.rule.R, ev.match)
        state = codom(kh)
        # update the remaining matches post rewrite
        for j in 1:length(events)
            events[j].match = postcompose_partial(kg, kh, events[j].match)
        end
        events = filter((e) -> !isnothing(e.match), events)
    end
    # write output
    out[t, :] = [nparts(state, x) for x in [:S, :I, :R]]
end

plot(
    (1:steps) * Δt,
    out,
    label=["S" "I" "R"],
    xlabel="Time",
    ylabel="Number"
)



# # --------------------------------------------------------
# match = homomorphism(SI, state)
# match = homomorphisms(SI, state)

# match = homomorphism(I, state)
# match = homomorphisms(I, state)

# # rewrite_match(L_recovery, R_recovery, match[1])

# # sample_matches(match, γ, Δt)

# state

# rewrite(L_recovery, R_recovery, state)

# # not sure how to update morphisms after the first one was applied.
# nparts(state, :R)
# state = apply_matches(state,match,L_recovery,R_recovery)
# nparts(state, :R)