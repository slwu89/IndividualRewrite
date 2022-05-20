using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

using Plots
using ProgressBars

include("utils.jl")

# Lotka-Volterra 
@present ThLV(FreeSchema) begin
    (Rabbit,Fox)::Ob
end

@acset_type LV(ThLV)

# birth of prey
L_birth = @acset LV begin Rabbit=1 end
I_birth = LV()
R_birth = @acset LV begin Rabbit=2 end

l_birth = ACSetTransformation(I_birth, L_birth)
r_birth = ACSetTransformation(I_birth, R_birth)

# predation
L_predation = @acset LV begin Rabbit=1; Fox=1 end
I_predation = @acset LV begin Fox=1 end
R_predation = @acset LV begin Fox=2 end

l_predation = ACSetTransformation(I_predation, L_predation; Fox=[1])
r_predation = ACSetTransformation(I_predation, R_predation; Fox=[1])

# death of predator
L_death = @acset LV begin Fox=1 end
I_death = LV()
R_death = LV()

l_death = ACSetTransformation(I_death, L_death)
r_death = ACSetTransformation(I_death, R_death)

# rules
rule_birth = Rule(l_birth, r_birth)
rule_predation = Rule(l_predation, r_predation)
rule_death = Rule(l_death, r_death)


# parameters
X = 50
Y = 10

Δt = 0.1
tmax = 365*2
steps = Int(tmax/Δt)

b = 1.0
mu = 0.01
β = 0.0025

# run model
state = @acset LV begin Rabbit=X; Fox=Y end

out = fill(-1, (steps, 2))

for t = ProgressBar(1:steps)

    # possible events
    birth = homomorphisms(L_birth, state, monic = true)
    predation = homomorphisms(L_predation, state, monic = true)
    death = homomorphisms(L_death, state, monic = true)

    # sample firings
    birth = sample_matches(birth, b, Δt)
    predation = sample_matches(predation, β, Δt)
    death = sample_matches(death, mu, Δt)

    # apply events
    events = MatchedRule[
        [MatchedRule(rule_birth, m) for m in birth]...,
        [MatchedRule(rule_predation, m) for m in predation]...,
        [MatchedRule(rule_death, m) for m in death]...
    ]
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
    out[t, :] = [nparts(state, x) for x in [:Rabbit, :Fox]]
end

plot(
    (1:steps) * Δt,
    out,
    label=["Rabbit" "Fox"],
    xlabel="Time",
    ylabel="Number"
)