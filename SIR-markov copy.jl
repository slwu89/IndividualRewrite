using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

using Plots

include("utils.jl")

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


# struct Rule
#     L::ACSetTransformation
#     R::ACSetTransformation
# end

# mutable struct MatchedRule
#     rule::Rule
#     match::Union{Nothing, ACSetTransformation}
# end


mutable struct rulematch
    m::Union{Nothing, ACSetTransformation}
end

struct queued_updates
    matches::Vector{rulematch}
    L::ACSetTransformation
    R::ACSetTransformation
end

# a very ugly function
function fire_events!(state::ACSet, events::Vector{queued_updates})
    # for each event type
    for i in 1:length(events)
        # for each match
        for j in 1:length(events[i].matches)
            if isnothing(events[i].matches[j])
                continue
            else
                # apply rewrite
                _, kg, _, kh = rewrite_match_maps(events[i].L, events[i].R, events[i].matches[j])
                state = codom(kh)
                # update remaining matches post rewrite
                if j < length(events[i].matches)                    
                    for k in j+1:length(events[i].matches)
                        events[i].matches[j] = postcompose_partial(kg, kh, events[i].matches[j])
                    end
                end
                # update remaining matches for other event types post rewrite
                if i < length(events)
                    for k in i+1:length(events)
                        for l in 1:length(events[k].matches)
                            events[k].matches[l] = postcompose_partial(kg, kh, events[k].matches[l])
                        end
                    end
                end
            end
        end
    end
end

# parameters
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