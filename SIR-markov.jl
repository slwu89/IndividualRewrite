using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

# Bernoulli sample of vector `m`
function sample_matches(m, r, Δt)
    p = cdf(Exponential(), r * Δt)
    randsubseq(m, p)
end

# apply matches
function apply_matches(state, m, L, R)
    if length(m) > 0
        for i = 1:length(m)
            state = rewrite_match(L, R, m[i])
        end
    end
    return state
end

N = 10
I0 = 3
S0 = N - I0

Δt = 0.1
tmax = 100
steps = Int(tmax/Δt)
γ = 1/10 # recovery rate
R0 = 2.5
β = R0 * γ # R0 for corresponding ODEs

@present ThSIR(FreeSchema) begin
    (S,I,R)::Ob
end

@acset_type SIR(ThSIR)

state = @acset SIR begin S=S0; I=I0; R=N-S0-I0 end

# infection rules
I = @acset SIR begin I=1 end 
I2 = @acset SIR begin I=2 end 
SI = @acset SIR begin S=1; I=1 end
L_infect =  ACSetTransformation(I, SI; I=[1]) # fn from dom I to codom SI, I is the injective function mapping stuff in dom to codom
R_infect = ACSetTransformation(I, I2; I=[1])

# recovery rules
R = @acset SIR begin R=1 end
Empty = @acset SIR begin end
L_recovery = ACSetTransformation(Empty, I)
R_recovery = ACSetTransformation(Empty, R)

match = homomorphism(SI, state)
match = homomorphisms(SI, state)

match = homomorphism(I, state)
match = homomorphisms(I, state)

# rewrite_match(L_recovery, R_recovery, match[1])

# sample_matches(match, γ, Δt)

state

rewrite(L_recovery, R_recovery, state)

# not sure how to update morphisms after the first one was applied.
nparts(state, :R)
state = apply_matches(state,match,L_recovery,R_recovery)
nparts(state, :R)