using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq, shuffle!
using StatsBase: sample

using Plots
using ProgressBars

include("utils.jl")

# schema
@present ThSIR(FreeSchema) begin
    (S,I,R,Agent)::Ob
    s::Hom(S,Agent)
    i::Hom(I,Agent)
    r::Hom(R,Agent)
end

@acset_type SIR(ThSIR)

# infection rule
inf_I = @acset SIR begin Agent=2; I=1; i=[1] end
inf_R = @acset SIR begin Agent=2; I=2; i=[1,2] end
inf_L = @acset SIR begin Agent=2; I=1; S=1; i=[1]; s=[2] end
inf_l =  ACSetTransformation(inf_I, inf_L; I=[1], Agent=[1,2])
inf_r = ACSetTransformation(inf_I, inf_R; I=[1], Agent=[1,2])

# recovery rule
rec_L = @acset SIR begin Agent=1; I=1; i=[1] end
rec_R = @acset SIR begin Agent=1; R=1; r=[1] end
rec_I = @acset SIR begin Agent=1 end
rec_l = ACSetTransformation(rec_I, rec_L; Agent=[1])
rec_r = ACSetTransformation(rec_I, rec_R; Agent=[1])

# parameters
N = 1000
I0 = 5
S0 = N - I0

Δt = 0.1
tmax = 100
steps = Int(tmax/Δt)
γ = 1/10 # recovery rate
R0 = 2.5
β = R0 * γ # R0 for corresponding ODEs

# vector of initial states
sir_index = vcat(fill.(["S", "I", "R"],[S0, I0, N-S0-I0])...)
shuffle!(sir_index)

# ACSet to store model state
state = @acset SIR begin Agent=N end
add_parts!(state, :S, sum(sir_index .== "S"))
add_parts!(state, :I, sum(sir_index .== "I"))
add_parts!(state, :R, sum(sir_index .== "R"))

set_subpart!(state, 1:nparts(state, :S), :s, findall(sir_index .== "S"))
set_subpart!(state, 1:nparts(state, :I), :i, findall(sir_index .== "I"))
set_subpart!(state, 1:nparts(state, :R), :r, findall(sir_index .== "R"))

# simulation loop
out = fill(-1, (steps, 3))

for t = ProgressBar(1:steps)
    # matches
    infections_m = homomorphisms(inf_L, state)
    recoveries_m = homomorphisms(rec_L, state)

    # queued updates objects
    infections = queued_updates(infections_m, inf_l, inf_r)
    recoveries = queued_updates(recoveries_m, rec_l, rec_r)

    # sample occurances
    sample_matches(infections, β/N, Δt)
    sample_matches(recoveries, γ, Δt)

    global state = fire_events(state, [infections, recoveries])

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