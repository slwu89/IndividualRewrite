using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

using Distributions: Exponential, Geometric, cdf
using Random: randsubseq
using StatsBase: sample

using Plots
using ProgressBars

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

state = @acset SIR begin S=S0; I=I0; R=N-S0-I0 end

# simulation loop
out = fill(-1, (steps, 3))

for t = ProgressBar(1:steps)
    # matches
    infections_m = homomorphisms(SI, state)
    recoveries_m = homomorphisms(I, state)

    # queued updates objects
    infections = queued_updates(infections_m, L_infect, R_infect)
    recoveries = queued_updates(recoveries_m, L_recovery, R_recovery)

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

# debug

# parameters
Δt = 5.0
tmax = 100
steps = Int(tmax/Δt)
γ = 1/10 # recovery rate
R0 = 2.5
β = R0 * γ # R0 for corresponding ODEs

state = @acset SIR begin S=2; I=3; R=1 end

# simulation loop
out = fill(-1, (steps, 3))

# for t = ProgressBar(1:steps)
# matches
infections_m = homomorphisms(SI, state)
recoveries_m = homomorphisms(I, state)

# queued updates objects
infections = queued_updates(infections_m, L_infect, R_infect)
recoveries = queued_updates(recoveries_m, L_recovery, R_recovery)

infections.valid .= true
recoveries.valid .= true

# sample occurances
# sample_matches(infections, β/N, Δt)
# sample_matches(recoveries, γ, Δt)

using Debugger

Debugger.@enter fire_events(state, [infections, recoveries])

#     # write output
#     out[t, :] = [nparts(state, x) for x in [:S, :I, :R]]
# end
