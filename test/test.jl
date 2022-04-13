using Revise

using IndividualRewrite

using Catlab.Present, Catlab.CategoricalAlgebra
using Plots, Random

Random.seed!(5);

# Define structure of a simulation state
@present ThSIR(FreeSchema) begin
    (S,I,R)::Ob
end
@acset_type SIR(ThSIR)

# Helper function for plotting SIR over time.
function plot_states(states::Vector{SIR})
    ss, is, rs = [[nparts(s, x) for s in states] for x in [:S,:I,:R]]
    x = 1:length(states)
    plot(x,ss; label="S")
    plot!(x, is; label="I")
    plot!(x, rs; label="R")
end

# ACSets for rewrite rules
I = @acset SIR begin I=1 end
R = @acset SIR begin R=1 end
I2 = @acset SIR begin I=2 end
SI = @acset SIR begin S=1; I=1 end

# Rewrite rule morphisms
L_infect =  homomorphism(I, SI)
R_infect = homomorphism(I, I2)
L_recov = homomorphism(SIR(), I)
R_recov = homomorphism(SIR(), R)

# Rewrite rules
rules = Dict(pairs((
  infect=Rule(L_infect,R_infect, 2, 1),
  recover=Rule(L_recov, R_recov, 1, 2)
))...)

# Initial state
init = @acset SIR begin S=10; I=3; R=1 end

# Run the simulation and visualize
states = sim(init, rules, 50; verbose=true)
plot_states(states)