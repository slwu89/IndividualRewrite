using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

N = 10
I0 = 3
S0 = N - I0
Δt = 0.1
tmax = 100
steps = Int(tmax/Δt)
γ = 1/10 # recovery rate
R0 = 2.5
β = R0 * γ # R0 for corresponding ODEs

initial_states = fill("S", N)
initial_states[rand(1:N, I0)] .= "I"
state_labels = ["S", "I", "R"];

@present ThSIR(FreeSchema) begin
    (S,I,R)::Ob
end

@acset_type SIR(ThSIR)



state = @acset SIR begin S=S0; I=I0; R=N-S0-I0 end

# infection rules
I = @acset SIR begin I=1 end 
I2 = @acset SIR begin I=2 end 
SI = @acset SIR begin S=1; I=1 end
L_infect =  ACSetTransformation(I, SI; I=[1])
R_infect = ACSetTransformation(I, I2; I=[1])

# recovery rules
R = @acset SIR begin R=1 end
Empty = @acset SIR begin end
L_recovery = ACSetTransformation(Empty, I)
R_recovery = ACSetTransformation(Empty, R)


rewrite_match(L_recovery, R_recovery, state)

match = homomorphism(SI, state)
match = homomorphisms(SI, state)

match = homomorphism(R, state)
rewrite(L_recovery, R_recovery, state)