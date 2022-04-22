using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq
using Random
using StatsBase: sample

using Distributions: Exponential, Geometric, cdf

using Plots

# age structure: fix age ------------------------------------------------------------
@present ThAgeSIR(FreeSchema) begin
    (S,I,R,Agent,Age)::Ob
     s::Hom(S, Agent)
     i::Hom(I, Agent)
     r::Hom(R, Agent)
     age::Hom(Age, Agent)

    AgeValue::AttrType
    agevalue::Attr(Age, AgeValue)
end

@acset_type AgeSIR(ThAgeSIR)

# setup a model...
sir_index = vcat(fill.(["S", "I", "R"],[3, 2, 1])...)
shuffle!(sir_index)

state = AgeSIR{Int64}()
add_parts!(state, :S, sum(sir_index .== "S"))
add_parts!(state, :I, sum(sir_index .== "I"))
add_parts!(state, :R, sum(sir_index .== "R"))

N = 6
add_parts!(state, :Agent, N)
add_parts!(state, :Age, N)

set_subpart!(state, 1:nparts(state, :S), :s, findall(sir_index .== "S"))
set_subpart!(state, 1:nparts(state, :I), :i, findall(sir_index .== "I"))
set_subpart!(state, 1:nparts(state, :R), :r, findall(sir_index .== "R"))

set_subpart!(state, 1:N, :age, 1:N)

ages = sample(1:80, N)
set_subpart!(state, 1:N, :agevalue, ages)

# infection rules
I = @acset AgeSIR{Int64} begin Agent=2; I=1; i=[2] end # need 2 agents otherwise have dangling edges
I2 = @acset AgeSIR{Int64} begin Agent=2; I=2; i=[1,2] end
SI = @acset AgeSIR{Int64} begin Agent=2; I=1; S=1; s=[1]; i=[2] end

L_infect = ACSetTransformation(I, SI; Agent = [1,2], I = [1])
R_infect = ACSetTransformation(I, I2; Agent = [1,2], I = [2])

# recovery rules
I1 = @acset AgeSIR{Int64} begin Agent=1; I=1; i=[1] end
R1 = @acset AgeSIR{Int64} begin Agent=1; R=1; r=[1] end
A1 = @acset AgeSIR{Int64} begin Agent=1 end

L_recover = ACSetTransformation(A1, I1; Agent = [1])
R_recover = ACSetTransformation(A1, R1; Agent = [1])

# test recovery
is_natural(L_recover)
is_natural(R_recover)

rec_ind = homomorphisms(I1, state)
can_pushout_complement(L_recover, rec_ind[1])

state_new = rewrite_match(L_recover, R_recover, rec_ind[2])

# test infection
# find matches
inf_pairs = homomorphisms(SI, state)

# get the ages of the S and I persons involved in each matching

# ages of S
vcat([state[collect(x[:S]), [:s, :age, :agevalue]] for x in inf_pairs]...)

# ages of I
vcat([state[collect(x[:I]), [:s, :age, :agevalue]] for x in inf_pairs]...)

# check
can_pushout_complement(L_infect, inf_pairs[1]) # fails b/c false
is_natural(L_infect)
is_natural(R_infect)

# apply the 2nd rewrite
state_new = rewrite_match(L_infect, R_infect, inf_pairs[2])

