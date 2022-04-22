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

# setup a model with 6 persons
N = 6
sir_index = vcat(fill.(["S", "I", "R"],[3, 2, 1])...)
shuffle!(sir_index)
ages = [50,23,51,49,16,19]

state = AgeSIR{Int64}()
add_parts!(state, :S, sum(sir_index .== "S"))
add_parts!(state, :I, sum(sir_index .== "I"))
add_parts!(state, :R, sum(sir_index .== "R"))

add_parts!(state, :Agent, N)
add_parts!(state, :Age, N)

set_subpart!(state, 1:nparts(state, :S), :s, findall(sir_index .== "S"))
set_subpart!(state, 1:nparts(state, :I), :i, findall(sir_index .== "I"))
set_subpart!(state, 1:nparts(state, :R), :r, findall(sir_index .== "R"))

set_subpart!(state, 1:N, :age, 1:N)
# set_subpart!(state, 1:N, :age, collect(reverse(1:N)))

set_subpart!(state, 1:N, :agevalue, ages)
# set_subpart!(state, 1:N, :agevalue, reverse(ages))

# test infection
# find matches
inf_pairs = homomorphisms(SI, state) 

# check
can_pushout_complement(L_infect, inf_pairs[2])
is_natural(L_infect)
is_natural(R_infect)

# apply the 2nd rewrite
state_new = rewrite_match(L_infect, R_infect, inf_pairs[2])

# age to ID
[incident(state, age, [:age, :agevalue]) for age in ages]
[incident(state_new, age, [:age, :agevalue]) for age in ages]
