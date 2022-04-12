using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories, Catlab.Present
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq


@present ThSIR(FreeSchema) begin
    (S,I,R)::Ob
end

"""This code will work for any choice of 'type graph'"""
function update_state(state::StructACSet{S}) where S  
    rules = []
    while !isempty(rules)
        state = apply_rule(state, pop!(rules))  
    end 
end 

@acset_type SIR(ThSIR)

# infect = Rule("name", L, R)

I = @acset SIR begin I=1 end 
I2 = @acset SIR begin I=2 end 
SI = @acset SIR begin S=1; I=1 end
# L_infect =  homomorphism(I, SI) # equivalent to below b/c only one morphism
L_infect =  ACSetTransformation(I, SI; I=[1])
R_infect = ACSetTransformation(I, I2; I=[1])

ex_state = @acset SIR begin S=10; I=3; R=1 end

# rewrite(L_infect, R_infect, ex_state) # execute immediately

match = homomorphism(SI, ex_state)


matches = [match_1, match_2] # pointing to ex_state


"""
L    I    R 
|
G <- K -> H (updated G)
   kg  kh 

"""


_, kg, _, kh = rewrite_match_maps(L_infect, R_infect, match)

# pretend state_2 is a homomorphism from state_1 -> state_2

matches = [match for match in matches]

state_3 = rewrite_match(L_infect, R_infect, match_2)

# type 2 --------------------------------------------------
@present ThSIR2(FreeSchema) begin
    State::Ob
    Agents::Ob
    state::Hom(Agents, State)
end

@acset_type SIR2(ThSIR2)

state0 = vcat(fill.(1:3, [5,4,1])...)

SIR2_state = @acset SIR2 begin State=3; Agents=10; state=state0 end

I = @acset SIR2 begin Agents=1; State=3; state=2 end 
I2 = @acset SIR2 begin Agents=2; State=3; state=[2,2] end 
SI = @acset SIR2 begin Agents=2; State=3; state=[1,2] end


# L_infect =  homomorphism(I, SI) # equivalent to below b/c only one morphism
homomorphism(I, SI)
homomorphism(I, I2)

L_infect =  ACSetTransformation(I, SI; Agents=[1])
R_infect = ACSetTransformation(I, I2)