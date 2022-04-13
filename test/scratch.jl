using CombinatorialSpaces.SimplicialSets: get_edge!
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



