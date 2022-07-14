using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

# schema
@present ThSIR(FreeSchema) begin
    (S,I,R,Agent)::Ob
    s::Hom(S,Agent)
    i::Hom(I,Agent)
    r::Hom(R,Agent)
end

@acset_type SIR(ThSIR)

# infection rules
inf_I = @acset SIR begin I=1 end
inf_R = @acset SIR begin I=2 end
inf_L = @acset SIR begin S=1; I=1 end
inf_l =  ACSetTransformation(inf_I, inf_L; I=[1]) # fn from dom I to codom SI, I is the injective function mapping stuff in dom to codom
inf_r = ACSetTransformation(inf_I, inf_R; I=[1])

# recovery rules
rec_L = @acset SIR begin Agent=1; I=1; i=[1] end
rec_R = @acset SIR begin Agent=1; R=1; r=[1] end
rec_l = ACSetTransformation(SIR(), rec_L)
rec_r = ACSetTransformation(SIR(), rec_R)

sir = @acset SIR begin Agent=3; I=2; S=1; i=[1,2]; s=[3] end

m = homomorphisms(rec_L, sir, monic = true)

can_pushout_complement(rec_l, m[1])

ik, kg = pushout_complement(rec_l, m[1])

# shouldn't work
rec_I1 = @acset SIR begin Agent=1 end
rec_l = ACSetTransformation(rec_I1, rec_L; Agent=[1])

m = homomorphisms(rec_L, sir, monic = true)

can_pushout_complement(rec_l, m[1])

ik, kg = pushout_complement(rec_l, m[1])
