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
inf_I = @acset SIR begin Agent=2; I=1; i=[1] end
inf_R = @acset SIR begin Agent=2; I=2; i=[1,2] end
inf_L = @acset SIR begin Agent=2; I=1; S=1; i=[1]; s=[2] end
inf_l =  ACSetTransformation(inf_I, inf_L; I=[1], Agent=[1,2])
inf_r = ACSetTransformation(inf_I, inf_R; I=[1], Agent=[1,2])

sir = @acset SIR begin Agent=4; I=1; S=2; R=1; i=[1]; s=[2,3]; r=[4] end

m = homomorphisms(inf_L, sir, monic = true)

ik, kg = pushout_complement(inf_l, m[1])
K = codom(ik)

sir_inf = rewrite_match(inf_l, inf_r, m[1])



# recovery rule
rec_L = @acset SIR begin Agent=1; I=1; i=[1] end
rec_R = @acset SIR begin Agent=1; R=1; r=[1] end
rec_I = @acset SIR begin Agent=1 end
rec_l = ACSetTransformation(rec_I, rec_L; Agent=[1])
rec_r = ACSetTransformation(rec_I, rec_R; Agent=[1])


# recovery rules 0
rec_L = @acset SIR begin Agent=1; I=1; i=[1] end
rec_R = @acset SIR begin Agent=1; R=1; r=[1] end
rec_l = ACSetTransformation(SIR(), rec_L)
rec_r = ACSetTransformation(SIR(), rec_R)

sir = @acset SIR begin Agent=3; I=2; S=1; i=[1,2]; s=[3] end

m = homomorphisms(rec_L, sir, monic = true)

can_pushout_complement(rec_l, m[1])

ik, kg = pushout_complement(rec_l, m[1])
K = codom(ik)

sir_rec = rewrite_match(rec_l, rec_r, m[1])

# recovery rules 1
rec_I1 = @acset SIR begin Agent=1 end
rec_l1 = ACSetTransformation(rec_I1, rec_L; Agent=[1])
rec_r1 = ACSetTransformation(rec_I1, rec_R; Agent=[1])

sir1 = @acset SIR begin Agent=3; I=2; S=1; i=[1,2]; s=[3] end

m = homomorphisms(rec_L, sir1, monic = true)

can_pushout_complement(rec_l1, m[1])

ik, kg = pushout_complement(rec_l1, m[1])
K1 = codom(ik)

sir1_rec = rewrite_match(rec_l1, rec_r1, m[1])