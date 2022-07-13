using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq

# schema ------------------------------------------------------------
@present ThVirSIR(FreeSchema) begin
    (S,I,R,Agent)::Ob
    s::Hom(S, Agent)
    i::Hom(I, Agent)
    r::Hom(R, Agent)

    Virulence::AttrType
    virulence::Attr(I, Virulence)
end

to_graphviz(ThVirSIR)

@acset_type VirSIR(ThVirSIR, index = [:s, :i, :r, :virulence])

state1 = @acset VirSIR{Float64} begin Agent=4; S=2; I=1; R=1; s=[1,4]; i=[3]; r=[2]; virulence=[3.5] end
incident(state1, 4, :s) # gives 2
state1[2, :s] # gives 4

# inconsistent state
state2 = @acset VirSIR{Float64} begin Agent=4; S=2; I=1; R=1; s=[1,2]; i=[3]; r=[2]; virulence=[3.5] end
# person 2 is both susceptible and recovered
incident(state2, 2, :s)
incident(state2, 2, :r)

# tests for SIR-age
@present ThAgeSIR(FreeSchema) begin
    (S,I,R,Agent,Age)::Ob
    s::Hom(S, Agent)
    i::Hom(I, Agent)
    r::Hom(R, Agent)
    age::Hom(Age, Agent)

    AgeValue::AttrType
    agevalue::Attr(Age, AgeValue)
end

@acset_type AgeSIR(ThAgeSIR, index = [:s, :i, :r, :age])

# i = 1
# j = 2
# L = @acset AgeSIR{Int64} begin Agent=2; I=1; S=1; Age=2; agevalue=[i,j]; age=[1,2]; s=[2]; i=[1] end
# I = @acset AgeSIR{Int64} begin Agent=2; I=1; Age=2; agevalue=[i,j]; age=[1,2]; i=[1] end
# R = @acset AgeSIR{Int64} begin Agent=2; I=2; Age=2; agevalue=[i,j]; age=[1,2]; i=[1,2] end
# l = ACSetTransformation(I, L; Agent = [1,2], I = [1], Age = [1,2])
# r = ACSetTransformation(I, R; Agent = [1,2], I = [1], Age = [1,2])

# G = @acset AgeSIR{Int64} begin Agent=4; S=2; I=1; R=1; Age=4; agevalue=[i,j,3,4]; age=[1,2,3,4]; s=[2,3]; i=[1]; r=[4] end

# m = homomorphisms(L, G)[1]

# i, g = pushout_complement(l, m)

# C = codom(i)

# rh, ch = pushout(r, i)

# H = codom(rh)

# infection rules
L = @acset AgeSIR{Int64} begin Agent=2; I=1; S=1; s=[1]; i=[2] end
I = @acset AgeSIR{Int64} begin Agent=2; I=1; i=[2] end # need 2 agents otherwise have dangling edges
R = @acset AgeSIR{Int64} begin Agent=2; I=2; i=[1,2] end

l = ACSetTransformation(I, L; Agent = [1,2], I = [1])
r = ACSetTransformation(I, R; Agent = [1,2], I = [2])

G = @acset AgeSIR{Int64} begin Agent=4; S=2; I=1; R=1; Age=4; agevalue=[1,2,3,4]; age=[1,2,3,4]; s=[2,3]; i=[1]; r=[4] end

m = homomorphisms(L, G)[1]

i, g = pushout_complement(l, m)

ik, kg, rh, kh = rewrite_match_maps(l, r, m)

# deletion 1
# G = @acset AgeSIR{Int64} begin Agent=2 end # works because no unknown context

L = @acset AgeSIR{Int64} begin Agent=1 end
I = AgeSIR{Int64}()
R = AgeSIR{Int64}()
l = ACSetTransformation(I, L)
r = ACSetTransformation(I, R)

m = homomorphisms(L, G; monic = true)

ik, kg, rh, kh = rewrite_match_maps(l, r, m[1])

# age and virulence
@present ThComplexSIR(FreeSchema) begin
    (S,I,R,Agent,Age)::Ob
    s::Hom(S, Agent)
    i::Hom(I, Agent)
    r::Hom(R, Agent)
    age::Hom(Age, Agent)

    AgeValue::AttrType
    agevalue::Attr(Age, AgeValue)

    Virulence::AttrType
    virulence::Attr(I, Virulence)
end

to_graphviz(ThComplexSIR)

@acset_type ComplexSIR(ThComplexSIR, index = [:s, :i, :r, :age])

state = @acset ComplexSIR{Int64, Float64} begin Agent=4; Age=4; S=2; I=1; R=1; s=[1,2]; i=[3]; r=[4]; age=[1,2,3,4]; agevalue=[15,39,42,57]; virulence=[0.32,0.13,0.95,0.37] end