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
    (S,I,R,Agent)::Ob
    s::Hom(S,Agent)
    i::Hom(I,Agent)
    r::Hom(R,Agent)
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