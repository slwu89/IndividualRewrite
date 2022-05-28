using Catlab.CategoricalAlgebra, Catlab.Graphs, Catlab.Present, Catlab.Graphics, Catlab.Theories
using Catlab.CategoricalAlgebra.FinCats: FinCatGraphEq
using Random
using StatsBase: sample

using Distributions: Exponential, Geometric, cdf

using Plots
using LinearAlgebra
using ProgressBars

include("utils.jl")

# so when someone goes from S->I, we want to schedule their I->R rewrite. But how do we store that?
# right now, we just update them to R in X time steps, regardless of their current state at the time.