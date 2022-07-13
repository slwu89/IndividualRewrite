import Graphs as SimpleGraphs
using Plots, GraphRecipes

using Catlab.Graphs.BasicGraphs

using Catlab.Graphics
using Catlab.Graphics.Graphviz

g = SimpleGraphs.erdos_renyi(10, 0.25)

graphplot(g, curves = false)

g1 = Graph(SimpleGraphs.nv(g))

# x, _ = Iterators.peel(SimpleGraphs.edges(g))

elist = [[SimpleGraphs.src(e), SimpleGraphs.dst(e)] for e in SimpleGraphs.edges(g)]
elist = hcat(elist...)

add_edges!(g1, elist[1, :], elist[2, :])
add_edges!(g1, elist[2, :], elist[1, :])

to_graphviz(g1)