using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra
using Catlab.Graphs, Catlab.WiringDiagrams
using Catlab.Graphics

G = Graph(5)
add_edges!(G,  [1,1,1,3,4,2],[4,3,2,2,5,5])

K = Graph(2)

L = Graph(3)
add_edges!(L, [1,1,3], [3,2,2])

R = Graph(4)
add_edges!(R, [1,1,4,3], [3,4,2,2])

m = homomorphisms(L, G)

l = CSetTransformation(K, L, V = [1,2])
r = CSetTransformation(K, R, V = [1,2])

H = rewrite_match(l, r, m[1])

to_graphviz(H)