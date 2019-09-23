# ==============================================================================
# Max-Plus Productive (flowshops).
# ==============================================================================

using LightGraphs

A = MP([1.0 -Inf; 0.0 2.0])
G1 = mpgraph(A)
G1.weights
@test nonzeros(G1.weights) == MP([1.0; 0.0; 2.0])
@test nv(G1) == 2
@test ne(G1) == 3

G2 = mpgraph(mpsparse(A))
@test G2.weights == G1.weights
@test nv(G2) == nv(G1)
@test ne(G2) == ne(G1)
