# ==============================================================================
# Max-Plus Productive (flowshops).
# Note: we use the Julia SimpleWeightedGraphs package for converting the flowshop
# to a grap. We use the Julia TikzPictures package instead of GraphPlot.jl for
# displaying them.
# ==============================================================================

using SimpleWeightedGraphs, TikzPictures

# ==============================================================================
# Save a max-plus matrice to a directed graph.

mpgraph(A::ArrMP) = SimpleWeightedDiGraph(A)
mpgraph(A::SpaMP) = SimpleWeightedDiGraph(A)

# ==============================================================================
# shift

# ==============================================================================
# From the matrix of time of tasks (machines, coins) build a max+ linear system

function flowshop(E::ArrMP{T}) where T
    nmach, npiece = size(E)

    A = mpzeros(T, nmach * npiece, nmach * npiece)
    B = mpzeros(T, nmach * npiece, nmach + npiece)
    C = mpzeros(T, nmach + npiece, nmach * npiece)
    D = mpzeros(T, nmach * npiece, nmach * npiece)
    X0 = mpzeros(T, nmach * npiece, 1)

    l = zeros(Int64, 1, nmach)
    d = zeros(Int64, 1, npiece)
    ij = 0
    for i in 1:npiece
        for j in 1:nmach
            ij = i + (j-1) * npiece
            if iszero(E[j, i])
                !iszero(l[j]) && (l[j] += 1)
                !iszero(d[i]) && (d[i] += 1)
            else
                !iszero(l[j]) && (D[ij, ij - l[j]] = E[j, i - l[j]])
                !iszero(d[i]) && (D[ij, ij - d[i] * npiece] = E[j - d[i],i])
                iszero(d[i])  && (B[ij, i] = mp1);  d[i] = 1
                iszero(l[j])  && (B[ij, j + npiece] = mp1); l[j] = 1
            end
        end
        C[i, ij - (d[i] - 1) * npiece] = E[nmach - d[i] + 1, i]
    end
    for j in 1:nmach
        C[j + npiece, j * npiece - l[j] + 1] = E[j, npiece - l[j] + 1]
    end

    mpsyslin(A, B, C, D, X0)
end

# ==============================================================================
#

#function flowshop_graph(E::arrMP, m::arrMP, p::arrMP)
#    (g, T, N)
#end

# ==============================================================================
#

#function flowshop_simu(E::arrMP, m::arrMP, p::arrMP, u::arrMP)
#    mpsimul()
#end

# ==============================================================================
#

#function critical_graph

# ==============================================================================

# SimpleWeightedDiGraph
function LaTeXG(PT::ArrMP, O::Float64)
    # Define colors and shapes for arrows.
    settings="\\tikzset{
              darrow/.style={->,draw=black,line width=1pt},
              barrow/.style={->,draw=blue,line width=1pt},
              yarrow/.style={->,draw=yellow,line width=1pt}}"

    # Define colors for places
    opts="entree/.style={circle,draw=blue,line width=1pt},
          sortie/.style={circle,draw=green,fill=green,line width=1pt},
          vide/.style={circle,draw=black,line width=1pt},
          pleine/.style={circle,draw=green,line width=1pt}"

    nm, np = size(PT)

    code = "%% Piece inputs -- Piece outputs\n"
    for i in 1:np
        code *= "\\node[entree] (pi$(i)) at ($(i),0) {}; \\node[sortie] (po$(i)) at ($(i),$(nm+1+O)) {};\n"
    end

    code *= "%% Machine inputs -- Machine outputs\n"
    for i in 1:nm
        code *= "\\node[entree] (mi$(i)) at (0,$(nm-i+1)) {}; \\node[sortie] (mo$(i)) at ($(np+1+O),$(nm-i+1)) {};\n"
    end

    code *= "%% Places\n"
    for j in 1:np
        code *= "% Row $(j)\n"
        for i in 1:nm
            if !iszero(PT[i,j])
                code *= "\\node[vide] (p$(i)$(j)) at ($(j+O),$(i+O)) {$(i)$(j)};\n"
            end
        end
        for i in 1:nm
            code *= "\\draw[yarrow] (mo$(i)) -- (mi$(i));\n"
        end
        code *= "\\draw[yarrow] (po$(j)) -- (pi$(j));\n"
    end

    code *= "%% Transitions\n"

    TikzPicture(code, options=opts, preamble=settings)
end
