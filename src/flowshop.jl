# ==============================================================================
# Max-Plus Productive (flowshops).
# ==============================================================================

using LightGraphs, SimpleWeightedGraphs

# ==============================================================================
#

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

function _save(filename::AbstractString, settings::AbstractString, options::AbstractString, code::AbstractString)
    tex = open("$(filename).tex", "w")
    println(tex, "\\documentclass{article}")
    println(tex, "\\usepackage{caption}")
    println(tex, "\\usepackage{tikz}")
    println(tex, "\\begin{document}")
    if length(settings) != 0
        print(tex, "\\tikzset{")
        println(tex, settings)
        println(tex, "}")
    end
    print(tex, "\\begin{tikzpicture}[")
    print(tex, options)
    println(tex, "]")
    println(tex, code)
    println(tex, "\\end{tikzpicture}")
    println(tex, "\\end{document}")
    close(tex)
end

tikzset="arrow/.style={->,draw=black,line width=1pt},
         barrow/.style={->,draw=blue,line width=1pt},
         yarrow/.style={->,draw=yellow,line width=1pt},"

options="input/.style={circle,draw=blue,line width=1pt},
         place/.style={circle,draw=black,line width=1pt},
         plein/.style={circle,draw=green,fill=green,line width=1pt},
         output/.style={circle,draw=green,line width=1pt},
         transition/.style={rectangle,draw=black!50,fill=black!20,thick}"

function LaTeX(PT::ArrMP, O::Float64)
    (nm,np)=size(PT)

    code = "%% Piece inputs -- Piece outputs\n"
    for i in 1:np
        code *= "\\node[input] (pi$(i)) at ($(i),0) {}; \\node[plein] (po$(i)) at ($(i),$(nm+1+O)) {};\n"
    end

    code *= "%% Machine inputs -- Machine outputs\n"
    for i in 1:nm
        code *= "\\node[input] (mi$(i)) at (0,$(nm-i+1)) {}; \\node[output] (mo$(i)) at ($(np+1+O),$(nm-i+1)) {};\n"
    end

    code *= "%% Places\n"
    for j in 1:np
        code *= "% Row $(j)"
        for i in 1:nm
            if PT[i,j] != -Inf #mpzero(T)
                code *= "\n\\node[place] (p$(i)$(j)) at ($(j+O),$(i+O)) {$(i)$(j)};"
            end
        end
        code *= "\n"
        for i in 1:nm
            code *= "\\draw[yarrow] (mo$(i)) -- (mi$(i));\n"
        end
        code *= "\\draw[yarrow] (po$(j)) -- (pi$(j));\n"
    end

    code *= "%% Transitions\n"

    _save("/tmp/graph", tikzset, options, code)
    #_compile()
end
