# ==============================================================================
# Max-Plus eigenvalues eigenvectors (Howard algorithm).
#
# The Max-Plus version of Howard's algorithm is described in the paper:
#    "Numerical computation of spectral elements in max-plus algebra"
# by Jean Cochet-Terrasson, Guy Cohen, Stephane Gaubert, Michael Mc Gettrick,
# Jean-Pierre Quadrat
#
# 1997: Stephane.Gaubert@inria.fr: Original C version for ScicosLab
#       http://www.cmap.polytechnique.fr/~gaubert/HOWARD2.html
# 2021: quentin.quadrat@gmail.com: Julia portage
# ==============================================================================

"""
    (IJ, A, nnodes, narcs) = spget(S::SparseMatrixCSC)

Sparse description of a matrix.

Inspired by Silab [ij,a,s]=spget(sparse([1 %0; 3 4])) but return (I,J)
on a vector, MaxPlus number for A are converted into classic number and
number of node and number of arcs instead of matrix dimension.

# Examples
```julia-repl
julia> S = sparse([1.0 mp0; 3 4])
2×2 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  1
  [2, 1]  =  3
  [2, 2]  =  4

julia> spget(S)
([1, 1, 2, 1, 2, 2], [1.0, 3.0, 4.0], 2, 3)
```
"""
function spget(S::SparseMatrixCSC)
    if ((size(S,1) != size(S,2)) || (size(S,1) == 0))
        throw(DimensionMismatch("Matrix shall be squared and not empty"))
    end
    (I,J,A) = findnz(S)
    return reshape([I'; J'], :), MaxPlus.plustimes(A), size(S,1), size(I,1)
end

# -----------------------------------------------------------------------------
# Howard results
# -----------------------------------------------------------------------------
struct HowardResult
    # Cycle time vector (dimension: nnodes)
    eigenvalues::Vector{MP}
    # Bias vectors (dimension: nnodes)
    eigenvectors::Vector{MP}
    # Optimal policy (dimension: nnodes)
    policy::Vector{Int64}
    # Number of connected components of the optimal policy which is returned.
    components::Int64
    # Number of iterations the algorithm took (for debug purpose)
    iterations::Int64
end

# -----------------------------------------------------------------------------
#function Base.show(io::IO, ::MIME"text/plain", R::HowardResult)
#    print(io, "Howard algorithm results:\n")
#    print(io, "  - Eigenvalues: ")
#    show(io, R.eigenvalues)
#    print(io, "\n  - Eigenvectors: ")
#    show(io, R.eigenvectors)
#    print(io, "\n  - Optimal policy: ")
#    show(io, R.policy)
#    print(io, "\n  - Number of connected components: ")
#    show(io, R.components)
#    print(io, "\n  - Algorithm iterations taken: ")
#    show(io, R.iterations)
#    print(io, "\n")
#end

# -----------------------------------------------------------------------------
# Howard context
# -----------------------------------------------------------------------------
struct Context
    # INPUTS: Sparse matrix given as input in howard() function
    ij::Vector{Int64}  # Index of the sparse matrix
    A::Vector{Float64} # Values in classic algebra of the sparse matrix
    nnodes::Int        # Number of nodes
    narcs::Int         # Number of arcs

    # INTERNALS
    newpi::Vector{Int64}
    piinv_idx::Vector{Int64}
    piinv_succ::Vector{Int64}
    piinv_elem::Vector{Int64}
    piinv_last::Vector{Int64}
    visited::Vector{Bool}
    component::Vector{Int64}
    c::Vector{Float64}
    newc::Vector{Float64}
    vaux::Vector{Float64}
    newchi::Vector{Float64}

    # OUTPUTS
    v::Vector{Float64}     # Bias vectors
    chi::Vector{Float64}   # Cycle time vector
    pi::Vector{Int64}      # Optimal policy

    function Context(S::SparseMatrixCSC{T,U}) where {T,U}
        (IJ, A, nnodes, narcs) = spget(dropzeros(S))
        new(# INPUTS
            IJ, A, nnodes, narcs,
            # INTERNALS
            zeros(Int64, nnodes),
            zeros(Int64, nnodes),
            zeros(Int64, nnodes),
            zeros(Int64, nnodes),
            zeros(Int64, nnodes),
            zeros(Bool, nnodes),
            zeros(Int64, nnodes),
            zeros(Float64, nnodes),
            zeros(Float64, nnodes),
            fill(-Inf, (nnodes,)),
            zeros(Float64, nnodes),
            # OUTPUTS
            zeros(Float64, nnodes),
            zeros(Float64, nnodes),
            zeros(Int64, nnodes))
    end
end

# -----------------------------------------------------------------------------
function show_inputs(con::Context)
    (@printf "narcs:%d   nnodes:%d\n" con.narcs con.nnodes)
    for i = 1:con.narcs
        (@printf "A[%d]: %d\n"  i  con.A[i])
    end
    for i = 1:2*con.narcs
        (@printf "ij[%d]: %d\n"  i  con.ij[i])
    end
end

# -----------------------------------------------------------------------------
function show_info(con::Context, iteration::Int)
    (@printf "*** ITERATION %d of Max Plus Howard Algorithm *** \n" iteration)
    (@printf "PI=")
    show(stdout, con.pi)
    (@printf "\nc=")
    show(stdout, con.c)
    (@printf "\nchi=")
    show(stdout, con.chi)
    (@printf "\nv=")
    show(stdout, con.v)
    (@printf "\n")
end

# -----------------------------------------------------------------------------
function sanity_checks(con::Context)
    # Check if nodes have successor
    u = zeros(Bool, con.nnodes)
    for i = 1:con.narcs
        u[con.ij[2i-1]] = true
    end

    for i = 1:con.nnodes
        if (u[i] == false)
            throw(InitError("Node " * string(i) * " has no successor"))
        end
    end
end

# -----------------------------------------------------------------------------
# This epsilon value is somehow arbitrary.
function epsilon(con::Context)
    m, M = extrema(con.A)
    (M - m + 1.0) * 0.000000001
end

# -----------------------------------------------------------------------------
function initial_policy(con::Context)
    for i = 1:con.narcs
        if (con.vaux[con.ij[2i-1]] <= con.A[i])
            con.pi[con.ij[2i-1]] = con.ij[2i]
            con.c[con.ij[2i-1]] = con.A[i]
            con.vaux[con.ij[2i-1]] = con.A[i]
        end
    end
end

# -----------------------------------------------------------------------------
function update_policy(con::Context)
    con.pi .= con.newpi
    con.c .= con.newc
    con.vaux .= con.v
end

# -----------------------------------------------------------------------------
function build_inverse(con::Context)
    con.piinv_idx .= -1
    con.piinv_last .= -1
    ptr = 1
    for i = 1:con.nnodes
        j = con.pi[i]
        if (con.piinv_idx[j] == -1)
            con.piinv_succ[ptr] = -1
            con.piinv_elem[ptr] = i
            con.piinv_last[j] = ptr
            con.piinv_idx[j] = ptr
        else
            con.piinv_succ[ptr] = -1
            con.piinv_elem[ptr] = i
            locus = con.piinv_last[j]
            con.piinv_succ[locus] = ptr
            con.piinv_last[j] = ptr
        end
        ptr += 1
    end
end

# -----------------------------------------------------------------------------
function depth_first_label(con::Context, lambda::Float64, color::Int, i::Int)
    a = con.piinv_idx[i]
    while ((a != -1) && (con.visited[con.piinv_elem[a]] == false))
        nexti = con.piinv_elem[a]
        con.visited[nexti] = true
        con.v[nexti] = -lambda + con.c[nexti] + con.v[i]
        con.component[nexti] = color
        con.chi[nexti]= lambda
        depth_first_label(con, lambda, color, nexti)
        a = con.piinv_succ[a]
    end
end

# -----------------------------------------------------------------------------
function visit_from(con::Context, initialpoint::Int, color::Int)
    index = initialpoint
    con.component[index] = color
    newindex = con.pi[index]
    while true
        con.component[newindex] = color
        index = newindex
        newindex = con.pi[index]
        (con.component[newindex] == 0) || break;
    end

    # A cycle has been detected, since newindex is already visited
    weight::Float64 = 0
    length::Float64 = 0
    i = index

    while true
        weight += con.c[i]
        length += 1
        i = con.pi[i]
        (i != index) || break
    end
    lambda = weight / length
    con.v[i] = con.vaux[i]
    con.chi[i] = lambda
    depth_first_label(con, lambda, color, index)
end

# -----------------------------------------------------------------------------
function value(con::Context)
    components = 0
    initialpoint = 1
    color = 1

    # init_depth_first()
    con.visited .= false
    con.component .= 0

    while true
        visit_from(con, initialpoint, color)
        while ((initialpoint <= con.nnodes) && (con.component[initialpoint] != 0))
            initialpoint += 1
        end
        color += 1
        (initialpoint <= con.nnodes) || break
    end

    components = color - 1
end

# -----------------------------------------------------------------------------
function first_order_improvement(con::Context)
    improved = false
    for i = 1:con.narcs
        if (con.chi[con.ij[2i]] > con.newchi[con.ij[2i-1]])
            improved = true;
            con.newpi[con.ij[2i-1]] = con.ij[2i]
            con.newchi[con.ij[2i-1]] = con.chi[con.ij[2i]]
            con.newc[con.ij[2i-1]] = con.A[i];
        end
    end
    return improved
end

# -----------------------------------------------------------------------------
function second_order_improvement(con::Context, components::Int, ϵ::Float64)
    improved = false
    if (components > 1)
        for i = 1:con.narcs
            if (con.chi[con.ij[2i]] == con.newchi[con.ij[2i-1]])
                w = con.A[i] + con.v[con.ij[2i]] - con.chi[con.ij[2i-1]]
                if (w > con.vaux[con.ij[2i-1]] + ϵ)
                    improved = true
                    con.vaux[con.ij[2i-1]] = w
                    con.newpi[con.ij[2 - 1]] = con.ij[2i]
                    con.newc[con.ij[2i-1]] = con.A[i]
                end
            end
        end
    else
        for i = 1:con.narcs
            w = con.A[i]+ con.v[con.ij[2i]] - con.chi[con.ij[2i-1]]
            if (w > con.vaux[con.ij[2i-1]] + ϵ)
                improved = true
                con.vaux[con.ij[2i-1]] = w
                con.newpi[con.ij[2i-1]] = con.ij[2i]
                con.newc[con.ij[2i-1]] = con.A[i]
            end
        end
    end
    return improved
end

# -----------------------------------------------------------------------------
function improve(con::Context, components::Int, ϵ::Float64)
    improved = false

    con.newchi .= con.chi
    con.vaux .= con.v
    con.newpi .= con.pi
    con.newc .= con.c

    if (components > 1)
        improved = first_order_improvement(con)
    end

    if (improved == false)
        improved = second_order_improvement(con, components, ϵ)
    end
    return improved
end

# -----------------------------------------------------------------------------
function howard(S::SparseMatrixCSC{MP}, max_iterations::Int64 = 1000) where MP
    improved::Bool = true
    iterations::Int64 = 0
    components::Int64 = 0

    con = Context(S)
    # show_inputs(con)
    sanity_checks(con)
    ϵ = epsilon(con)
    initial_policy(con)
    build_inverse(con)

    while improved
        components = value(con)
        #show_info(con, iterations)
        improved = improve(con, components, ϵ)
        update_policy(con)
        build_inverse(con)
        iterations += 1

        if (iterations >= max_iterations)
            throw(AssertionError("Max iterations reached"))
        end
    end

    if (components > 0)
        return HowardResult(MP(con.chi), MP(con.v), con.pi, components, iterations)
    end
    return HowardResult(MP([]), MP([]), [], 0, iterations)
end

# -----------------------------------------------------------------------------
function mpeigen(S::SparseMatrixCSC{MP})
    r = howard(S)
    return r.eigenvalues, r.eigenvectors
end

# -----------------------------------------------------------------------------
function mpeigen(S::Matrix{MP})
    r = howard(sparse(S))
    return r.eigenvalues, r.eigenvectors
end
