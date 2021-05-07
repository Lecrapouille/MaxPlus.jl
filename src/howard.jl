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
# 2021: quentin.quadrat@gmail.com: Julia portage of the C version.

"""
    sparse_ij(S::SpaMP)

Concat indices of a sparse matrix line by line.

# Examples
```julia-repl
julia> A = [1 2; 3 4; 5 6]
3×2 Matrix{Int64}:
 1  2
 3  4
 5  6

julia> vcat(A'...)
6-element Vector{Int64}:
 1
 2
 3
 4
 5
 6

julia> sparse_ij(mpsparse([1.0 2; 3 4]))
```
"""
function sparse_ij(S::SparseMatrixCSC)
    (I,J,) = findnz(S)
    return vcat([J I]'...)
end

# -----------------------------------------------------------------------------
# Howard context
# -----------------------------------------------------------------------------
struct Context
    # INPUTS: Sparse matrix given as input in howard() function
    ij::Vector{Int64}  # Index of the sparse matrix to egein
    A::Vector{Float64} # Values in classic algebra of the sparse matrix
    nnodes::Int      # sparse matrix dimension
    narcs::Int       # sparse matrix dimension

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
    pi::Vector{Int64}  # Optimal policy

    # FIXME transposed matrix ???
    function Context(S::SparseMatrixCSC{T,U}) where {T,U}
        A = convert(SparseMatrixCSC{Float64,U},transpose(S)).nzval
        nnodes = size(S,1)
        narcs = size(A,1)
        new(## INPUTS: Scilab: [ij,a,s]=spget(sparse(M))
            sparse_ij(S),
            A,
            nnodes,
            narcs,
            # INTERNALS
            zeros(Int, nnodes),
            zeros(Int, nnodes),
            zeros(Int, nnodes),
            zeros(Int, nnodes),
            zeros(Int, nnodes),
            zeros(Bool, nnodes),
            zeros(Int, nnodes),
            zeros(Float64, nnodes),
            zeros(Float64, nnodes),
            fill(-Inf, (nnodes,)), # epsilon value
            zeros(Float64, nnodes),
            # OUTPUTS
            zeros(Float64, nnodes),
            zeros(Float64, nnodes),
            fill(673720360, nnodes))
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
function show_info_improve_bias(con::Context, i::Int)
    (@printf "type 2 improvement\n")
    I = con.ij[2i - 1]
    J = con.ij[2i]
    (@printf "Improvement of the BIAS at node %d\n" I)
    (@printf "A[%d] + v[%d] - chi[%d] - v[%d]= %f + %f - %f - %f = %f > 0\n" i J I I con.A[i] con.v[J] con.chi[I] con.v[I] con.A[i]+con.v[J]-con.chi[I]-con.v[I])
end

# -----------------------------------------------------------------------------
function sanity_checks(con::Context)
    # Check matrix dimension
    if ((con.nnodes < 1) || (con.narcs < 1))
        throw(DimensionMismatch("Sparse matrix shall not be empty"))
    end

    # Check if nodes have successor
    u = zeros(Bool, con.nnodes)
    #(@printf "sanity_checks %d:\n" con.narcs)
    for i = 1:con.narcs
        #(@printf "  %d\n" i)
        #(@printf " => ij[%d]: %d\n"  2i-1  con.ij[2i - 1])
        u[con.ij[2i - 1]] = true
    end

    for i = 1:con.nnodes
        if (u[i] == false)
            throw(InitError("Node X has no successor"))
        end
    end
    #(@printf "U=")
    #show(stdout, u)
    #(@printf "\n")
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
        if (con.vaux[con.ij[2i - 1]] <= con.A[i])
            con.pi[con.ij[2i - 1]] = con.ij[2i]
            con.c[con.ij[2i - 1]] = con.A[i]
            con.vaux[con.ij[2i - 1]] = con.A[i]
        end
    end

    #(@printf "PI=")
    #show(stdout, con.pi)
    #(@printf "\nC=")
    #show(stdout, con.c)
    #(@printf "\nvaux=")
    #show(stdout, con.vaux)
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

    #for i = 1:con.nnodes
    #    (@printf "%d:\n" i)
    #    (@printf "  piinv_succ=%d\n" con.piinv_succ[i])
    #    (@printf "  piinv_elem=%d\n" con.piinv_elem[i])
    #    (@printf "  piinv_last=%d\n" con.piinv_last[i])
    #    (@printf "  piinv_idx=%d\n" con.piinv_idx[i])
    #end
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
    #(@printf "Visiting from node %d color=%d\n" initialpoint color)
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
function value(con::Context, components::Int)
    initialpoint = 1
    color = 1

    # init_depth_first()
    con.visited .= false
    con.component .= 0

    while true
        visit_from(con, initialpoint, color);
        while ((initialpoint <= con.nnodes) && (con.component[initialpoint] != 0))
            initialpoint += 1
        end
        color += 1
        (initialpoint <= con.nnodes) || break
    end

    components = color - 1
end

# -----------------------------------------------------------------------------
function first_order_improvement(con::Context, improved::Bool)
    for i = 1:con.narcs
        if (con.chi[con.ij[2i]] > con.newchi[con.ij[2i - 1]])
            show_info_improve_chi(i)
            improved = true;
            con.newpi[con.ij[2i - 1]] = con.ij[2i]
            con.newchi[con.ij[2i - 1]] = con.chi[con.ij[2i]]
            con.newc[con.ij[2i - 1]] = con.A[i];
        end
    end
end

# -----------------------------------------------------------------------------
function second_order_improvement(con::Context, components::Int, improved::Bool, ϵ::Float64)
    if (components > 1)
        for i = 1:con.narcs
            if (con.chi[con.ij[2i]] == con.newchi[con.ij[2i - 1]])
                w = a[i] + v[con.ij[2i]] - con.chi[con.ij[2i - 1]]
                if (w > con.vaux[con.ij[2i - 1]] + ϵ)
                    show_info_improve_bias(con, i)
                    improved = true
                    con.vaux[ij[2i - 1]] = w
                    con.newpi[ij[2 - 1]] = con.ij[2i]
                    con.newc[ij[2i - 1]] = con.A[i]
                end
            end
        end
    else
        for i = 1:con.narcs
            w = con.A[i]+ con.v[con.ij[2i]] - con.chi[con.ij[2i - 1]]
            if (w > con.vaux[con.ij[2i - 1]] + ϵ)
                show_info_improve_bias(con, i)
                improved = true
                con.vaux[con.ij[2i - 1]] = w
                con.newpi[con.ij[2i - 1]] = con.ij[2i]
                con.newc[con.ij[2i - 1]] = con.A[i]
            end
        end
    end
end

# -----------------------------------------------------------------------------
function improve(con::Context, components::Int, improved::Bool, ϵ::Float64)
    improved = false

    con.newchi .= con.chi
    con.vaux .= con.v
    con.newpi .= con.pi
    con.newc .= con.c

    if (components > 1)
        first_order_improvement(con, improved, ϵ)
    end

    if (improved == false)
        second_order_improvement(con, components, improved, ϵ)
    end
end

# -----------------------------------------------------------------------------
"""
    howard(S::SpaMP)

TODO

# Examples
```julia-repl
julia> S = mpsparse([1 2; 3 4])
2×2 Max-Plus Sparse Matrix with 4 stored entries:
  1   2
  3   4

julia> l,v = howard(S)
(MP{Float64}[4, 4], MP{Float64}[2, 4])

julia> (A * v) == (l[1] * v)
true

julia> S = mpsparse([mp0 2 mp0; mp1 mp0 mp0; mp0 mp0 2])
3×3 Max-Plus Sparse Matrix with 3 stored entries:
  .   2   .
  0   .   .
  .   .   2

julia> l,v = howard(S)
(MP{Float64}[1, 1, 2], MP{Float64}[1, 0, 2])

TODO finish
```
"""
function howard(S::SpaMP)
    improved::Bool = false
    iteration::Int = 0
    components::Int = 0
    MAX_ITERATIONS::Int = 1000

    con = Context(plustimes(S))
    #(@printf "IJ=")
    #show(stdout, con.ij)
    #(@printf "\nA=")
    #show(stdout, con.A)
    #(@printf "\n")
    sanity_checks(con)
    initial_policy(con)
    ϵ = epsilon(con)
    build_inverse(con)

    while true
        value(con, components)
        #show_info(con, iteration)
        improve(con, components, improved, ϵ)
        update_policy(con)
        build_inverse(con)
        iteration += 1
        ((improved != false) && (iteration < MAX_ITERATIONS)) || break
    end

    if (iteration >= MAX_ITERATIONS)
        throw(AssertiobError("Max iteration reached"))
    end

    MP(con.chi), MP(con.v)
end

# -----------------------------------------------------------------------------
"""
"""
howard(M::ArrMP) = howard(mpsparse(M))
