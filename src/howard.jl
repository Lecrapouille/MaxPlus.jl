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
    (IJ, A, nnodes, narcs) = spget(S::SpaMP)

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
        new(IJ, A, nnodes, narcs,
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
function show_info_improve_chi(con::Context, i::Int)
    I = con.ij[2i-1]
    J = con.ij[2i]
    (@printf "Improvement of the CYCLE TIME at node %d\n" I)
    (@printf "arc %d: %d--->%d chi[%d]-chi[%d]=%f-%f=%g>0\n" i I J I J con.chi[J] con.chi[I] con.chi[J]-con.chi[I])
end

# -----------------------------------------------------------------------------
function show_info_improve_bias(con::Context, i::Int)
    #(@printf "type 2 improvement\n")
    I = con.ij[2i-1]
    J = con.ij[2i]
    #(@printf "Improvement of the BIAS at node %d\n" I)
    #(@printf "A[%d] + v[%d] - chi[%d] - v[%d]= %f + %f - %f - %f = %f > 0\n" i J I I con.A[i] con.v[J] con.chi[I] con.v[I] con.A[i]+con.v[J]-con.chi[I]-con.v[I])
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
        #(@printf " => ij[%d]: %d\n"  2i-1  con.ij[2i-1])
        u[con.ij[2i-1]] = true
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
        #(@printf "\ncon.vaux[%d]=%f A[%d]=%f\n" con.ij[2i-1] con.vaux[con.ij[2i-1]] i con.A[i])
        if (con.vaux[con.ij[2i-1]] <= con.A[i])
            #(@printf "=> pi[ij[2i-1]] = ij[2i] => pi[%d] = %d\n" con.ij[2i-1] con.ij[2i])
            con.pi[con.ij[2i-1]] = con.ij[2i]
            con.c[con.ij[2i-1]] = con.A[i]
            con.vaux[con.ij[2i-1]] = con.A[i]
            #(@printf "=> c[ij[2i-1]] = A[i] => c[%d] = %f\n" con.ij[2i-1] con.A[i])
            #(@printf "=> vaux[ij[2i-1]] = A[i] => vaux[%d] = %f\n" con.ij[2i-1] con.A[i])
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
            #show_info_improve_chi(con, i)
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
                    #show_info_improve_bias(con, i)
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
                #show_info_improve_bias(con, i)
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
"""
    howard(S::SpaMP)

TODO

# Example 1: Irreducible matrix, only 1 eigenvalue.
```julia-repl
julia> using SparseArrays

julia> S = sparse(MP([1 2; 3 4]))
2×2 Max-Plus Sparse Matrix with 4 stored entries:
  1   2
  3   4

julia> λ,v = howard(S)
(MP[4, 4], MP[2, 4])

# λ is constant
julia> (S * v) == (λ[1] * v)
true
```

# Example 2: Two blocks diagonal matrix.
```julia-repl
julia> using SparseArrays, LinearAlgebra

julia> S = sparse([mp0 2 mp0; mp1 mp0 mp0; mp0 mp0 2])
3×3 Max-Plus Sparse Matrix with 3 stored entries:
  .   2   .
  0   .   .
  .   .   2

julia> λ,v = howard(S)
(MP[1, 1, 2], MP[1, 0, 2])

# The entries of λ take two values
julia> (S / Matrix(Diagonal(λ))) * v == v
true
```

# Example 3: Block triangular matrix with 2 eigenvalues.
```julia-repl
julia> S = sparse([1 1; mp0 2])
2×2 Max-Plus sparse matrix with 3 stored entries:
  1   1
  .   2

julia> λ,v = howard(S)
(MP[2, 2], MP[1, 2])

julia> (S * v) == (λ[1] * v)
true

# But MP(1) is also eigen value
S * [0; %0] == MP(1) * [0; mp0]
```

# Example 4: Block triangular matrix with 1 eigenvalue
```julia-repl
julia> using SparseArrays, LinearAlgebra

julia> S = sparse([2 1; mp0 mp1])
2×2 Max-Plus sparse matrix with 3 stored entries:
  2   1
  .   0

julia> λ,v = howard(S)
(MP[2, 0], MP[2, 0])

# λ is not constant λ[1] is eigen value
# with eigen vector [v(1);0]
julia> (S / Matrix(Diagonal(λ))) * v == v
true
```
"""
function howard(S::SparseMatrixCSC{MP}) where MP
    improved::Bool = false
    iterations::Int = 0
    components::Int = 0
    MAX_ITERATIONS::Int = 1000

    con = Context(S)
    #(@printf "IJ=")
    #show(stdout, con.ij)
    #(@printf "\nA=")
    #show(stdout, con.A)
    #(@printf "\n")
    sanity_checks(con)
    ϵ = epsilon(con)
    initial_policy(con)
    build_inverse(con)

    while true
        components = value(con)
        #show_info(con, iterations)
        improved = improve(con, components, ϵ)
        update_policy(con)
        build_inverse(con)
        iterations += 1
        ((improved != false) && (iterations < MAX_ITERATIONS)) || break
    end

    if (iterations >= MAX_ITERATIONS)
        throw(AssertionError("Max iterations reached"))
    end

    MP(con.chi), MP(con.v)
end
