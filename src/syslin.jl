# ==============================================================================

# Implicit dynamic linear maxplus system.
struct MPSysLin{T}
    A::ArrMP{T}
    B::ArrMP{T}
    C::ArrMP{T}
    D::ArrMP{T}
    x0::ArrMP{T}
end

# ==============================================================================

function Base.:(==)(x::MPSysLin, y::MPSysLin)
    (x.A == y.A) && (x.B == y.B) && (x.C == y.C) && (x.D == y.D) && (x.x0 == y.x0)
end

# ==============================================================================

function Base.show(io::IO, S::MPSysLin{T}) where T
    (@printf io "Implicit dynamic linear maxplus system:\n")
    (@printf io "  x(n) = D*x(n) + A*x(n-1) + B*u(n)\n  y(n) = C*x(n)\n  x(0) = x0\n\nwith:")
    (@printf io "\n  D  = "); Base.show(io, S.D)
    (@printf io "\n  A  = "); Base.show(io, S.A)
    (@printf io "\n  B  = "); Base.show(io, S.B)
    (@printf io "\n  C  = "); Base.show(io, S.C)
    (@printf io "\n  x0 = "); Base.show(io, S.x0)
    (@printf io "\n")
end

# ==============================================================================

function LaTeX(io::IO, S::MPSysLin)
    (@printf io "\\begin{array}{lcl}\n")
    (@printf io "x_n & = & ")
    LaTeX(io, S.D)
    (@printf io " x_n \\oplus ")
    LaTeX(io, S.A)
    (@printf io " x_{n-1} \\oplus ")
    LaTeX(io, S.B)
    (@printf io " u_n\\\\")
    (@printf io "y_n & = & ")
    LaTeX(io, S.C)
    (@printf io " x_n\\\\")
    (@printf io "x_0 & = & ")
    LaTeX(io, S.x0)
    (@printf io "\\end{array}")
end

# ==============================================================================
# Utility functions against vectors! Force seeing a vector like a matrix.

# Vector 2 matrix conversion
v2m(v::Vector) = reshape(v, length(v), 1)
v2m(v::Array) = v

# Fix Vector size returning a single param
size2(A::SpaMP) = (size(A, 1), size(A, 2))
size2(A::Array) = (size(A, 1), size(A, 2))

# ==============================================================================
# Dense

function mpsyslin(A::ArrMP{T}, B::ArrMP{T}, C::ArrMP{T}, D::ArrMP{T}, x0::ArrMP{T}) where T
    (ma,na) = size2(A)
    (ma != na) && error("Matrix A shall be squared")
    (mb,nb) = size2(B)
    ((mb != na) && (mb != 0)) && error("The row dimension of B does not agree with dimensions of A")
    (mc,nc) = size2(C)
    ((nc != na) && (mc != 0)) && error("The column dimension of C $(mc) $(nc) does not agree with dimensions of A $(ma) $(na)")
    (mx0,nx0) = size2(x0)
    ((mx0 != na) || (nx0 != min(na, 1))) && error("Dimensions of x0 do not agree with dimensions of A")
    (md,nd) = size2(D)
    ((md != na) || (nd != na)) && error("The column dimension of D does not agree with dimensions of A")
    MPSysLin(A, v2m(B), v2m(C), D, v2m(x0))
end

function mpsyslin(A::ArrMP{T}, B::ArrMP{T}, C::ArrMP{T}, D::ArrMP{T}) where T
    mpsyslin(A, B, C, D, mpfzeros(T, size(A,2), 1))
end

function mpsyslin(A::ArrMP{T}, B::ArrMP{T}, C::ArrMP{T}) where T
    na = size(A,2)
    mpsyslin(A, B, C, mpeye(T, na, na), mpfzeros(T, na, 1))
end

# ==============================================================================
# Sparse

function mpsyslin(A::SpaMP{T}, B::SpaMP{T}, C::SpaMP{T}, D::SpaMP{T}, x0::SpaMP{T}) where T
    mpsyslin(full(A), full(B), full(C), full(D), full(x0))
end

function mpsyslin(A::SpaMP{T}, B::SpaMP{T}, C::SpaMP{T}, D::SpaMP{T}) where T
    mpsyslin(full(A), full(B), full(C), full(D), mpfzeros(T, size(A,2), 1))
end

function mpsyslin(A::SpaMP{T}, B::SpaMP{T}, C::SpaMP{T}) where T
    na = size(A,2)
    mpsyslin(full(A), full(B), full(C), mpeye(T, na, na), mpfzeros(T, na, 1))
end

# ==============================================================================
# %mpls_a_mpls.sci
# Parallel composition

function Base.:(+)(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    MPSysLin([x.A mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.A],
             [x.B; y.B],
             [x.C y.C],
             [x.D mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_m_mpls.sci
# Series composition.

function Base.:(*)(y::MPSysLin{T}, x::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc2 = size(y.C, 1)
    MPSysLin([x.A mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.A],
             [x.B; mpfzeros(T, n2, nb1)],
             [mpfzeros(T, nc2, n1) y.C],
             [x.D mpfzeros(T, n1, n2); y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_g_mpls.sci
# Diagonal composition

function Base.:(|)(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nb2 = size(y.B, 2)
    nc1 = size(x.C, 1)
    nc2 = size(y.C, 1)
    MPSysLin([x.A mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.A],
             [x.B mpfzeros(T, n1, nb2); mpfzeros(T, n2, nb1) y.B],
             [x.C mpfzeros(T, nc1, n2); mpfzeros(T, nc2, n1) y.C],
             [x.D mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_m_mpls.sci

# ==============================================================================
# %mpls_f_mpls.sci
# computes [S1;S2] that is :  inputs in common,  concatenation of outputs

function Base.:vcat(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nc1 = size(x.C, 1)
    nc2 = size(y.C, 1)
    MPSysLin([x.A mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.A],
             [x.B; y.B],
             [x.C mpfzeros(T, nc1, n2); mpfzeros(T, nc2, n1) y.C],
             [x.D mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_c_mpls.sci
# computes [S1,S2] that is  concatenation of inputs, addition of outputs
# TODO [S1 S2] works but not [S1, S2]

function Base.:hcat(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nb2 = size(y.B, 2)
    MPSysLin([x.A mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.A],
             [x.B mpfzeros(T, n1, nb2); mpfzeros(T, n2, nb1) y.B],
             [x.C y.C],
             [x.D mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_v_mpls.sci
# Feedback composition: computes star(S1*S2)*S1 in state-space form.

function Base.:(/)(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc1 = size(y.C, 1)
    MPSysLin([x.A mpfzeros(T, n1, n2); mpfzeros(T, n2, n1) y.A],
             [x.B; mpfzeros(T, n2, nb1)],
             [x.C mpfzeros(T, nc1, n2)],
             [x.D x.B * y.C; y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_m_talg.sci
# S1 * k

function Base.:(*)(S::MPSysLin{T}, k::U) where {T,U}
    MPSysLin(S.A, S.B * k, S.C, S.D, S.x0)
end

function Base.:(*)(k::U, S::MPSysLin{T}) where {T,U}
    MPSysLin(S.A, S.B, S.C * k, S.D, S.x0)
end

# %mpls_m_talg.sci
function Base.:(*)(S::MPSysLin{T}, M::SpaMP{T}) where T
    MPSysLin(S.A, S.B * M, S.C, S.D, S.x0)
end

function Base.:(*)(M::SpaMP{T}, S::MPSysLin{T}) where T
    MPSysLin(S.A, M * S.B, S.C, S.D, S.x0)
end

# ==============================================================================
# mplssize.sci
function Base.:size(S::MPSysLin{T}) where T
    [size(S.B, 1), size(S.B, 2), size(S.C, 1)]
end

# ==============================================================================
# explicit.sci
#mpstar ko pour identity

function mpexplicit(S::MPSysLin{T}) where T
    ds = mpstar(S.D)
    bs = ds * S.B
    ac = [ds * S.A; S.C]
    zerocol = map(x -> Bool(x.λ), mpones(T, 1, size(ac, 1)) * (ac .!= mpzero(T)))
    keep = findall(zerocol[1,:])
    c = [S.C[i] for i in keep]
    mpsyslin([ac[i, j] for i in keep, j in keep],
             [bs[i] for i in keep],
             reshape(c, 1, length(c)))
end

# ==============================================================================

#TODO
"""
    mpsyslin(A::SpaMP{T}, x0::Vector{MP{T}}, k::Int64; history=false)

Compute states X of an autonomous linear max-plus system:
`x(n+1) = A x(n)`  for n = 0 .. k
where:
- `A` is a system matrix
- `x0` is a initial state vector,
- `k` is the number of iterations
- when history is set to true save all computed states, else return the last one.

# Examples
```julia-repl
julia> S1 = mpsyslin(MP([1.0 2; 3 4]), MP([0.0; 0]), MP([0.0 0]), mpeye(Float64, 2,2));
julia> mpsimul(S1, MP(1:10))
1×10 Array{MP{Float64},2}:
 1  5  9  13  17  21  25  29  33  37
julia> mpsimul(S1, MP(1:10), history=false)
1×1 Array{MP{Float64},2}:
 37
```
"""
function mpsimul(S::MPSysLin{T}, u::ArrMP{T}, history::Bool) where T
    x = S.x0
    k = size(u, 1)
    if history
        Y = mpones(T, size(S.C, 1), k)
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
            Y[:,i] = S.C * x
        end
        Y
    else
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
        end
        Y = S.C * x
        Y[:,1]
    end
end

function mpsimul(S::MPSysLin{T}, u::VecMP{U}, history::Bool) where {T,U}
    mpsimul(S, map(x -> MP(T(x.λ)), reshape(u, length(u), 1)), history)
end
