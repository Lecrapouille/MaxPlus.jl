# ==============================================================================

# Implicit dynamic linear maxplus system
struct MPSysLin{T}
 A::SpaMP{T}
 B::SpaMP{T}
 C::SpaMP{T}
 D::SpaMP{T}
 x0::SpaMP{T}
end

# ==============================================================================

function Base.show(io::IO, S::MPSysLin{T}) where T
    (@printf io "Implicit dynamic linear maxplus system:\n")
    (@printf io "  x(n) = D*x(n) + A*x(n-1) + B*u(n)\n  y(n) = C*x(n)\n  x(0) = x0\n\nwith:")
    (@printf io "\n  D  = "); Base.show(io, dense(S.D))
    (@printf io "\n  A  = "); Base.show(io, dense(S.A))
    (@printf io "\n  B  = "); Base.show(io, dense(S.B))
    (@printf io "\n  C  = "); Base.show(io, dense(S.C))
    (@printf io "\n  x0 = "); Base.show(io, dense(S.x0))
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
size2(A::SpaMP{T}) where T = (size(A, 1), size(A, 2))
size2(A::Array{T}) where T = (size(A, 1), size(A, 2))

# ==============================================================================
# Sparse

function mpsyslin(A::SpaMP{T}, B::SpaMP{T}, C::SpaMP{T}, D::SpaMP{T}, x0::SpaMP{T}) where T
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
    MPSysLin(A,B,C,D,x0)
end

# S4 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]), mpsparse(mpeye(Int64, 3,3)), mpzeros(Int64, 3,1))
# S1 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(mpeye(Int64, 2,2)), mpzeros(Int64, 2,1))

# ==============================================================================
# Sparse

function mpsyslin(A::SpaMP{T}, B::SpaMP{T}, C::SpaMP{T}, D::SpaMP{T}) where T
    mpsyslin(A, B, C, D, mpzeros(T, size(A,2), 1))
end

# S4 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]), mpsparse(mpeye(Int64, 3,3)))
# S1 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(mpeye(Int64, 2,2)))

# ==============================================================================
# Sparse

function mpsyslin(A::SpaMP{T}, B::SpaMP{T}, C::SpaMP{T}) where T
    na = size(A,2)
    mpsyslin(A, B, C, mpzeros(T, na, na), mpzeros(T, na, 1))
end

# S4 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]))
# S1 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]))

# ==============================================================================
# Dense

function mpsyslin(A::ArrMP{T}, B::ArrMP{T}, C::ArrMP{T}, D::ArrMP{T}, x0::ArrMP{T}) where T
    mpsyslin(mpsparse(A), mpsparse(B), mpsparse(C), mpsparse(D), mpsparse(x0))
end

# S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), mpeye(Int64, 3,3), MP([0; 0; 0]))
# S1 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), mpeye(Int64, 2,2), MP([0; 0]))

function mpsyslin(A::ArrMP{T}, B::ArrMP{T}, C::ArrMP{T}, D::ArrMP{T}) where T
    mpsyslin(mpsparse(A), mpsparse(B), mpsparse(C), mpsparse(D), mpzeros(T, size(A,2), 1))
end

# S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), mpeye(Int64, 3,3))
# S1 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), mpeye(Int64, 2,2))

function mpsyslin(A::ArrMP{T}, B::ArrMP{T}, C::ArrMP{T}) where T
    na = size(A,2)
    mpsyslin(mpsparse(A), mpsparse(B), mpsparse(C), mpzeros(T, na, na), mpzeros(T, na, 1))
end

# S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]))
# S1 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]))

# ==============================================================================
# %mpls_a_mpls.sci
# Parallel composition

function Base.:(+)(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    MPSysLin([x.A mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.A],
             [x.B; y.B],
             [x.C y.C],
             [x.D mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_m_mpls.sci
# Series composition
# FIXME param x and y inverted ?

function Base.:(*)(x::MPSysLin{T}, y::MPSysLin{T}) where T
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc2 = size(y.C, 1)
    MPSysLin([x.A mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.A],
             [x.B; mpzeros(T, n2, nb1)],
             [mpzeros(T, nc2, n1) y.C],
             [x.D mpzeros(T, n1, n2); y.B * x.C y.D],
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
    MPSysLin([x.A mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.A],
             [x.B mpzeros(T, n1, nb2); mpzeros(T, n2, nb1) y.B],
             [x.C mpzeros(T, nc1, n2); mpzeros(T, nc2, n1) y.C],
             [x.D mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.D],
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
    MPSysLin([x.A mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.A],
             [x.B; y.B],
             [x.C mpzeros(T, nc1, n2); mpzeros(T, nc2, n1) y.C],
             [x.D mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.D],
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
    MPSysLin([x.A mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.A],
             [x.B mpzeros(T, n1, nb2); mpzeros(T, n2, nb1) y.B],
             [x.C y.C],
             [x.D mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.D],
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
    MPSysLin([x.A mpzeros(T, n1, n2); mpzeros(T, n2, n1) y.A],
             [x.B; mpzeros(T, n2, nb1)],
             [x.C mpzeros(T, nc1, n2)],
             [x.D x.B * y.C; y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# %mpls_m_talg.sci
# S1 * k

function Base.:(*)(S::MPSysLin{T}, k::T) where T
    MPSysLin(S.A, S.B * k, S.C, S.D, S.x0)
end

function Base.:(*)(k::T, S::MPSysLin{T}) where T
    MPSysLin(S.A, S.B, S.C * MP(k), S.D, S.x0)
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
# %mpls_e.sci
# extraction

# ==============================================================================
function full(S::MPSysLin{T}) where T
    (full(S.A), full(S.B), full(S.C), full(S.D), full(S.x0))
end

# ==============================================================================
# explicit.sci
#mpstar ko pour identity

function mpexplicit(S::MPSysLin{T}) where T
    (A,B,C,D,_) = full(S)
    ds = mpstar(D)
    bs = ds * B
    ac = [ds * A; C]
    zerocol = map(x -> Bool(x.λ), mpones(T, 1, size(ac, 1)) * (ac .!= mpzero(T)))
    keep = findall(zerocol[1,:])
    c = [C[i] for i in keep]
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
julia> A = MP([1.0 2.0; 3.0 4.0])
x0 = MP([1.0; 2.0])
x1 = mpsyslin(A, x0, 1)
x2 = mpsyslin(A, x1, 1)
X = mpsyslin(A, x0, 2, history=true)
x1 == A * x0
x2 == A * x1
[x0 x1 x2] == X
```
"""
function mpsimul(S::MPSysLin{T}, u::ArrMP{T}; history=true) where T
    x = full(S.x0)
    k = size(u, 2)
    if history
        Y = mpones(T, size(C, 1), k+1)
        Y[:,1] = S.C * x
        for i = 1:k
            x = S.A * x + S.B * u[:,i]
            Y[:,i+1] = S.C * x
        end
        Y
    else
        for i = 1:k
            x = S.A * x + S.B * u[:,i]
        end
        S.C * x
    end
end

# u=MP([1 2 3])
