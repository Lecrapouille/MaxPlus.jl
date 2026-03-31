################################################################################
# State-space representation for Max-Plus linear systems
# Max-Plus implicit state-space form:
#     X(n) = D*X(n) ⨁ A*X(n-1) ⨁ B*U(n)
#     Y(n) = C*X(n)
################################################################################

# ==============================================================================
# Utility functions: normalize vectors/matrices for consistent 2D shapes

function _col(v::SparseVector{T}) where {T<:Tropical}
    n = length(v)
    I = SparseArrays.nonzeroinds(v)
    V = SparseArrays.nonzeros(v)
    SparseArrays.sparse(I, fill(1, length(I)), V, n, 1)
end
_col(v::AbstractVector{<:Tropical}) = reshape(v, :, 1)
_col(M::AbstractMatrix{<:Tropical}) = M

function _row(v::SparseVector{T}) where {T<:Tropical}
    n = length(v)
    I = SparseArrays.nonzeroinds(v)
    V = SparseArrays.nonzeros(v)
    SparseArrays.sparse(fill(1, length(I)), I, V, 1, n)
end
_row(v::AbstractVector{<:Tropical}) = reshape(v, 1, :)
_row(M::AbstractMatrix{<:Tropical}) = M
_size2(A::AbstractVecOrMat) = (size(A, 1), size(A, 2))

# Use sparse D / x0 when building from sparse data (flowshop, mpshift, …).
function _symlin_want_sparse(As::MPAbstractVecOrMat...)
    return any(A -> A isa SparseMatrixCSC{MP} || A isa SparseVector{MP}, As)
end

# ==============================================================================
# Main struct for Max-Plus System Linear (MPSysLin)

struct MPSysLin
    A::AbstractMatrix{MP}
    B::AbstractMatrix{MP}
    C::AbstractMatrix{MP}
    D::AbstractMatrix{MP}
    x0::AbstractMatrix{MP}

    # Constructor: A, B, C, D, x0 (all Max-Plus matrices or vectors)
    function MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat, D::MPAbstractVecOrMat, x0::MPAbstractVecOrMat)
        Am = _col(A)
        Bm = _col(B)
        Cm = _row(C)
        Dm = _col(D)
        x0m = _col(x0)

        (ma, na) = _size2(Am)
        if ma != na
            error("Matrix A must be square (found $(ma)x$(na))")
        end

        (mb, nb) = _size2(Bm)
        if mb != na
            error("B has $(mb) rows, but should match A dimension $(na)")
        end

        (mc, nc) = _size2(Cm)
        if nc != na
            error("C has $(nc) columns, but should match A dimension $(na)")
        end

        (mx0, nx0) = _size2(x0m)
        if mx0 != na
            error("x0 has $(mx0) rows, but should match A dimension $(na)")
        end

        (md, nd) = _size2(Dm)
        if md != na || nd != na
            error("D has dimension $(md)x$(nd), should be square and match A ($(na)x$(na))")
        end

        new(Am, Bm, Cm, Dm, x0m)
    end

    # Constructor: A, B, C, D (x0 defaults to zeros / sparse zeros)
    function MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat, D::MPAbstractVecOrMat)
        na = size(A, 2)
        x0 = _symlin_want_sparse(A, B, C, D) ? spzeros(MP, na, 1) : zeros(MP, na, 1)
        MPSysLin(A, B, C, D, x0)
    end

    # Constructor: A, B, C (D defaults to identity, x0 defaults to zeros / sparse zeros)
    function MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat)
        na = size(A, 2)
        sp = _symlin_want_sparse(A, B, C)
        Ddef = sp ? speye(MP, na, na) : eye(MP, na, na)
        x0 = sp ? spzeros(MP, na, 1) : zeros(MP, na, 1)
        MPSysLin(A, B, C, Ddef, x0)
    end
end # struct MPSysLin

# ==============================================================================
# Equality between two MPSysLin systems

function Base.:(==)(x::MPSysLin, y::MPSysLin)
    x.A == y.A && x.B == y.B && x.C == y.C && x.D == y.D && x.x0 == y.x0
end

# ==============================================================================
# Return size info [number of states, inputs, outputs]

function Base.size(S::MPSysLin)
    [size(S.B, 1), size(S.B, 2), size(S.C, 1)]
end

# ==============================================================================
# Scalar ("gain") multiplication, Max-plus style

Base.:(*)(a::MP, S::MPSysLin) = MPSysLin(S.A, S.B, a .* S.C, S.D, S.x0)
Base.:(*)(S::MPSysLin, a::MP) = MPSysLin(S.A, a .* S.B, S.C, S.D, S.x0)
Base.:(*)(a::Real, S::MPSysLin) = MP(a) * S
Base.:(*)(S::MPSysLin, a::Real) = S * MP(a)

# ==============================================================================
# Parallel composition. S1 + S2
# From ScicosLab file: %mpls_a_mpls.sci

function Base.:+(x::MPSysLin, y::MPSysLin)
    n1, n2 = size(x.A, 1), size(y.A, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; y.B],
             [x.C y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Series composition. S1 * S2
# From ScicosLab file: %mpls_m_mpls.sci

function Base.:(*)(y::MPSysLin, x::MPSysLin)
    n1, n2 = size(y.A, 1), size(x.A, 1)
    nb1 = size(y.B, 2)
    nc2 = size(x.C, 1)
    MPSysLin([y.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) x.A],
             [y.B; spzeros(MP, n2, nb1)],
             [spzeros(MP, nc2, n1) x.C],
             [y.D spzeros(MP, n1, n2); x.B * y.C x.D],
             [y.x0; x.x0])
end

# ==============================================================================
# Diagonal composition. S1 | S2
# From ScicosLab file: %mpls_g_mpls.sci

function Base.:|(x::MPSysLin, y::MPSysLin)
    n1, n2 = size(x.A, 1), size(y.A, 1)
    nb1, nb2 = size(x.B, 2), size(y.B, 2)
    nc1, nc2 = size(x.C, 1), size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B spzeros(MP, n1, nb2); spzeros(MP, n2, nb1) y.B],
             [x.C spzeros(MP, nc1, n2); spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Inputs in common. [S1; S2]
# From ScicosLab file: %mpls_f_mpls.sci

function Base.vcat(x::MPSysLin, y::MPSysLin)
    n1, n2 = size(x.A, 1), size(y.A, 1)
    nc1, nc2 = size(x.C, 1), size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; y.B],
             [x.C spzeros(MP, nc1, n2); spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Outputs addition. [S1 S2]
# From ScicosLab file: %mpls_c_mpls.sci

function Base.hcat(x::MPSysLin, y::MPSysLin)
    n1, n2 = size(x.A, 1), size(y.A, 1)
    nb1, nb2 = size(x.B, 2), size(y.B, 2)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B spzeros(MP, n1, nb2); spzeros(MP, n2, nb1) y.B],
             [x.C y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Feedback composition. S1 / S2. computes star(S1*S2)*S1 in state-space form.
# From ScicosLab file: %mpls_v_mpls.sci

function Base.:/(x::MPSysLin, y::MPSysLin)
    n1, n2 = size(x.A, 1), size(y.A, 1)
    nb1 = size(x.B, 2)
    nc1 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; spzeros(MP, n2, nb1)],
             [x.C spzeros(MP, nc1, n2)],
             [x.D x.B * y.C; y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# From ScicosLab file: explicit.sci  (sélection des états utiles via [A; C])

function explicit(S::MPSysLin)
    ds = star(S.D)
    a = ds * S.A
    b = ds * S.B
    c = S.C
    ac = [a; c]
    n = size(a, 2)
    keep = Int[j for j in 1:n if any(r -> ac[r, j] != mp0, 1:size(ac, 1))]
    if isempty(keep)
        error("explicit: no active state (all columns of [star(D)*A; C] are zero)")
    end
    Ak = a[keep, keep]
    Bk = b[keep, :]
    Ck = c[:, keep]
    nk = length(keep)
    Dk = spzeros(MP, nk, nk)
    x0k = S.x0[keep, :]
    MPSysLin(Ak, Bk, Ck, Dk, x0k)
end

# ==============================================================================
function implicit(S::MPSysLin)
    n = size(S.A, 1)
    ds = star(S.D)
    A2 = ds * S.A
    B2 = ds * S.B
    sp = _symlin_want_sparse(S.A, S.B, S.C, S.D)
    Dnew = sp ? speye(MP, n, n) : eye(MP, n, n)
    return MPSysLin(A2, B2, S.C, Dnew, S.x0)
end

# ==============================================================================
# Simulation: Recurrence for x = A*x + B*u (explicit form) or X = D*X + A*x_prev + B*u

function simul(S::MPSysLin, u::AbstractMatrix{MP}, history::Bool)
    ds = star(S.D)
    Aex = ds * S.A
    Bex = ds * S.B
    x = copy(S.x0)
    k = size(u, 1)
    ny = size(S.C, 1)
    if history
        Y = Matrix{MP}(undef, k, ny)
        @inbounds for i = 1:k
            ui = reshape(u[i, :], :, 1)
            x = Aex * x + Bex * ui
            yi = S.C * x
            Y[i, :] .= vec(yi)
        end
        return Y
    else
        @inbounds for i = 1:k
            ui = reshape(u[i, :], :, 1)
            x = Aex * x + Bex * ui
        end
        return vec(S.C * x)
    end
end

# ==============================================================================
# Overload for input u as vector (convert to matrix)

function simul(S::MPSysLin, u::AbstractVector{MP}, history::Bool)
    simul(S, reshape(u, length(u), 1), history)
end

# ==============================================================================
# mpfull / mpsparse : docstrings dans docstrings/syslin.jl

full(S::MPSysLin) = MPSysLin(dense(S.A), dense(S.B), dense(S.C), dense(S.D), dense(S.x0))
SparseArrays.sparse(S::MPSysLin) = MPSysLin(SparseArrays.sparse(S.A), SparseArrays.sparse(S.B), SparseArrays.sparse(S.C), SparseArrays.sparse(S.D), SparseArrays.sparse(S.x0))

# Aliases for Scilab compatibility
const mpfull = full
const mpsparse = SparseArrays.sparse