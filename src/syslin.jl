# ##############################################################################
# State space representation of Max-Plus linear systems.
# Creation of Max-Plus dynamical linear systems in implicit state form:
#    X(n) = D.X(n) ⨁ A.X(n-1) ⨁ B.U(n),
#    Y(n) = C.X(n)
# ##############################################################################

# ==============================================================================
# Utility function: return an unified size for a matrix or a vector
size2(A::AbstractVecOrMat) = (size(A, 1), size(A, 2))

# ==============================================================================
# Implicit dynamic linear Max-Plus system.
struct MPSysLin
    A::MPAbstractVecOrMat
    B::MPAbstractVecOrMat
    C::MPAbstractVecOrMat
    D::MPAbstractVecOrMat
    x0::MPAbstractVecOrMat

    # Constructor with Max-Plus matrices: A, B, C, D and x0
    function MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat, D::MPAbstractVecOrMat, x0::MPAbstractVecOrMat)
        (ma,na) = size2(A)
        (ma != na) &&
            error("Matrix A shall be squared")

        (mb,nb) = size2(B)
        (mb != na) &&
            error("The row dimension of B: $(mb) is not in accordance with dimensions of A: $(na)")

        (mc,nc) = size2(C)
        (nc != na) &&
            error("The column dimension of C: $(nc) is not in accordance with dimensions of A: $(na)")

        (mx0,nx0) = size2(x0)
        (mx0 != na) &&
            error("Dimensions of x0: $(mx0) do not in accordance with dimensions of A: $(na)")

        (md,nd) = size2(D)
        ((md != na) || (nd != na)) &&
            error("The column dimension of D: $(md)×$(nd) is not in accordance with dimensions of A: $(na)")

        new(A, B, C, D, x0)
    end

    # Constructor with Max-Plus matrices: A, B, C, D (and implicit x0 set as zeros)
    function MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat, D::MPAbstractVecOrMat)
        na = size(A,2)
        MPSysLin(A, B, C, D, spzeros(MP, na, 1))
    end

    # Constructor with Max-Plus matrices: A, B, C (and implicit D and x0 set as zeros)
    function MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat)
        na = size(A,2)
        MPSysLin(A, B, C, spzeros(MP, na, na), spzeros(MP, na, 1))
    end
end # MPSysLin

# ==============================================================================
function Base.:(==)(x::MPSysLin, y::MPSysLin)
    (x.A == y.A) && (x.B == y.B) && (x.C == y.C) && (x.D == y.D) && (x.x0 == y.x0)
end

# ==============================================================================
# From ScicosLab file: mplssize.sci
function Base.:size(S::MPSysLin)
    [size(S.B, 1), size(S.B, 2), size(S.C, 1)]
end

# ==============================================================================
# Parallel composition.
# From ScicosLab file: %mpls_a_mpls.sci

function Base.:(+)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; y.B],
             [x.C y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Series composition.
# From ScicosLab file: %mpls_m_mpls.sci

function Base.:(*)(y::MPSysLin, x::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc2 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; spzeros(MP, n2, nb1)],
             [spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Diagonal composition.
# From ScicosLab file: %mpls_g_mpls.sci

function Base.:(|)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nb2 = size(y.B, 2)
    nc1 = size(x.C, 1)
    nc2 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B spzeros(MP, n1, nb2); spzeros(MP, n2, nb1) y.B],
             [x.C spzeros(MP, nc1, n2); spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# computes [S1;S2] that is : inputs in common, concatenation of outputs.
# From ScicosLab file: %mpls_f_mpls.sci

function Base.:vcat(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nc1 = size(x.C, 1)
    nc2 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; y.B],
             [x.C spzeros(MP, nc1, n2); spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# computes [S1,S2] that is concatenation of inputs, addition of outputs
# From ScicosLab file: %mpls_c_mpls.sci

function Base.:hcat(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nb2 = size(y.B, 2)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B spzeros(MP, n1, nb2); spzeros(MP, n2, nb1) y.B],
             [x.C y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Feedback composition: computes star(S1*S2)*S1 in state-space form.
# From ScicosLab file: %mpls_v_mpls.sci

function Base.:(/)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc1 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; spzeros(MP, n2, nb1)],
             [x.C spzeros(MP, nc1, n2)],
             [x.D x.B * y.C; y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# From ScicosLab file: explicit.sci

function mpexplicit(S::MPSysLin)
    ds = mpstar(S.D)
    A = ds * S.A
    B = ds * S.B
    AC = [A; S.C]
    MPSysLin(AC, B, S.C,
             spzeros(MP, size(AC,1), size(AC,2)),
             S.x0)
end

# ==============================================================================
function mpsimul(S::MPSysLin, u::MPAbstractVecOrMat, history::Bool)
    x = S.x0
    k = size(u, 1)
    if history
        Y = ones(MP, k, size(S.C, 1))
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
            Y[i,:] = S.C * x
        end
        Y
    else
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
        end
        Y = S.C * x
        Y[1,:]
    end
end

# ==============================================================================
function mpsimul(S::MPSysLin, u::AbstractVector{MP}, history::Bool)
    mpsimul(S, map(x -> MP(x.λ), reshape(u, length(u), 1)), history)
end
