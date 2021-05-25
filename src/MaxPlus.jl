# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# In the way of the ScicosLab Max-Plus toolbox.
#
# The documentation of functions for the REPL are hiden in docstrings.jl
# ==============================================================================

module MaxPlus

using
    LinearAlgebra, SparseArrays, PrettyTables, Printf

export # Max-Plus core
    MP, SpaMP, SpvMP, ArrMP, VecMP, mpzero, mpone, mp0, mp1, mptop, ε, mpe, mpI,
    mpeye, mpzeros, mpones, full, dense, array, plustimes, inv, mpsparse_map,
    mptrace, mpnorm, mpastarb, mpstar, mpplus, howard, mp_change_display,
    mpshow, LaTeX

export # Max-Plus Linear system
    MPSysLin, mpsimul, mpexplicit

export # Max-Plus flowhop
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

# ==============================================================================
# Types

# Max-Plus type
struct MP <: Real
    λ::Float64

    function MP(x::Real)
        # Avoid pathological cases: MP(Nan) made with mp0 * mptop
        isnan(x) ? new(typemin(Float64)) : new(Float64(x))
    end
end

# Sparse matrix (classic algebra) shorter name
const Sparse{T,U} = SparseMatrixCSC{T,U}

# Sparse matrix (max-plus algebra) shorter name
const SpaMP{U} = SparseMatrixCSC{MP,U}

# Sparse vector (max-plus algebra) shorter name
const SpvMP{U} = SparseVector{MP,U}

# Dense matrix (max-plus algebra) shorter name
const ArrMP{N} = Array{MP,N}

# Dense vector (max-plus algebra) shorter name
const VecMP{N} = Vector{MP}

# ==============================================================================

include("fallbacks.jl") # Fixes for Julia regressions
include("show.jl") # Display scalars and matrices

# ==============================================================================
# Julia type promotions and conversions to/from Max-Plus number

Base.promote_rule(::Type{MP}, ::Type{T}) where T = MP
Base.convert(::MP, x::Number) = MP(x)
Base.float(x::MP) = x.λ

# ==============================================================================
# Constructors

# Copy constructor
MP(x::MP) = MP(x.λ)

# Constructor needed for I operator building identity matrices
MP(x::Bool) = x ? mpone() : mpzero()

# Constructor from non Max-Plus dense matrix
MP(A::Array) = map(MP, A)

# Constructor from non Max-Plus sparse matrix
MP(S::Sparse{T,U}) where {T,U} = convert(SpaMP{U}, S)

# Constructor like SparseArrays.sparse
MP(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {Tv,Ti<:Integer} = sparse(I, J, MP(V))

# Constructor from non Max-Plus sparse vector
MP(V::SparseVector{T,U}) where {T, U} = convert(SparseVector{MP,U}, V)

# Constructor from non Max-Plus range
MP(x::UnitRange) = Vector{MP}(x)
MP(x::StepRangeLen) = Vector{MP}(x)

# ==============================================================================
# Algebra redefinition for Max-Plus: zero and one

Base.zero(::Type{MP}) = MP(typemin(Float64))
Base.zero(x::MP) = zero(typeof(x))

Base.one(::Type{MP}) = MP(zero(Float64))
Base.one(x::MP) = one(typeof(x))

mpzero() = MP(typemin(Float64))
mpone() = MP(zero(Float64))

# ==============================================================================
# Max-Plus and Min-Plus constants

const global mp0 = mpzero()
const global ε = mp0
const global mp1 = mpone()
const global mpe = mp1
const global mptop = MP(typemax(Float64))

# ==============================================================================
# Max-Plus core overriden operators

Base.:(+)(x::MP, y::MP) = MP(max(x.λ, y.λ))
Base.:(+)(x::MP, y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP) = MP(max(x, y.λ))

Base.:(*)(x::MP, y::MP) = MP(x.λ + y.λ)
Base.:(*)(x::MP, y::Real) = MP(x.λ + y)
Base.:(*)(x::Real, y::MP) = MP(x + y.λ)

Base.:(/)(x::MP, y::MP) = MP(x.λ - y.λ)
Base.:(/)(x::MP, y::Real) = MP(x.λ - b)
Base.:(/)(x::Real, y::MP) = MP(a - y.λ)

Base.:(-)(x::MP, y::MP) = MP(x.λ - y.λ)
Base.:(-)(x::MP, y::Real) = ((x == mp0) && (y == typemin(Real))) ? mp0 : MP(x.λ - y)
Base.:(-)(x::Real, y::MP) = ((x == typemin(Real)) && (y == mp0)) ? mp0 : MP(x - y.λ)
Base.:(-)(x::MP) = (x == mpzero()) ? mpzero() : MP(-x.λ)
Base.sign(x::MP) = Base.sign(x.λ)

function Base.inv(A::ArrMP)
    isempty(A) && return A

    Z = transpose(-A)

    # Checks
    (A * Z != mpeye(size(A,1), size(Z,2))) && error("The matrix cannot be inversed")
    if (size(A,1) == size(A,2))
        (Z * A != mpeye(size(Z,1), size(A,2))) && error("The matrix cannot be inversed")
    end

    return Z
end

Base.:(\)(A::ArrMP, b::ArrMP) = inv(A) * b

# Needed to accept MP(1)^0.5
Base.:(^)(x::MP, y::Float64) = MP(x.λ * y)

Base.:min(x::MP, y::MP) = (x * y) / (x + y)
Base.:min(x::MP, y::Real) = min(x, MP(y))
Base.:min(x::Real, y::MP) = min(MP(x), y)
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

# ==============================================================================
# Comparator

Base.:(<=)(x::MP,   y::MP)   = (x.λ <= y.λ)
Base.:(<=)(x::MP,   y::Real) = (x.λ <= y)
Base.:(<=)(x::Real, y::MP)   = (x   <= y.λ)

Base.:(<)(x::MP,   y::MP)   = (x.λ < y.λ)
Base.:(<)(x::MP,   y::Real) = (x.λ < y)
Base.:(<)(x::Real, y::MP)   = (x   < y.λ)

# ==============================================================================
# Algebra conversion: Max-Plus to Min-Plus or Max-Plus to classic algebra

plustimes(n::MP) = n.λ
plustimes(A::ArrMP) = map(x -> x.λ, A)
plustimes(S::SpaMP) = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(x -> x.λ, S.nzval))
full(S::SpaMP) = Matrix(S)
dense(S::SpaMP) = Matrix(S)
array(S::SpaMP) = Matrix(map(x -> x.λ, S))

# ==============================================================================
# Matrices constructions

const global mpI = UniformScaling(mpone())

# Identity dense matrices
mpeye(n::Int64) = Matrix{MP}(mpI, n, n)
mpeye(m::Int64, n::Int64) = Matrix{MP}(mpI, m, n)
mpeye(A::Array) = Matrix{MP}(mpI, size(A,1), size(A,2))
mpeye(A::ArrMP) = Matrix{MP}(mpI, size(A,1), size(A,2))

# Zero sparse matrices
mpzeros(n::Int64) = spzeros(MP, n)
mpzeros(m::Int64, n::Int64) = spzeros(MP, m, n)
mpzeros(A::Array) = spzeros(MP, size(A,1), size(A,2))
mpzeros(A::ArrMP) = spzeros(MP, size(A,1), size(A,2))

# One dense matrices
mpones(n::Int64) = ones(MP, n)
mpones(m::Int64, n::Int64) = ones(MP, m, n)
mpones(A::Array) = mpones(size(A,1), size(A,2))
mpones(A::ArrMP) = mpones(size(A,1), size(A,2))

# ==============================================================================
# These functions fix Julia bugs

Base.abs(x::MP) = MP(abs(x.λ))
Base.abs2(x::MP) = x.λ + x.λ
Base.round(x::MP, n::RoundingMode) = MP(round(x.λ, n))
Base.big(x::MP) = Base.big(x.λ)

# ==============================================================================
# Max-Plus matrices operations

# Map function f to a sparse matrice
function mpsparse_map(f, M::SpaMP)
    SparseMatrixCSC(M.m, M.n, M.colptr, M.rowval, map(f, M.nzval))
end

# Eigenvalues and eigenvectors
include("howard.jl")

# Trace
mptrace(A::ArrMP) = isempty(A) ? mp0 : sum(diag(A))
mptrace(A::Array) = isempty(A) ? mp0 : sum(MP(diag(A)))
mptrace(S::SpaMP) = isempty(S) ? mp0 : sum(diag(S))
mptrace(S::Sparse) = sum(MP(diag(S)))

# Norm
mpnorm(A::Array) = MP(A[argmax(A)].λ - A[argmin(A)].λ)
mpnorm(S::SpaMP) = MP(S[argmax(S)].λ - S[argmin(S)].λ)

# A^+
function mpplus(A::ArrMP)
    n = size(A, 1)
    (n != size(A, 2)) && error("Matrix shall be squared")
    C = A
    for k in 1:n
        t = (C[k,k].λ <= zero(Float64)) ? zero(Float64) : typemax(Float64);
        for j in 1:n, i in 1:n
            C[i,j] = MP(max(C[i,j].λ, C[i,k].λ + C[k,j].λ + t))
        end
    end
    C
end

mpplus(x::MP) = mpplus([x])[1,1]

# A^*
mpstar(A::ArrMP) = mpeye(A) + mpplus(A)
mpstar(x::MP) = mpstar([x])[1,1]

# A^* b
mpastarb(A::ArrMP, b::ArrMP) = mpstar(A) * b

# A^* and A^+
Base.:(^)(A::Matrix{MP},::typeof(*)) = mpplus(A)
Base.:(^)(A::Matrix{MP},::typeof(+)) = mpstar(A)

# ==============================================================================
# Max-Plus domain specific toolboxes

include("syslin.jl") # Dynamic Linear system manipulation helpers
#include("flowshop.jl")

# ==============================================================================
# Max-Plus REPL documentation

include("docstrings.jl")

end # MaxPlus module
