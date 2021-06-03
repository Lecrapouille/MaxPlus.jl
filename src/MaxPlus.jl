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
    Trop,
    MP, SpaMP, SpvMP, ArrMP, VecMP,
    MI, SpaMI, SpvMI, ArrMI, VecMI,
    mp0, mp1, mptop, ε, mpe, mpI,
    mpeye, mpzeros, mpones, full, dense, array, plustimes, inv, mpsparse_map,
    mptrace, mpnorm, mpastarb, mpstar, mpplus, howard, mp_change_display,
    mpshow, LaTeX

export # Max-Plus Linear system
    MPSysLin, mpsimul, mpexplicit

export # Max-Plus flowhop
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

# ==============================================================================
# Tropical types

# Fake templates to make difference between Min-Plus and Max-Plus numbers
struct Min end
struct Max end
const MM = Union{Min, Max}

# Base class for Max-Plus and Min-Plus structures
struct Trop{T <: MM} <: Number
    λ::Float64
end

# Max-Plus structure
const MP = Trop{Max}

# Min-Plus structure
const MI = Trop{Min}

# Max-Plus or Min-Plus structure
const Tropical = Union{MP, MI}

# Sparse matrix (classic algebra) shorter name
const SparseMatrix{T,U} = SparseMatrixCSC{T,U}

# Sparse matrix (max-plus algebra) shorter name
const SpaTrop{U} = SparseMatrix{Tropical,U}
const SpaMP{U} = SparseMatrix{MP,U}
const SpaMI{U} = SparseMatrix{MI,U}

# Sparse vector (max-plus algebra) shorter name
const SpvTrop{U} = SparseVector{Tropical,U}
const SpvMP{U} = SparseVector{MP,U}
const SpvMI{U} = SparseVector{MI,U}

# Dense matrix (max-plus algebra) shorter name
const ArrTrop{N} = Array{Tropical,N}
const ArrMP{N} = Array{MP,N}
const ArrMI{N} = Array{MI,N}

# Dense vector (max-plus algebra) shorter name
const VecTrop{N} = Vector{Tropical}
const VecMP{N} = Vector{MP}
const VecMI{N} = Vector{MI}

# ==============================================================================

include("fallbacks.jl") # Fixes for Julia regressions
include("show.jl") # Display scalars and matrices

# ==============================================================================
# Julia type promotions and conversions to/from Max-Plus number

Base.promote_rule(::Type{Tropical}, ::Type{T}) where T = Tropical

Base.convert(::Trop{T}, x::Number) where T <: MM = Trop{T}(x)
Base.convert(::Type{MI}, x::MP) = MI(x.λ)
Base.convert(::Type{MP}, x::MI) = MP(x.λ)

Base.float(x::Trop) = x.λ

# ==============================================================================
# Constructors

# Copy constructor (i.e. MP(MI(42)))
Trop{T}(x::Trop) where T <: MM = Trop{T}(x.λ)

# Fake constructor with Not-A-Number check
Trop(::Type{T}, n::Float64) where T <: MM = isnan(n) ? zero(Trop{T}) : Trop{T}(n)
MPnan(n::Float64) = Trop(Max, n)
MInan(n::Float64) = Trop(Min, n)

# Constructor needed for I operator building identity matrices
MP(x::Bool) = x ? one(MP) : zero(MP)
MI(x::Bool) = x ? one(MI) : zero(MI)

# Constructor from non Max-Plus dense matrix
Trop{T}(A::Array) where T <: MM = map(Trop{T}, A)

# Constructor from non Max-Plus sparse matrix
Trop{T}(S::SparseMatrix{T1,U}) where {T<:MM,T1,U} = convert(SpaMP{U}, S)

# Constructor like SparseArrays.sparse
Trop{T}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MM,Tv,Ti<:Integer} = sparse(I, J, Trop{T}(V))

# Constructor from non Max-Plus sparse vector
Trop{T}(V::SparseVector{T1,U}) where {T<:MM,T1, U} = convert(SparseVector{Trop{T},U}, V)

# Constructor from non Max-Plus range
Trop{T}(x::UnitRange) where T <: MM = Vector{Trop{T}}(x)
Trop{T}(x::StepRangeLen) where T <: MM = Vector{Trop{T}}(x)

# ==============================================================================
# Algebra redefinition for Max-Plus: zero and one

Base.zero(::Type{MP}) = MP(typemin(Float64))
Base.zero(x::MP) = zero(typeof(x))
Base.one(::Type{MP}) = MP(zero(Float64))
Base.one(x::MP) = one(typeof(x))

# ==============================================================================
# Algebra redefinition for Min-Plus: zero and one

Base.zero(::Type{MI}) = MI(typemax(Float64))
Base.zero(x::MI) = zero(typeof(x))
Base.one(::Type{MI}) = MI(zero(Float64))
Base.one(x::MI) = one(typeof(x))

# ==============================================================================
# Max-Plus and Min-Plus constants

#TODO
const global mp0 = zero(MP)
const global ε = mp0
const global mp1 = one(MP)
const global mpe = mp1
const global mptop = MI(typemax(Float64))

# ==============================================================================
# Julia already manage it but I prefer to force it
# to be sure it will use the best algorithm
Base.:(^)(x::Trop{T}, y::Real) where T <: MM = Trop{T}(x.λ * y)

# ==============================================================================
# Max-Plus core overriden operators

Base.:(+)(x::MP, y::MP) = MP(max(x.λ, y.λ))
Base.:(+)(x::MP, y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP) = MP(max(x, y.λ))

Base.:(*)(x::MP, y::MP) = MPnan(x.λ + y.λ)
Base.:(*)(x::MP, y::Real) = MPnan(x.λ + y)
Base.:(*)(x::Real, y::MP) = MPnan(x + y.λ)

Base.:(+)(x::MI, y::MI) = MI(min(x.λ, y.λ))
Base.:(+)(x::MI, y::Real) = MI(min(x.λ, y))
Base.:(+)(x::Real, y::MI) = MI(min(x, y.λ))

Base.:(*)(x::MI, y::MI) = MInan(x.λ + y.λ)
Base.:(*)(x::MI, y::Real) = MInan(x.λ + y)
Base.:(*)(x::Real, y::MI) = MInan(x + y.λ)

#Base.:(-)(x::MP, y::MP) = MP(x.λ - y.λ)
#Base.:(-)(x::MP, y::Real) = ((x == mp0) && (y == typemin(Real))) ? mp0 : MP(x.λ - y)
#Base.:(-)(x::Real, y::MP) = ((x == typemin(Real)) && (y == mp0)) ? mp0 : MP(x - y.λ)
#Base.:(-)(x::MP) = (x == zero(MP)()) ? zero(MP)() : MP(-x.λ)
#Base.sign(x::MP) = Base.sign(x.λ)

# ==============================================================================
# Residu Max-Plus

Base.inv(A::ArrMP) = MI(transpose((map(x -> -x.λ, A))))

Base.:(/)(x::MP, y::MP) = MPnan(x.λ - y.λ)
Base.:(/)(A::ArrMP, B::ArrMP) = MP(inv(A) * MI(B))
Base.:(\)(x::MP, y::MP) = MPnan(y.λ - x.λ)
Base.:(\)(A::ArrMP, B::ArrMP) = MP(MI(A) * inv(B))

# ==============================================================================
# 

Base.:min(x::MP, y::MP) = (x * y) / (x + y)
Base.:min(x::MP, y::Real) = min(x, MP(y))
Base.:min(x::Real, y::MP) = min(MP(x), y)
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

# ==============================================================================
# Comparator

Base.:(<=)(x::Trop, y::Trop) = (x.λ <= y.λ)
Base.:(<=)(x::Trop, y::Real) = (x.λ <= y)
Base.:(<=)(x::Real, y::Trop) = (x   <= y.λ)

Base.:(<)(x::Trop, y::Trop) = (x.λ < y.λ)
Base.:(<)(x::Trop, y::Real) = (x.λ < y)
Base.:(<)(x::Real, y::Trop) = (x   < y.λ)

# ==============================================================================
# Algebra conversion: Max-Plus to Min-Plus or Max-Plus to classic algebra

plustimes(n::Trop) = n.λ
plustimes(A::ArrMP) = map(x -> x.λ, A)
plustimes(S::SpaMP) = SparseMatrix(S.m, S.n, S.colptr, S.rowval, map(x -> x.λ, S.nzval))
full(S::SpaMP) = Matrix(S)
dense(S::SpaMP) = Matrix(S)
array(S::SpaMP) = Matrix(map(x -> x.λ, S))

# ==============================================================================
# Matrices constructions

const global mpI = UniformScaling(one(MP))

# Identity dense matrices
mpeye(::Type{Trop{T}}, n::Int64) where T <: MM = Matrix{Trop{T}}(I, n, n)
mpeye(::Type{Trop{T}}, m::Int64, n::Int64) where T <: MM = Matrix{Trop{T}}(I, m, n)

# Sparse matrices of Max-Plus zeros
const DimOrInd = Union{Integer, AbstractUnitRange}
mpzeros(::Type{Trop{T}}, dims::DimOrInd...) where T <: MM = spzeros(Trop{T}, dims...)
mpzeros(::Type{Trop{T}}, dims::Tuple{Vararg{DimOrInd}}) where T <: MM  = spzeros(Trop{T}, dims)
mpzeros(::Type{Trop{T}}, dims::NTuple{N, Union{Integer, Base.OneTo}}) where {T<:MM,N}  = spzeros(Trop{T}, dims)
mpzeros(::Type{Trop{T}}, dims::NTuple{N, Integer}) where {T<:MM,N}  = spzeros(Trop{T}, dims)
mpzeros(::Type{Trop{T}}, dims::Tuple{}) where T <: MM = spzeros(Trop{T}, dims)

mpones(::Type{Trop{T}}, dims::DimOrInd...) where T <: MM = ones(Trop{T}, dims...)
mpones(::Type{Trop{T}}, dims::Tuple{Vararg{DimOrInd}}) where T <: MM  = ones(Trop{T}, dims)
mpones(::Type{Trop{T}}, dims::NTuple{N, Union{Integer, Base.OneTo}}) where {T<:MM,N}  = ones(Trop{T}, dims)
mpones(::Type{Trop{T}}, dims::NTuple{N, Integer}) where {T<:MM,N}  = ones(Trop{T}, dims)
mpones(::Type{Trop{T}}, dims::Tuple{}) where T <: MM = ones(Trop{T}, dims)

# ==============================================================================
# Needed Julia functions

Base.abs(x::Trop{T}) where T <: MM = Trop{T}(abs(x.λ))
Base.abs2(x::Trop{T}) where T <: MM = Trop{T}(x.λ + x.λ)
Base.round(x::Trop{T}, n::RoundingMode) where T <: MM = Trop{T}(round(x.λ, n))
Base.big(x::Trop{T}) where T <: MM = Base.big(x.λ)

# ==============================================================================
# Max-Plus matrices operations

# Map function f to a sparse matrice
function mpsparse_map(f, M::SparseMatrix)
    SparseMatrix(M.m, M.n, M.colptr, M.rowval, map(f, M.nzval))
end

# Eigenvalues and eigenvectors
include("howard.jl")

# Trace
mptrace(A::ArrMP) = isempty(A) ? mp0 : sum(diag(A))
mptrace(S::SpaMP) = isempty(S) ? mp0 : sum(diag(S))

# Norm
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
