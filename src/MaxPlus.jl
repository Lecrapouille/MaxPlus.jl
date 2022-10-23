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
    # Struct
    Tropical, MP, MI,
    MPAbstractArray, MPArray, MPMatrix, MPVector,
    MPAbstractSparseArray, MPSparseMatrix, MPSparseVector,
    MPAbstractVecOrMat,
    MIAbstractArray, MIArray, MIMatrix, MIVector,
    MIAbstractSparseArray, MISparseMatrix, MISparseVector,
    MIAbstractVecOrMat,
    # Algebra
    zero, one,
    # Constants
    mp0, ε,  mp1, mpe, mi0, mi1, mie, mptop, mitop,
    # Matrix
    mpI, miI, eye, speye, full, dense, array, plustimes, inv,
    tr, norm, astarb, star, plus, howard,
    # Display
    set_tropical_display, tropshow, LaTeX

export # Max-Plus Linear system
    MPSysLin, mpsimul, mpexplicit

export # Max-Plus flowhop
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

# Fake templates to make difference between Min-Plus and Max-Plus numbers.
# These structures are not exported.
struct Min end
struct Max end
const MinOrMax = Union{Min, Max}

# Definition of neutral and absorbing elements for (max,+) and (min,+).
# Needed for Tropical{T} inner constructors.
Base.zero(::Type{Max}) = typemin(Float64)
Base.zero(::Type{Min}) = typemax(Float64)
Base.one(::Type{Max}) = zero(Float64)
Base.one(::Type{Min}) = zero(Float64)

# Base class for Max-Plus and Min-Plus structures
struct Tropical{T <: MinOrMax} <: Number
    λ::Float64

    Tropical{T}(n::Int64) where {T<:MinOrMax} = new(Float64(n))

    # Constructor with conversion of Float scalar in classic (+,*) algebra into
    # a (max,+) or (min,+) number. The scalar can be a 'Not a Number' value.
    Tropical{T}(n::Float64) where {T<:MinOrMax} = new(isnan(n) ? zero(T) : n)

    # Constructor with conversion from Bool into a (max,+) or (min,+). Bool is not
    # in tropical algebra but a hack to support the Julia generic identity matrix
    # operator I, replacing the Matlab eye() function, and that uses internaly Bool.
    Tropical{T}(n::Bool) where {T<:MinOrMax} = new(n ? one(T) : zero(T))

    # Constructor with "tropicality" (max,+) from/to (min,+) conversion
    Tropical{T}(n::Tropical) where {T<:MinOrMax} = new(n.λ)
end

# Max-Plus structure
const MP = Tropical{Max}
const MPAbstractArray{N} = AbstractArray{MP,N}
const MPArray{N} = Array{MP,N}
const MPMatrix = Matrix{MP}
const MPVector = Vector{MP}
const MPAbstractSparseArray{Ti,N} = AbstractSparseArray{MP,Ti,N}
const MPSparseMatrix = SparseMatrixCSC{MP}
const MPSparseVector = SparseVector{MP}
const MPAbstractVecOrMat = AbstractVecOrMat{MP}

# Min-Plus structure
const MI = Tropical{Min}
const MIAbstractArray{N} = AbstractArray{MI,N}
const MIArray{N} = Array{MI,N}
const MIMatrix = Matrix{MI}
const MIVector = Vector{MI}
const MIAbstractSparseArray{Ti,N} = AbstractSparseArray{MI,Ti,N}
const MISparseMatrix = SparseMatrixCSC{MI}
const MISparseVector = SparseVector{MI}
const MIAbstractVecOrMat = AbstractVecOrMat{MI}

# Definition of neutral and absorbing elements for (max,+) and (min,+).
Base.zero(::Type{MP}) = MP(typemin(Float64))
Base.zero(::Type{MI}) = MI(typemax(Float64))
Base.zero(x::Tropical) = zero(typeof(x))
Base.one(x::Tropical) = one(typeof(x))
Base.one(::Type{Tropical{T}}) where {T<:MinOrMax} = Tropical{T}(zero(Float64))

# Shorter names for neutral and absorbing elements
const global mp0 = zero(MP)
const global ε = zero(MP) # FIXME ε is also == zero(MI) how to deal ?
const global mp1 = one(MP)
const global mpe = one(MP)
const global mptop = MP(Inf)
const global mi0 = zero(MI)
const global mi1 = one(MI)
const global mie = one(MI)
const global mitop = MI(-Inf)

Tropical(::Type{T}, n::Float64) where T = isnan(n) ? zero(T) : T(n)
# Constructor converting a dense/sparse vector or matrix from classic (+,*) algebra
# into a (max,+) or (min,+) dense/sparse vector or matrix.
Tropical{T}(A::AbstractArray) where {T<:MinOrMax} = map(Tropical{T}, A)
sparse_map(f, S::SparseMatrixCSC) = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(f, S.nzval))
sparse_map(f, V::SparseVector) = SparseVector(V.n, V.nzind, map(f, V.nzval))
Tropical{T}(S::SparseMatrixCSC) where {T<:MinOrMax} = sparse_map(Tropical{T}, S)
Tropical{T}(S::SparseVector) where {T<:MinOrMax} = sparse_map(Tropical{T}, S)
# Same than sparse(I,J,Tv) but with dropping zero() elements
Tropical{T}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MinOrMax,Tv,Ti} = dropzeros(sparse(I, J, Tropical{T}(V)))
Tropical{T}(I::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MinOrMax,Tv,Ti} = dropzeros(sparsevec(I, Tropical{T}(V)))
full(S::AbstractArray{Tropical{T}}) where {T<:MinOrMax} = S
full(S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax} = Matrix(S)
full(S::SparseVector{Tropical{T}}) where {T<:MinOrMax} = Vector(S)
dense(S::AbstractArray{Tropical{T}}) where {T<:MinOrMax} = S
dense(S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax} = Matrix(S)
dense(S::SparseVector{Tropical{T}}) where {T<:MinOrMax} = Vector(S)

# TODO
# Constructor converting a range (i.e [1.. 3]) from classic (+,*) algebra into
# (max,+) or (min,+) range.
Tropical{T}(x::UnitRange) where {T<:MinOrMax} = Vector{Tropical{T}}(x)
Tropical{T}(x::StepRangeLen) where {T<:MinOrMax} = Vector{Tropical{T}}(x)

# Promotion: any classic Julia Int or Float values are implicitly converted into
# a (max,+) or (min,+) number since we are suppose to work either in (max,+) or
# in (min,+) algebra and no more classic (+,*) algebra. Therefore for simplying
# the writting/reading of the code source, any Julia numbers are considered as
# tropical numbers to the current algebra. The context of the algebra shall be
# known. Therefore, this toolbox forbids mixture of operations on (max,+) and
# (min,+) numbers by throwing an error. This will force the developper commuting
# of algebra. Some calculus are easier made depending on the tropical algebra and
# the context may switch.
Base.promote_rule(::Type{Tropical{T2}}, ::Type{T1}) where {T1,T2<:MinOrMax} = Tropical{T2}
Base.promote_rule(::Type{MP}, ::Type{MI}) = error("Cannot promote (max,+) to (min,+)")
Base.promote_rule(::Type{MI}, ::Type{MP}) = error("Cannot promote (min,+) to (max,+)")
Base.convert(::Tropical{T}, x::Number) where {T<:MinOrMax} = Tropical{T}(x)

# Convert a (max,+) or (min,+) scalar into classic (+,*) algebra.
Base.float(x::Tropical) = x.λ
Float64(x::Tropical) = x.λ
plustimes(x::Tropical) = x.λ

# Convert a (max,+) or (min,+) vector or matrix, dense or sparse into its
# equivalent in the classic (+,*).
plustimes(A::AbstractArray{Tropical{T}}) where {T<:MinOrMax} = map(x -> x.λ, A)

# Define the (max,+) and (min,+) alegbra operators. Since any Julia number shall
# be considered as number of the current tropical, we have to force the implicit
# conversion.
Base.:(+)(x::MP, y::MP) = MP(max(x.λ, y.λ))
Base.:(+)(x::MP, y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP) = MP(max(x, y.λ))
Base.:(+)(x::MI, y::MI) = MI(min(x.λ, y.λ))
Base.:(+)(x::MI, y::Real) = MI(min(x.λ, y))
Base.:(+)(x::Real, y::MI) = MI(min(x, y.λ))
Base.:(*)(x::MP, y::MP) = Tropical(MP, x.λ + y.λ)
Base.:(*)(x::MP, y::Real) = Tropical(MP, x.λ + y)
Base.:(*)(x::Real, y::MP) = Tropical(MP, x + y.λ)
Base.:(*)(x::MI, y::MI) = Tropical(MI, x.λ + y.λ)
Base.:(*)(x::MI, y::Real) = Tropical(MI, x.λ + y)
Base.:(*)(x::Real, y::MI) = Tropical(MI, x + y.λ)
Base.:(-)(x::MP, y::MP) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::MP, y::Real) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::Real, y::MP) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::MI, y::MI) = error("Minus operator does not exist in (min,+) algebra")
Base.:(-)(x::MI, y::Real) = error("Minus operator does not exist in (min,+) algebra")
Base.:(-)(x::Real, y::MI) = error("Minus operator does not exist in (min,+) algebra")

# Comparaison operations. Since any Julia number shall be considered as number
# of the current tropical, we have to force the implicit conversion.
Base.:(==)(x::Tropical{T}, y::Tropical{T}) where {T<:MinOrMax} = (x.λ == y.λ)
Base.:(==)(x::Tropical, y::Real) = (x.λ == y)
Base.:(==)(x::Real, y::Tropical) = (x == y.λ)
Base.:(<=)(x::Tropical, y::Tropical) = (x.λ <= y.λ)
Base.:(<=)(x::Tropical, y::Real) = (x.λ <= y)
Base.:(<=)(x::Real, y::Tropical) = (x <= y.λ)
Base.:(<)(x::Tropical, y::Tropical) = (x.λ < y.λ)
Base.:(<)(x::Tropical, y::Real) = (x.λ < y)
Base.:(<)(x::Real, y::Tropical) = (x < y.λ)
Base.isless(x::Tropical, y::Tropical) = x.λ < y.λ
Base.isless(x::Tropical, y::Real) = x.λ < y
Base.isless(x::Real, y::Tropical) = x < y.λ

# Power, division, residuation
Base.literal_pow(::typeof(^), x::Tropical{T}, ::Val{0}) where {T<:MinOrMax} = one(Tropical{T})
Base.:(^)(x::Tropical{T}, y::Int) where {T<:MinOrMax} = Tropical{T}(x.λ * y)
Base.:(^)(x::Tropical{T}, y::Float64) where {T<:MinOrMax} = Tropical{T}(x.λ * y)
Base.:(-)(x::Tropical{T}) where {T<:MinOrMax} = (x == zero(Tropical{T})) ? zero(Tropical{T}) : Tropical{T}(-x.λ)
Base.sign(x::Tropical) = Base.sign(x.λ)
Base.inv(x::Tropical{T}) where {T<:MinOrMax} = Tropical{T}(-x.λ)
Base.inv(A::Matrix{MP}) = Matrix(transpose(-A))
Base.inv(S::SparseMatrixCSC{MP}) = SparseMatrixCSC(transpose(-S))
Base.:(\)(x::Tropical{T}, y::Tropical{T}) where {T<:MinOrMax} = inv(inv(y) * x)
Base.:(\)(A::AbstractMatrix{Tropical{T}}, B::AbstractMatrix{Tropical{T}}) where {T<:MinOrMax} = inv(inv(B) * A)
Base.:(\)(A::AbstractMatrix{Tropical{T}}, B::Tropical{T}) where {T<:MinOrMax} = inv(inv(B) * A)
Base.:(\)(A::Tropical{T}, B::AbstractMatrix{Tropical{T}}) where {T<:MinOrMax} = inv(inv(B) * A)
Base.:(/)(x::Tropical{T}, y::Tropical{T}) where {T<:MinOrMax} = inv(y * inv(x))
Base.:(/)(A::AbstractMatrix{Tropical{T}}, B::AbstractMatrix{Tropical{T}}) where {T<:MinOrMax} = inv(B * inv(A))
Base.:(/)(A::AbstractMatrix{Tropical{T}}, B::Tropical{T}) where {T<:MinOrMax} = inv(B * inv(A))
Base.:(/)(A::Tropical{T}, B::AbstractMatrix{Tropical{T}}) where {T<:MinOrMax} = inv(B * inv(A))
Base.:min(x::MP, y::MP) = (x * y) / (x + y)
Base.:min(x::AbstractMatrix{MP}, y::AbstractMatrix{MP}) = (x * y) / (x + y)
Base.:min(x::MP, y::Real) = min(x, MP(y))
Base.:min(x::Real, y::MP) = min(MP(x), y)

# Compatibility with Julia builtin operators.
Base.isnan(x::Tropical) = isnan(x.λ)
Base.isinf(x::Tropical) = isinf(x.λ)
Base.abs(x::Tropical{T}) where {T<:MinOrMax} = Tropical{T}(abs(x.λ))
Base.abs2(x::Tropical{T}) where {T<:MinOrMax} = Tropical{T}(x.λ + x.λ)
Base.round(x::Tropical{T}, n::RoundingMode=RoundNearest) where {T<:MinOrMax} = Tropical{T}(round(x.λ, n))
Base.trunc(x::Tropical{T}) where {T<:MinOrMax} = Tropical{T}(round(x.λ, RoundToZero))
Base.floor(x::Tropical{T}) where {T<:MinOrMax} = Tropical{T}(round(x.λ, RoundDown))
Base.ceil(x::Tropical{T}) where {T<:MinOrMax} = Tropical{T}(round(x.λ, RoundUp))
Base.big(x::Tropical) = Base.big(x.λ)

# Since Julia has remove the Matlab eye() function, making identity matrix, now
# replaced by the Julia I operator, we reintroduce it.
eye(::Type{Tropical{T}}, n::Int64) where {T<:MinOrMax} = Matrix{Tropical{T}}(I, n, n)
eye(::Type{Tropical{T}}, m::Int64, n::Int64) where {T<:MinOrMax} = Matrix{Tropical{T}}(I, m, n)
eye(A::AbstractMatrix{Tropical{T}}) where {T<:MinOrMax,N} = Matrix{Tropical{T}}(mpI, size(A))
speye(::Type{Tropical{T}}, n::Int64) where {T<:MinOrMax} = SparseMatrixCSC{Tropical{T}}(I, n, n)
speye(::Type{Tropical{T}}, m::Int64, n::Int64) where {T<:MinOrMax} = SparseMatrixCSC{Tropical{T}}(I, m, n)
speye(A::AbstractMatrix{Tropical{T}}) where {T<:MinOrMax,N} = SparseMatrixCSC{Tropical{T}}(mpI, size(A))

# Julia builtin zero, zeros, spzero, and their equivalent functions for ones are
# defacto good for building scalar and matrices. We complete some functions.
Base.zeros(A::AbstractArray{Tropical{T},N}) where {T<:MinOrMax,N} = zeros(Tropical{T}, size(A))
SparseArrays.spzeros(A::Array{Tropical{T},N}) where {T<:MinOrMax,N} = spzeros(Tropical{T}, size(A))
Base.ones(A::AbstractArray{Tropical{T},N}) where {T<:MinOrMax,N} = ones(Tropical{T}, size(A))

# Equivalent to I but for (max,+) and (min,+)
const global mpI = UniformScaling(one(MP))
const global miI = UniformScaling(one(MI))

# Redefine LinearAlgebra trace and norm.
LinearAlgebra.tr(A::AbstractMatrix{MP}) = sum(diag(A))
LinearAlgebra.norm(S::AbstractMatrix{MP}) = MP(S[argmax(S)].λ - S[argmin(S)].λ)
LinearAlgebra.tr(A::SparseMatrixCSC{MP}) = sum(diag(A))
LinearAlgebra.norm(S::SparseMatrixCSC{MP}) = MP(S[argmax(S)].λ - S[argmin(S)].λ)

# Do the A*, A+, Ax+b operations
squared_size(A::AbstractVector) = error("Matrix shall be squared")
squared_size(A::AbstractMatrix) = (n = size(A, 1); (n != size(A, 2)) && error("Matrix shall be squared") || n)
@inline diag_map!(f, A::AbstractMatrix) = for i = 1:size(A,1) A[i,i] = f(A[i,i]) end
star(A::AbstractVector{MP}) = error("Matrix shall be squared")
function star(A::AbstractMatrix{MP})
    s = Int(round(log2(squared_size(A)))) # FIXME float ok ?
    M = A
    f = x -> x < 0.0 ? mp1 : x > 0.0 ? mptop : 0.0
    diag_map!(f, M)
    for i = 0:s
        diag_map!(f, M)
        M *= M
    end
    M
end
star(x::MP) = star(hcat(x))[1,1]
plus(A::AbstractVector{MP}) = error("Matrix shall be squared")
plus(A::AbstractMatrix{MP}) = A * star(A)
plus(x::MP) = x * star(x)
astarb(A::AbstractMatrix{MP}, b::AbstractVector{MP}) = star(A) * b

include("fallbacks.jl")
include("howard.jl")
include("syslin.jl")
include("show.jl")
include("docstrings.jl")

end # MaxPlus module
