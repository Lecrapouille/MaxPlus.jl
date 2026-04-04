# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.10
# In the way of the ScicosLab Max-Plus toolbox.
#
# The documentation of functions for the REPL are hidden in docstrings/
# ==============================================================================

module MaxPlus

using LinearAlgebra, SparseArrays, Printf

export
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
    ε, mp0, mp1, mpe, mi0, mi1, mie, mptop, mitop,
    # Matrix
    mpI, miI, eye, speye, spzeros, full, dense, plustimes,
    inv, tr, norm, astarb, star, plus, ones, zeros,
    # Spectral
    howard, mpeigen, semihoward,
    # Display
    set_tropical_display, tropshow, LaTeX

export # Max-Plus Linear system
    MPSysLin, simul, explicit, implicit, mpsimul, mpexplicit, mpimplicit, mpfull, mpsparse,
    flowshop, flowshop_graph, flowshop_simu, mpshift

# Fake templates to make difference between Min-Plus and Max-Plus numbers.
# These structures are not exported.
struct Min end
struct Max end
const MinOrMax = Union{Min, Max}

# Internal helper: absorbing element value for each algebra and float type.
# Max-Plus absorbing = -Inf, Min-Plus absorbing = +Inf.
_absorbing(::Type{Max}, ::Type{T}) where T <: AbstractFloat = typemin(T)
_absorbing(::Type{Min}, ::Type{T}) where T <: AbstractFloat = typemax(T)

# Base class for Max-Plus and Min-Plus structures.
# T is AbstractFloat to support Float32 (GPU/cuASR) and Float64 (default).
# Integer types are intentionally excluded: there is no exact -Inf in integers,
# so the absorbing element would overflow under tropical multiplication (x.λ + y.λ).
struct Tropical{Sense <: MinOrMax, T <: AbstractFloat} <: Real
    λ::T

    # Constructor with conversion of Float scalar in classic (+,*) algebra into
    # a (max,+) or (min,+) number. The scalar can be a 'Not a Number' value.
    Tropical{Sense,T}(n::AbstractFloat) where {Sense <: MinOrMax, T <: AbstractFloat} =
        new(isnan(n) ? _absorbing(Sense, T) : T(n))

    Tropical{Sense,T}(n::Integer) where {Sense <: MinOrMax, T <: AbstractFloat} =
        new(T(n))

    # Constructor with conversion from Bool into a (max,+) or (min,+). Bool is not
    # in tropical algebra but a hack to support the Julia generic identity matrix
    # operator I, replacing the Matlab eye() function, and that uses internaly Bool.
    Tropical{Sense,T}(b::Bool) where {Sense <: MinOrMax, T <: AbstractFloat} =
        new(b ? zero(T) : _absorbing(Sense, T))

    # Constructor with "tropicality" (max,+) from/to (min,+) conversion
    Tropical{Sense,T}(n::Tropical{<:MinOrMax,T}) where {Sense <: MinOrMax, T <: AbstractFloat} =
        new(n.λ)
end

# Outer constructor: captures T from the literal so MP(3.0f0) → Tropical{Sense,Float32},
# MP(3.0) → Tropical{Sense,Float64}. This makes MP and MI generic constructors rather
# than fixed Float64 aliases.
Tropical{Sense}(x::T) where {Sense <: MinOrMax, T <: AbstractFloat} = Tropical{Sense,T}(x)
# Bool before Integer: Bool <: Integer but must use tropical identity/absorbing (I matrix), not λ=1
Tropical{Sense}(b::Bool) where {Sense <: MinOrMax} = Tropical{Sense,Float64}(b)
Tropical{Sense}(x::Integer) where {Sense <: MinOrMax} = Tropical{Sense,Float64}(Float64(x))
Tropical{Sense}(n::Tropical{Sense,T}) where {Sense <: MinOrMax, T <: AbstractFloat} = n
# Copy λ between (max,+) and (min,+) (same float; algebra context changes)
Tropical{Max}(x::Tropical{Min,T}) where {T <: AbstractFloat} = Tropical{Max,T}(x.λ)
Tropical{Min}(x::Tropical{Max,T}) where {T <: AbstractFloat} = Tropical{Min,T}(x.λ)

# Max-Plus scalar type — generic constructor, T inferred from argument
const MP = Tropical{Max}
const MPAbstractArray{N} = AbstractArray{<:Tropical{Max},N}
const MPArray{N} = Array{<:Tropical{Max},N}
const MPMatrix = Matrix{<:Tropical{Max}}
const MPVector = Vector{<:Tropical{Max}}
const MPAbstractSparseArray{Ti,N} = AbstractSparseArray{<:Tropical{Max},Ti,N}
const MPSparseMatrix = SparseMatrixCSC{<:Tropical{Max}}
const MPSparseVector = SparseVector{<:Tropical{Max}}
const MPAbstractVecOrMat = AbstractVecOrMat{<:Tropical{Max}}

# Min-Plus structure (long name and short alias)
const MI = Tropical{Min}
const MIAbstractArray{N} = AbstractArray{<:Tropical{Min},N}
const MIArray{N} = Array{<:Tropical{Min},N}
const MIMatrix = Matrix{<:Tropical{Min}}
const MIVector = Vector{<:Tropical{Min}}
const MIAbstractSparseArray{Ti,N} = AbstractSparseArray{<:Tropical{Min},Ti,N}
const MISparseMatrix = SparseMatrixCSC{<:Tropical{Min}}
const MISparseVector = SparseVector{<:Tropical{Min}}
const MIAbstractVecOrMat = AbstractVecOrMat{<:Tropical{Min}}

# Legacy ScicosLab-style: Tropical(MP, x) / Tropical(MI, x) — one method (typeof(MP) === typeof(MI) as UnionAll)
function Tropical(which::UnionAll, x::Number)
    if which === Tropical{Max}
        return Tropical{Max}(x)
    elseif which === Tropical{Min}
        return Tropical{Min}(x)
    else
        throw(ArgumentError("first argument must be MP or MI (i.e. Tropical{Max} or Tropical{Min})"))
    end
end

# Definition of neutral and absorbing elements for (max,+) and (min,+).
# Fallback to Float64 when T is not specified.
Base.zero(::Type{Tropical{Max,T}}) where T <: AbstractFloat = Tropical{Max,T}(_absorbing(Max, T))
Base.zero(::Type{Tropical{Min,T}}) where T <: AbstractFloat = Tropical{Min,T}(_absorbing(Min, T))
Base.zero(::Type{Tropical{Max}}) = Tropical{Max,Float64}(_absorbing(Max, Float64))
Base.zero(::Type{Tropical{Min}}) = Tropical{Min,Float64}(_absorbing(Min, Float64))
Base.zero(x::Tropical) = zero(typeof(x))
Base.one(x::Tropical) = one(typeof(x))
Base.one(::Type{Tropical{Sense,T}}) where {Sense <: MinOrMax, T <: AbstractFloat} = Tropical{Sense,T}(zero(T))
Base.one(::Type{Tropical{Sense}}) where {Sense <: MinOrMax} = Tropical{Sense,Float64}(0.0)

# TropicalZero: universal absorbing element that adapts to both MP and MI.
# Follows the same pattern as LinearAlgebra.I (UniformScaling): a typed singleton
# that resolves to the correct zero once the algebraic context is known.
struct TropicalZero end
const ε = TropicalZero()

Base.:(+)(a::Tropical{S,T}, ::TropicalZero) where {S,T} = a
Base.:(+)(::TropicalZero, a::Tropical{S,T}) where {S,T} = a
Base.:(+)(::TropicalZero, ::TropicalZero)               = ε
Base.:(*)(::Tropical{S,T}, ::TropicalZero) where {S,T}  = zero(Tropical{S,T})
Base.:(*)(::TropicalZero, ::Tropical{S,T}) where {S,T}  = zero(Tropical{S,T})
Base.:(*)(::TropicalZero, ::TropicalZero)               = ε
Base.:(==)(::TropicalZero, ::TropicalZero)              = true
Base.:(==)(a::Tropical{S,T}, ::TropicalZero) where {S,T} = (a == zero(Tropical{S,T}))
Base.:(==)(::TropicalZero, a::Tropical{S,T}) where {S,T} = (a == zero(Tropical{S,T}))
Base.isless(::TropicalZero, ::TropicalZero) = false
Base.isless(::TropicalZero, ::Tropical{Max,T}) where T = true   # -Inf < tout
Base.isless(::Tropical{Max,T}, ::TropicalZero) where T = false
Base.isless(::TropicalZero, ::Tropical{Min,T}) where T = false  # +Inf > tout
Base.isless(::Tropical{Min,T}, ::TropicalZero) where T = true

# Shorter names for neutral and absorbing elements (all Float64 defaults)
# mp0 = -Inf (absorbing element for max), mp1 = 0 (neutral element for +)
# mi0 = +Inf (absorbing element for min), mi1 = 0 (neutral element for +)
const mp0 = zero(Tropical{Max,Float64})
const mp1 = one(Tropical{Max,Float64})
const mpe = mp1
const mptop = Tropical{Max,Float64}(Inf)
const mi0 = zero(Tropical{Min,Float64})
const mi1 = one(Tropical{Min,Float64})
const mie = mi1
const mitop = Tropical{Min,Float64}(-Inf)

# Promotion: any classic Julia Int or Float values are implicitly converted into
# a (max,+) or (min,+) number since we are suppose to work either in (max,+) or
# in (min,+) algebra and no more classic (+,*) algebra. Therefore for simplying
# the writting/reading of the code source, any Julia numbers are considered as
# tropical numbers to the current algebra. The context of the algebra shall be
# known. Therefore, this toolbox forbids mixture of operations on (max,+) and
# (min,+) numbers by throwing an error. This will force the developper commuting
# of algebra. Some calculus are easier made depending on the tropical algebra and
# the context may switch.
Base.promote_rule(::Type{Tropical{S,T1}}, ::Type{Tropical{S,T2}}) where {S <: MinOrMax, T1 <: AbstractFloat, T2 <: AbstractFloat} =
    Tropical{S, promote_type(T1,T2)}
Base.promote_rule(::Type{Tropical{Max,T}}, ::Type{Tropical{Min,T}}) where {T <: AbstractFloat} =
    error("Cannot promote (max,+) to (min,+)")
Base.promote_rule(::Type{Tropical{Min,T}}, ::Type{Tropical{Max,T}}) where {T <: AbstractFloat} =
    error("Cannot promote (min,+) to (max,+)")
# Classic reals and same-algebra tropicals; `Tropical <: Real` requires disambiguation from mixed MP/MI rules
function Base.promote_rule(::Type{Tropical{S,T}}, ::Type{R}) where {S <: MinOrMax, T <: AbstractFloat, R <: Real}
    if R <: Tropical{Max} && S === Min
        error("Cannot promote (max,+) to (min,+)")
    elseif R <: Tropical{Min} && S === Max
        error("Cannot promote (min,+) to (max,+)")
    end
    Tropical{S,T}
end
Base.convert(::Type{Tropical{S,T}}, x::Real) where {S <: MinOrMax, T <: AbstractFloat} =
    Tropical{S,T}(T(x))

# Convert a (max,+) or (min,+) scalar into classic (+,*) algebra.
Base.float(x::Tropical) = x.λ
Base.Float64(x::Tropical) = Float64(x.λ)
plustimes(x::Tropical) = x.λ

# Convert a (max,+) or (min,+) vector or matrix, dense or sparse into its
# equivalent in the classic (+,*).
plustimes(A::AbstractArray{<:Tropical}) = map(x -> x.λ, A)

# Define the (max,+) and (min,+) alegbra operators. Since any Julia number shall
# be considered as number of the current tropical, we have to force the implicit
# conversion.
Base.:(+)(x::Tropical{Max,T}, y::Tropical{Max,T}) where T = Tropical{Max,T}(max(x.λ, y.λ))
Base.:(+)(x::Tropical{Min,T}, y::Tropical{Min,T}) where T = Tropical{Min,T}(min(x.λ, y.λ))
Base.:(*)(x::Tropical{S,T},   y::Tropical{S,T})   where {S <: MinOrMax, T} = Tropical{S,T}(x.λ + y.λ)
Base.:(-)(x::Tropical{Max,T}, y::Tropical{Max,T}) where T = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::Tropical{Min,T}, y::Tropical{Min,T}) where T = error("Minus operator does not exist in (min,+) algebra")

# Comparaison operations. Since any Julia number shall be considered as number
# of the current tropical, we have to force the implicit conversion.
Base.:(==)(x::Tropical{S,T}, y::Tropical{S,T}) where {S,T} = (x.λ == y.λ)
Base.:(<=)(x::Tropical, y::Tropical) = (x.λ <= y.λ)
Base.:(<)(x::Tropical,  y::Tropical) = (x.λ <  y.λ)
Base.isless(x::Tropical, y::Tropical) = x.λ < y.λ
Base.:(==)(x::Tropical{S,T}, y::Real)  where {S,T} = x.λ == T(y)
Base.:(==)(x::Real, y::Tropical{S,T})  where {S,T} = T(x) == y.λ
Base.isless(x::Tropical, y::Real)  = x.λ < y
Base.isless(x::Real, y::Tropical)  = x < y.λ
Base.:(<=)(x::Tropical, y::Real)   = x.λ <= y
Base.:(<=)(x::Real, y::Tropical)   = x <= y.λ
Base.:(<)(x::Tropical,  y::Real)   = x.λ < y
Base.:(<)(x::Real,  y::Tropical)   = x < y.λ

# Power, division, residuation
Base.literal_pow(::typeof(^), x::Tropical{S,T}, ::Val{0}) where {S,T} = one(Tropical{S,T})
Base.:(^)(x::Tropical{S,T}, y::Int)     where {S,T} = Tropical{S,T}(x.λ * y)
Base.:(^)(x::Tropical{S,T}, y::Float64) where {S,T} = Tropical{S,T}(x.λ * y)
Base.:(-)(x::Tropical{S,T}) where {S,T} =
    (x == zero(Tropical{S,T})) ? zero(Tropical{S,T}) : Tropical{S,T}(-x.λ)
Base.sign(x::Tropical) = Base.sign(x.λ)
Base.inv(x::Tropical{S,T}) where {S,T} = Tropical{S,T}(-x.λ)
Base.inv(A::Matrix{<:Tropical{S}})            where S = Matrix(transpose(-A))
Base.inv(S::SparseMatrixCSC{<:Tropical{Se}})  where Se = SparseMatrixCSC(transpose(-S))
Base.:(\)(x::Tropical{S,T}, y::Tropical{S,T}) where {S,T} = inv(inv(y) * x)
Base.:(\)(A::AbstractMatrix{<:Tropical{S}}, B::AbstractMatrix{<:Tropical{S}}) where S = inv(inv(B) * A)
Base.:(\)(A::AbstractMatrix{<:Tropical{S}}, B::Tropical{S,T})                 where {S,T} = inv(inv(B) * A)
Base.:(\)(A::Tropical{S,T}, B::AbstractMatrix{<:Tropical{S}})                 where {S,T} = inv(inv(B) * A)
Base.:(/)(x::Tropical{S,T}, y::Tropical{S,T}) where {S,T} = inv(y * inv(x))
Base.:(/)(A::AbstractMatrix{<:Tropical{S}}, B::AbstractMatrix{<:Tropical{S}}) where S = inv(B * inv(A))
Base.:(/)(A::AbstractMatrix{<:Tropical{S}}, B::Tropical{S,T})                 where {S,T} = inv(B * inv(A))
Base.:(/)(A::Tropical{S,T}, B::AbstractMatrix{<:Tropical{S}})                 where {S,T} = inv(B * inv(A))
Base.:min(x::Tropical{Max,T}, y::Tropical{Max,T}) where T = (x * y) / (x + y)
Base.:min(x::AbstractMatrix{<:Tropical{Max}}, y::AbstractMatrix{<:Tropical{Max}}) = (x * y) / (x + y)

# Compatibility with Julia builtin operators.
Base.isnan(x::Tropical)  = isnan(x.λ)
Base.isinf(x::Tropical)  = isinf(x.λ)
Base.abs(x::Tropical{S,T})   where {S,T} = Tropical{S,T}(abs(x.λ))
Base.abs2(x::Tropical{S,T})  where {S,T} = Tropical{S,T}(x.λ + x.λ)
Base.round(x::Tropical{S,T}, n::RoundingMode=RoundNearest) where {S,T} = Tropical{S,T}(round(x.λ, n))
Base.trunc(x::Tropical{S,T}) where {S,T} = Tropical{S,T}(round(x.λ, RoundToZero))
Base.floor(x::Tropical{S,T}) where {S,T} = Tropical{S,T}(round(x.λ, RoundDown))
Base.ceil(x::Tropical{S,T})  where {S,T} = Tropical{S,T}(round(x.λ, RoundUp))
Base.big(x::Tropical)  = Base.big(x.λ)

# Since Julia has remove the Matlab eye() function, making identity matrix, now
# replaced by the Julia I operator, we reintroduce it.
eye(::Type{Tropical{S,T}}, n::Int64)         where {S,T} = Matrix{Tropical{S,T}}(I, n, n)
eye(::Type{Tropical{S,T}}, m::Int64, n::Int64) where {S,T} = Matrix{Tropical{S,T}}(I, m, n)
eye(A::AbstractMatrix{Tropical{S,T}})          where {S,T} = Matrix{Tropical{S,T}}(mpI, size(A))
# MP / MI are the same UnionAll `typeof`; dispatch with === on Tropical{Max} / Tropical{Min}
function eye(Tmod::UnionAll, n::Int64)
    Tmod === Tropical{Max} && return eye(Tropical{Max,Float64}, n)
    Tmod === Tropical{Min} && return eye(Tropical{Min,Float64}, n)
    throw(ArgumentError("eye: first argument must be MP or MI (Tropical{Max} / Tropical{Min})"))
end
function eye(Tmod::UnionAll, m::Int64, n::Int64)
    Tmod === Tropical{Max} && return eye(Tropical{Max,Float64}, m, n)
    Tmod === Tropical{Min} && return eye(Tropical{Min,Float64}, m, n)
    throw(ArgumentError("eye: first argument must be MP or MI (Tropical{Max} / Tropical{Min})"))
end
speye(::Type{Tropical{S,T}}, n::Int64)         where {S,T} = SparseMatrixCSC{Tropical{S,T}}(I, n, n)
speye(::Type{Tropical{S,T}}, m::Int64, n::Int64) where {S,T} = SparseMatrixCSC{Tropical{S,T}}(I, m, n)
speye(A::AbstractMatrix{Tropical{S,T}})          where {S,T} = SparseMatrixCSC{Tropical{S,T}}(mpI, size(A))
function speye(Tmod::UnionAll, n::Int64)
    Tmod === Tropical{Max} && return speye(Tropical{Max,Float64}, n)
    Tmod === Tropical{Min} && return speye(Tropical{Min,Float64}, n)
    throw(ArgumentError("speye: first argument must be MP or MI"))
end
function speye(Tmod::UnionAll, m::Int64, n::Int64)
    Tmod === Tropical{Max} && return speye(Tropical{Max,Float64}, m, n)
    Tmod === Tropical{Min} && return speye(Tropical{Min,Float64}, m, n)
    throw(ArgumentError("speye: first argument must be MP or MI"))
end

# Julia builtin zero, zeros, spzero, and their equivalent functions for ones are
# defacto good for building scalar and matrices. We complete some functions.
Base.zeros(A::AbstractArray{Tropical{S,T},N}) where {S,T,N} = zeros(Tropical{S,T}, size(A))
SparseArrays.spzeros(A::Array{Tropical{S,T},N}) where {S,T,N} = spzeros(Tropical{S,T}, size(A))
Base.ones(A::AbstractArray{Tropical{S,T},N})  where {S,T,N} = ones(Tropical{S,T}, size(A))

function Base.zeros(Tmod::UnionAll, dims::Integer...)
    Tmod === Tropical{Max} && return zeros(Tropical{Max,Float64}, dims...)
    Tmod === Tropical{Min} && return zeros(Tropical{Min,Float64}, dims...)
    throw(ArgumentError("zeros: first argument must be MP or MI"))
end
function Base.ones(Tmod::UnionAll, dims::Integer...)
    Tmod === Tropical{Max} && return ones(Tropical{Max,Float64}, dims...)
    Tmod === Tropical{Min} && return ones(Tropical{Min,Float64}, dims...)
    throw(ArgumentError("ones: first argument must be MP or MI"))
end
function SparseArrays.spzeros(Tmod::UnionAll, dims::Integer...)
    Tmod === Tropical{Max} && return spzeros(Tropical{Max,Float64}, dims...)
    Tmod === Tropical{Min} && return spzeros(Tropical{Min,Float64}, dims...)
    throw(ArgumentError("spzeros: first argument must be MP or MI"))
end

# Equivalent to I but for (max,+) and (min,+)
const global mpI = UniformScaling(one(Tropical{Max,Float64}))
const global miI = UniformScaling(one(Tropical{Min,Float64}))

# Constructor converting a dense/sparse vector or matrix from classic (+,*) algebra
# into a (max,+) or (min,+) dense/sparse vector or matrix.
# Accepts Integer/Real entries (not only AbstractFloat) so e.g. MP([1, 2, 3]) and
# MP(sparse([1 2; 0 4])) work.
Tropical{Sense}(A::AbstractArray{<:Real}) where {Sense <: MinOrMax} =
    map(x -> Tropical{Sense}(x), A)
sparse_map(f, S::SparseMatrixCSC) = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(f, S.nzval))
sparse_map(f, V::SparseVector)    = SparseVector(V.n, V.nzind, map(f, V.nzval))
Tropical{Sense}(S::SparseMatrixCSC{<:Real}) where {Sense <: MinOrMax} =
    sparse_map(x -> Tropical{Sense}(x), S)
Tropical{Sense}(S::SparseVector{<:Real}) where {Sense <: MinOrMax} =
    sparse_map(x -> Tropical{Sense}(x), S)
# Same than sparse(I,J,Tv) but with dropping zero() elements
Tropical{Sense}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{<:Real}) where {Sense <: MinOrMax, Ti} =
    dropzeros(sparse(I, J, map(x -> Tropical{Sense}(x), V)))
Tropical{Sense}(I::AbstractVector{Ti}, V::AbstractVector{<:Real}) where {Sense <: MinOrMax, Ti} =
    dropzeros(sparsevec(I, map(x -> Tropical{Sense}(x), V)))

full(S::AbstractArray{<:Tropical})      = S
full(S::SparseMatrixCSC{<:Tropical})    = Matrix(S)
full(S::SparseVector{<:Tropical})       = Vector(S)
dense(S::AbstractArray{<:Tropical})     = S
dense(S::SparseMatrixCSC{<:Tropical})   = Matrix(S)
dense(S::SparseVector{<:Tropical})      = Vector(S)

# TODO
# Constructor converting a range (i.e [1.. 3]) from classic (+,*) algebra into
# (max,+) or (min,+) range.
Tropical{Sense}(x::UnitRange) where {Sense <: MinOrMax} = [Tropical{Sense}(i) for i in x]
Tropical{Sense}(x::StepRangeLen) where {Sense <: MinOrMax} = [Tropical{Sense}(Float64(i)) for i in x]

# Redefine LinearAlgebra trace and norm.
LinearAlgebra.tr(A::AbstractMatrix{<:Tropical{Max}})          = sum(diag(A))
LinearAlgebra.norm(S::AbstractMatrix{Tropical{Max,T}}) where T = Tropical{Max,T}(S[argmax(S)].λ - S[argmin(S)].λ)
LinearAlgebra.tr(A::SparseMatrixCSC{<:Tropical{Max}})         = sum(diag(A))
LinearAlgebra.norm(S::SparseMatrixCSC{Tropical{Max,T}}) where T = Tropical{Max,T}(S[argmax(S)].λ - S[argmin(S)].λ)

# Do the A*, A+, Ax+b operations
squared_size(A::AbstractVector) = error("Matrix shall be squared")
squared_size(A::AbstractMatrix) = (n = size(A, 1); (n != size(A, 2)) && error("Matrix shall be squared") || n)
@inline diag_map!(f, A::AbstractMatrix) = for i = 1:size(A,1) A[i,i] = f(A[i,i]) end
star(A::AbstractVector{<:Tropical{Max}}) = error("Matrix shall be squared")
function star(A::AbstractMatrix{Tropical{Max,T}}) where T
    s = Int(round(log2(squared_size(A)))) # FIXME float ok ?
    M = A
    f = x -> x < 0.0 ? one(Tropical{Max,T}) : x > 0.0 ? mptop : 0.0
    diag_map!(f, M)
    for i = 0:s
        diag_map!(f, M)
        M *= M
    end
    M
end
star(x::Tropical{Max,T}) where T = star(hcat(x))[1,1]
plus(A::AbstractVector{<:Tropical{Max}}) = error("Matrix shall be squared")
plus(A::AbstractMatrix{<:Tropical{Max}}) = A * star(A)
plus(x::Tropical{Max,T}) where T = x * star(x)
astarb(A::AbstractMatrix{<:Tropical{Max}}, b::AbstractVector{<:Tropical{Max}}) = star(A) * b

include("fallbacks.jl")

include("howard.jl")
include("syslin.jl")
include("flowshop.jl")
include("show.jl")
include("docstrings/mp.jl")
include("docstrings/mi.jl")
include("docstrings/syslin.jl")

end # MaxPlus module