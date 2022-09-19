# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# A portage of the ScicosLab Max-Plus toolbox http://www.scicoslab.org/
# License: public domain
#
# Note: the documentation of functions for the REPL are placed in docstrings.jl
# ==============================================================================
# This file defines (max,+) and (min,+) core functions.
# ==============================================================================

using
    LinearAlgebra

export
    mp0, mp1, mi0, mi1, mptop, mitop,
    eye, full, dense, array, plustimes, inv,
    tr, norm, mpastarb, star, plus

Tropical{T}(x::Bool) where {T<:MM} = x ? one(Tropical{T}) : zero(Tropical{T})
Tropical{T}(x::Tropical) where {T<:MM} = Tropical{T}(x.v)
Tropical(::Type{T}, n::Float64) where T = isnan(n) ? zero(T) : T(n)
Tropical{T}(A::Array) where {T<:MM} = map(Tropical{T}, A)
Tropical{T}(S::SparseMatrixCSC{T1,U}) where {T<:MM,T1,U} = convert(TropicalSparseMatrix{T,U}, S)
Tropical{T}(V::SparseVector{T1}) where {T<:MM,T1} = convert(SparseVector{Tropical{T}}, V)
Tropical{T}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MM,Tv,Ti<:Integer} = sparse(I, J, Tropical{T}(V)) # FIXME: to be removed
Tropical{T}(x::UnitRange) where {T<:MM} = Vector{Tropical{T}}(x)
Tropical{T}(x::StepRangeLen) where {T<:MM} = Vector{Tropical{T}}(x)
Base.promote_rule(::Type{Tropical{T2}}, ::Type{T1}) where {T1,T2} = Tropical{T2}
Base.promote_rule(::Type{MaxPlus}, ::Type{MinPlus}) = error("Cannot promote (max,+) to (min,+) number")
Base.promote_rule(::Type{MinPlus}, ::Type{MaxPlus}) = error("Cannot promote (min,+) to (max,+) number")
Base.convert(::Tropical{T}, x::Number) where {T<:MM} = Tropical{T}(x)
Base.float(x::Tropical) = x.v
Base.isnan(x::Tropical) = isnan(x.v)
Base.zero(::Type{MaxPlus}) = MaxPlus(typemin(Float64))
Base.zero(::Type{MinPlus}) = MinPlus(typemax(Float64))
Base.zero(x::Tropical) = zero(typeof(x))
Base.one(x::Tropical) = one(typeof(x))
Base.one(::Type{Tropical{T}}) where {T<:MM} = Tropical{T}(zero(Float64))
Base.:(+)(x::MaxPlus, y::MaxPlus) = MaxPlus(max(x.v, y.v))
Base.:(+)(x::MaxPlus, y::Real) = MaxPlus(max(x.v, y))
Base.:(+)(x::Real, y::MaxPlus) = MaxPlus(max(x, y.v))
Base.:(+)(x::MinPlus, y::MinPlus) = MinPlus(min(x.v, y.v))
Base.:(+)(x::MinPlus, y::Real) = MinPlus(min(x.v, y))
Base.:(+)(x::Real, y::MinPlus) = MinPlus(min(x, y.v))
Base.:(*)(x::MaxPlus, y::MaxPlus) = Tropical(MaxPlus, x.v + y.v)
Base.:(*)(x::MaxPlus, y::Real) = Tropical(MaxPlus, x.v + y)
Base.:(*)(x::Real, y::MaxPlus) = Tropical(MaxPlus, x + y.v)
Base.:(*)(x::MinPlus, y::MinPlus) = Tropical(MinPlus, x.v + y.v)
Base.:(*)(x::MinPlus, y::Real) = Tropical(MinPlus, x.v + y)
Base.:(*)(x::Real, y::MinPlus) = Tropical(MinPlus, x + y.v)
Base.:(-)(x::MaxPlus, y::MaxPlus) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::MaxPlus, y::Real) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::Real, y::MaxPlus) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::MinPlus, y::MinPlus) = error("Minus operator does not exist in (min,+) algebra")
Base.:(-)(x::MinPlus, y::Real) = error("Minus operator does not exist in (min,+) algebra")
Base.:(-)(x::Real, y::MinPlus) = error("Minus operator does not exist in (min,+) algebra")
Base.literal_pow(::typeof(^), x::Tropical, ::Val{0}) = one(x)
Base.:(^)(x::Tropical{T}, y::Int) where {T<:MM} = Tropical{T}(x.v * y)
Base.:(^)(x::Tropical{T}, y::Float64) where {T<:MM} = Tropical{T}(x.v * y)
Base.:(-)(x::Tropical{T}) where {T<:MM} = (x == zero(Tropical{T})) ? zero(Tropical{T}) : Tropical{T}(-x.v)
Base.sign(x::Tropical) = Base.sign(x.v)
Base.inv(x::Tropical{T}) where {T<:MM} = Tropical{T}(-x.v)
Base.inv(A::ArrMaxPlus) = Matrix(transpose(-A))
Base.inv(S::SpaMaxPlus) = SparseMatrixCSC(transpose(-S))
Base.:(\)(x::Tropical{T}, y::Tropical{T}) where {T<:MM} = inv(inv(y) * x)
Base.:(\)(A::TropicalArray{T}, B::TropicalArray{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(\)(A::TropicalArray{T}, B::Tropical{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(\)(A::Tropical{T}, B::TropicalArray{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(/)(x::Tropical{T}, y::Tropical{T}) where {T<:MM} = inv(y * inv(x))
Base.:(/)(A::TropicalArray{T}, B::TropicalArray{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::TropicalArray{T}, B::Tropical{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::Tropical{T}, B::TropicalArray{T}) where {T<:MM} = inv(B * inv(A))

#FIXME
Base.:(/)(A::TropicalSparseMatrix{T,U}, B::TropicalArray{T}) where {T<:MM,U} = inv(B * inv(A))

Base.:(==)(x::Tropical{T}, y::Tropical{T}) where {T<:MM} = (x.v == y.v)
Base.:(==)(x::Tropical, y::Real) = (x.v == y)
Base.:(==)(x::Real, y::Tropical) = (x == y.v)
Base.:(<=)(x::Tropical, y::Tropical) = (x.v <= y.v)
Base.:(<=)(x::Tropical, y::Real) = (x.v <= y)
Base.:(<=)(x::Real, y::Tropical) = (x <= y.v)
Base.:(<)(x::Tropical, y::Tropical) = (x.v < y.v)
Base.:(<)(x::Tropical, y::Real) = (x.v < y)
Base.:(<)(x::Real, y::Tropical) = (x < y.v)
Base.isless(x::Tropical, y::Tropical) = x.v < y.v
Base.isless(x::Tropical, y::Real) = x.v < y
Base.isless(x::Real, y::Tropical) = x < y.v
Base.:min(x::MaxPlus, y::MaxPlus) = (x * y) / (x + y)
Base.:min(x::MaxPlus, y::Real) = min(x, MaxPlus(y))
Base.:min(x::Real, y::MaxPlus) = min(MaxPlus(x), y)
#Base.:min(A::ArrMaxPlus, B::ArrMaxPlus) = map(Base.:min, A, B)
#Base.:min(A::SpaMaxPlus, B::SpaMaxPlus) = map(Base.:min, A, B)
Base.abs(x::Tropical{T}) where {T<:MM} = Tropical{T}(abs(x.v))
Base.abs2(x::Tropical{T}) where {T<:MM} = Tropical{T}(x.v + x.v)
Base.round(x::Tropical{T}, n::RoundingMode) where {T<:MM} = Tropical{T}(round(x.v, n))
Base.big(x::Tropical) = Base.big(x.v)
LinearAlgebra.tr(A::TropicalArray) = isempty(A) ? mp0 : sum(diag(A))
LinearAlgebra.tr(S::TropicalSparseMatrix) = isempty(S) ? mp0 : sum(diag(S))
LinearAlgebra.norm(S::TropicalArray) = MaxPlus(S[argmax(S)].v - S[argmin(S)].v) # FIXME MaxPlus + factorise
LinearAlgebra.norm(S::TropicalSparseMatrix) = MaxPlus(S[argmax(S)].v - S[argmin(S)].v)
eye(::Type{Tropical{T}}, n::Int64) where {T<:MM} = Matrix{Tropical{T}}(I, n, n)
eye(::Type{Tropical{T}}, m::Int64, n::Int64) where {T<:MM} = Matrix{Tropical{T}}(I, m, n)
plustimes(n::Tropical) = n.v
plustimes(A::TropicalArray) = map(x -> x.v, A)
plustimes(S::TropicalSparseMatrix) = sparse_map(x -> x.v, S)
full(S::TropicalSparseMatrix) = Matrix(S)
dense(S::TropicalSparseMatrix) = Matrix(S)
array(S::TropicalSparseMatrix) = Matrix(map(x -> x.v, S))

const global mp0 = zero(MaxPlus)
const global mp1 = one(MaxPlus)
const global mptop = MaxPlus(Inf)
const global mi0 = zero(MinPlus)
const global mi1 = one(MinPlus)
const global mitop = MinPlus(-Inf)

squared_size(A::AbstractArray) = (n = size(A, 1); (n != size(A, 2)) && error("Matrix shall be squared") || n)
sparse_map(f, S::SparseMatrixCSC) = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(f, S.nzval))
@inline diag_map!(f, A::AbstractArray) = for i = 1:size(A,1) A[i,i] = f(A[i,i]) end

function star(A::Array{MaxPlus})
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

star(x::MaxPlus) = star(hcat(x))[1,1]
plus(A::Array{MaxPlus}) = A * star(A)
plus(x::MaxPlus) = plus(hcat(x))[1,1]
astarb(A::ArrMaxPlus, b::ArrMaxPlus) = star(A) * b
