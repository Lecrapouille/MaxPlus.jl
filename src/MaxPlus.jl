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
    mp0, mp1, mi0, mi1, mptop, mitop,
    eye, full, dense, array, plustimes, inv,
    tr, norm, mpastarb, star, plus, howard, mp_change_display,
    mpshow, LaTeX

export # Max-Plus Linear system
    MPSysLin, mpsimul, mpexplicit

export # Max-Plus flowhop
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

include("types.jl")
include("fallbacks.jl")
include("show.jl")

Trop{T}(x::Bool) where {T<:MM} = x ? one(Trop{T}) : zero(Trop{T})
Trop{T}(x::Trop) where {T<:MM} = Trop{T}(x.v)
Trop(::Type{T}, n::Float64) where T = isnan(n) ? zero(T) : T(n)
Trop{T}(A::Array) where {T<:MM} = map(Trop{T}, A)
Trop{T}(S::SparseMatrix{T1,U}) where {T<:MM,T1,U} = convert(SpaTrop{T,U}, S)
Trop{T}(V::SparseVector{T1}) where {T<:MM,T1} = convert(SparseVector{Trop{T}}, V)
Trop{T}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MM,Tv,Ti<:Integer} = sparse(I, J, Trop{T}(V)) # FIXME: to be removed
Trop{T}(x::UnitRange) where {T<:MM} = Vector{Trop{T}}(x)
Trop{T}(x::StepRangeLen) where {T<:MM} = Vector{Trop{T}}(x)
Base.promote_rule(::Type{Trop{T2}}, ::Type{T1}) where {T1,T2} = Trop{T2}
Base.promote_rule(::Type{MP}, ::Type{MI}) = error("Cannot promote Max-Plus to Min-Plus")
Base.promote_rule(::Type{MI}, ::Type{MP}) = error("Cannot promote Min-Plus to Max-Plus")
Base.convert(::Trop{T}, x::Number) where {T<:MM} = Trop{T}(x)
Base.float(x::Trop) = x.v
Base.isnan(x::Trop) = isnan(x.v)
Base.zero(::Type{MP}) = MP(typemin(Float64))
Base.zero(::Type{MI}) = MI(typemax(Float64))
Base.zero(x::Trop) = zero(typeof(x))
Base.one(x::Trop) = one(typeof(x))
Base.one(::Type{Trop{T}}) where {T<:MM} = Trop{T}(zero(Float64))
Base.:(+)(x::MP, y::MP) = MP(max(x.v, y.v))
Base.:(+)(x::MP, y::Real) = MP(max(x.v, y))
Base.:(+)(x::Real, y::MP) = MP(max(x, y.v))
Base.:(+)(x::MI, y::MI) = MI(min(x.v, y.v))
Base.:(+)(x::MI, y::Real) = MI(min(x.v, y))
Base.:(+)(x::Real, y::MI) = MI(min(x, y.v))
Base.:(*)(x::MP, y::MP) = Trop(MP, x.v + y.v)
Base.:(*)(x::MP, y::Real) = Trop(MP, x.v + y)
Base.:(*)(x::Real, y::MP) = Trop(MP, x + y.v)
Base.:(*)(x::MI, y::MI) = Trop(MI, x.v + y.v)
Base.:(*)(x::MI, y::Real) = Trop(MI, x.v + y)
Base.:(*)(x::Real, y::MI) = Trop(MI, x + y.v)
Base.:(-)(x::MP, y::MP) = error("Minus operator does not exist in max+ algebra")
Base.:(-)(x::MP, y::Real) = error("Minus operator does not exist in max+ algebra")
Base.:(-)(x::Real, y::MP) = error("Minus operator does not exist in max+ algebra")
Base.:(-)(x::MI, y::MI) = error("Minus operator does not exist in min+ algebra")
Base.:(-)(x::MI, y::Real) = error("Minus operator does not exist in min+ algebra")
Base.:(-)(x::Real, y::MI) = error("Minus operator does not exist in min+ algebra")
Base.literal_pow(::typeof(^), x::Trop, ::Val{0}) = one(x)
Base.:(^)(x::Trop{T}, y::Int) where {T<:MM} = Trop{T}(x.v * y)
Base.:(^)(x::Trop{T}, y::Float64) where {T<:MM} = Trop{T}(x.v * y)
Base.:(-)(x::Trop{T}) where {T<:MM} = (x == zero(Trop{T})) ? zero(Trop{T}) : Trop{T}(-x.v)
Base.sign(x::Trop) = Base.sign(x.v)
Base.inv(x::Trop{T}) where {T<:MM} = Trop{T}(-x.v)
Base.inv(A::ArrMP) = Matrix(transpose(-A))
Base.inv(S::SpaMP) = SparseMatrixCSC(transpose(-S))
Base.:(\)(x::Trop{T}, y::Trop{T}) where {T<:MM} = inv(inv(y) * x)
Base.:(\)(A::ArrTrop{T}, B::ArrTrop{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(\)(A::ArrTrop{T}, B::Trop{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(\)(A::Trop{T}, B::ArrTrop{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(/)(x::Trop{T}, y::Trop{T}) where {T<:MM} = inv(y * inv(x))
Base.:(/)(A::ArrTrop{T}, B::ArrTrop{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::ArrTrop{T}, B::Trop{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::Trop{T}, B::ArrTrop{T}) where {T<:MM} = inv(B * inv(A))

#FIXME
Base.:(/)(A::SpaTrop{T,U}, B::ArrTrop{T}) where {T<:MM,U} = inv(B * inv(A))

Base.:(==)(x::Trop{T}, y::Trop{T}) where {T<:MM} = (x.v == y.v)
Base.:(==)(x::Trop, y::Real) = (x.v == y)
Base.:(==)(x::Real, y::Trop) = (x == y.v)
Base.:(<=)(x::Trop, y::Trop) = (x.v <= y.v)
Base.:(<=)(x::Trop, y::Real) = (x.v <= y)
Base.:(<=)(x::Real, y::Trop) = (x <= y.v)
Base.:(<)(x::Trop, y::Trop) = (x.v < y.v)
Base.:(<)(x::Trop, y::Real) = (x.v < y)
Base.:(<)(x::Real, y::Trop) = (x < y.v)
Base.isless(x::Trop, y::Trop) = x.v < y.v
Base.isless(x::Trop, y::Real) = x.v < y
Base.isless(x::Real, y::Trop) = x < y.v
Base.:min(x::MP, y::MP) = (x * y) / (x + y)
Base.:min(x::MP, y::Real) = min(x, MP(y))
Base.:min(x::Real, y::MP) = min(MP(x), y)
#Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
#Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)
Base.abs(x::Trop{T}) where {T<:MM} = Trop{T}(abs(x.v))
Base.abs2(x::Trop{T}) where {T<:MM} = Trop{T}(x.v + x.v)
Base.round(x::Trop{T}, n::RoundingMode) where {T<:MM} = Trop{T}(round(x.v, n))
Base.big(x::Trop) = Base.big(x.v)
LinearAlgebra.tr(A::ArrTrop) = isempty(A) ? mp0 : sum(diag(A))
LinearAlgebra.tr(S::SpaTrop) = isempty(S) ? mp0 : sum(diag(S))
LinearAlgebra.norm(S::ArrTrop) = MP(S[argmax(S)].v - S[argmin(S)].v) # FIXME MP + factorise
LinearAlgebra.norm(S::SpaTrop) = MP(S[argmax(S)].v - S[argmin(S)].v)
eye(::Type{Trop{T}}, n::Int64) where {T<:MM} = Matrix{Trop{T}}(I, n, n)
eye(::Type{Trop{T}}, m::Int64, n::Int64) where {T<:MM} = Matrix{Trop{T}}(I, m, n)
plustimes(n::Trop) = n.v
plustimes(A::ArrTrop) = map(x -> x.v, A)
plustimes(S::SpaTrop) = sparse_map(x -> x.v, S)
full(S::SpaTrop) = Matrix(S)
dense(S::SpaTrop) = Matrix(S)
array(S::SpaTrop) = Matrix(map(x -> x.v, S))

const global mp0 = zero(MP)
const global mp1 = one(MP)
const global mptop = MP(Inf)
const global mi0 = zero(MI)
const global mi1 = one(MI)
const global mitop = MI(-Inf)

squared_size(A::AbstractArray) = (n = size(A, 1); (n != size(A, 2)) && error("Matrix shall be squared") || n)
sparse_map(f, S::SparseMatrix) = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(f, S.nzval))
@inline diag_map!(f, A::AbstractArray) = for i = 1:size(A,1) A[i,i] = f(A[i,i]) end

function star(A::Array{MP})
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
plus(A::Array{MP}) = A * star(A)
plus(x::MP) = plus(hcat(x))[1,1]
astarb(A::ArrMP, b::ArrMP) = star(A) * b

include("howard.jl")
include("syslin.jl")
#include("flowshop.jl")
include("docstrings.jl")

end # MaxPlus module
