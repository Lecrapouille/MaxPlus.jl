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
Trop{T}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MM,Tv,Ti<:Integer} = sparse(I, J, Trop{T}(V))
Trop{T}(x::UnitRange) where {T<:MM} = Vector{Trop{T}}(x)
Trop{T}(x::StepRangeLen) where {T<:MM} = Vector{Trop{T}}(x)
Base.promote_rule(::Type{Trop{T2}}, ::Type{T}) where {T,T2} = Trop{T2}
Base.promote_rule(::Type{MP}, ::Type{MI}) = error("Cannot promote Max-Plus to Min-Plus")
Base.promote_rule(::Type{MI}, ::Type{MP}) = error("Cannot promote Min-Plus to Max-Plus")
Base.convert(::Trop{T}, x::Number) where {T<:MM} = Trop{T}(x)
Base.float(x::Trop) = x.v
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
Base.:(*)(x::MI, y::MI) = MInan(x.v + y.v)
Base.:(*)(x::MI, y::Real) = MInan(x.v + y)
Base.:(*)(x::Real, y::MI) = MInan(x + y.v)
Base.:(^)(x::Trop{T}, y::Int) where {T<:MM} = Trop{T}(x.v * y)
Base.:(^)(x::Trop{T}, y::Float64) where {T<:MM} = Trop{T}(x.v * y)
Base.:(-)(x::Trop{T}) where {T<:MM} = (x == zero(Trop{T})) ? zero(Trop{T}) : Trop{T}(-x.v)
Base.sign(x::Trop) = Base.sign(x.v)
Base.inv(A::ArrMP) = MI(transpose(-A))
Base.:(/)(x::MP, y::MP) = Trop(MP, x.v - y.v)
Base.:(/)(A::ArrMP, B::ArrMP) = MP(inv(A) * MI(B))
Base.:(\)(x::MP, y::MP) = Trop(MP, y.v - x.v)
Base.:(\)(A::ArrMP, B::ArrMP) = MI(MP(A) * inv(B))
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
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)
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
    s = Int(round(log2(squared_size(A))))
    M = A
    f = x -> x < 0.0 ? mp1 : x > 0.0 ? mptop : 0.0
    diag_map!(f, M)
    for i = 0:s
        diag_map!(f, M)
        M *= M
    end
    M
end

star(x::MP) = star(reshape([x], :, 1))[1,1]



# A^+
function plus(A::ArrMP)
    n = size(A, 1)
    (n != size(A, 2)) && error("Matrix shall be squared")
    C = A
    for k in 1:n
        t = (C[k,k].v <= zero(Float64)) ? zero(Float64) : typemax(Float64);
        for j in 1:n, i in 1:n
            C[i,j] = MP(max(C[i,j].v, C[i,k].v + C[k,j].v + t))
        end
    end
    C
end

plus(x::MP) = plus([x])[1,1]

# A^* b
mpastarb(A::ArrMP, b::ArrMP) = mpstar(A) * b

include("howard.jl")
include("syslin.jl")
#include("flowshop.jl")
include("docstrings.jl")

# Map function f to a sparse matrice


end # MaxPlus module
