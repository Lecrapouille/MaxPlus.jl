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
    mp0, ε,  mp1, mpe, mi0, mi1, mptop, mitop, mpI,
    eye, full, dense, array, plustimes, inv, sparse_map,
    tr, norm, astarb, star, plus, howard, mp_change_display,
    mpshow, LaTeX

export # Max-Plus Linear system
    MPSysLin, mpsimul, mpexplicit

export # Max-Plus flowhop
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

include("types.jl")
include("fallbacks.jl")
include("show.jl")

# Constructor with conversion from Bool into a (max,+) or (min,+). Bool is not
# in tropical algebra but a hack to support the Julia generic identity matrix
# operator I, replacing the Matlab eye() function, and that uses internaly Bool.
Trop{T}(x::Bool) where {T<:MM} = x ? one(Trop{T}) : zero(Trop{T})

# Constructor with "tropicality" (max,+) from/to (min,+) conversion
Trop{T}(x::Trop) where {T<:MM} = Trop{T}(x.λ)

# Constructor with conversion of Float scalar in classic (+,*) algebra into a
# (max,+) or (min,+) number. The scalar can be a 'Not a Number' value.
Trop(::Type{T}, n::Float64) where T = isnan(n) ? zero(T) : T(n)

# Constructor converting a dense vector or matrix from classic (+,*) algebra
# into a (max,+) or (min,+) dense vector or matrix.
Trop{T}(A::Array) where {T<:MM} = map(Trop{T}, A)

# Constructor converting a sparse vector or matrix from classic (+,*) algebra
# into a (max,+) or (min,+) sparse vector or matrix.
Trop{T}(S::SparseMatrixCSC{T1,U}) where {T<:MM,T1,U} = convert(SpaTrop{T,U}, S)
Trop{T}(V::SparseVector{T1}) where {T<:MM,T1} = convert(SparseVector{Trop{T}}, V)
Trop{T}(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {T<:MM,Tv,Ti<:Integer} = sparse(I, J, Trop{T}(V)) # FIXME: to be removed

# Constructor converting a range (i.e [1.. 3]) from classic (+,*) algebra into
# (max,+) or (min,+) range.
Trop{T}(x::UnitRange) where {T<:MM} = Vector{Trop{T}}(x)
Trop{T}(x::StepRangeLen) where {T<:MM} = Vector{Trop{T}}(x)

# Promotion: any classic Julia Int or Float values are implicitly converted into
# a (max,+) or (min,+) number since we are suppose to work either in (max,+) or
# in (min,+) algebra and no more classic (+,*) algebra. Therefore for simplying
# the writting/reading of the code source, any Julia numbers are considered as
# tropical numbers to the current algebra. The context of the algebra shall be
# known. Therefore, this toolbox forbids mixture of operations on (max,+) and
# (min,+) numbers by throwing an error. This will force the developper commuting
# of algebra. Some calculus are easier made depending on the tropical algebra and
# the context may switch.
Base.promote_rule(::Type{Trop{T2}}, ::Type{T1}) where {T1,T2} = Trop{T2}
Base.promote_rule(::Type{MP}, ::Type{MI}) = error("Cannot promote (max,+) to (min,+)")
Base.promote_rule(::Type{MI}, ::Type{MP}) = error("Cannot promote (min,+) to (max,+)")
Base.convert(::Trop{T}, x::Number) where {T<:MM} = Trop{T}(x)

# Convert a (max,+) or (min,+) scalar into classic (+,*) algebra.
Base.float(x::Trop) = x.λ
plustimes(n::Trop) = n.λ

# Map operator for sparse matrix or vector of (max,+) or (min,+) numbers.
# Introduce because map operator on sparse seemed missing.
sparse_map(f, S::SparseMatrixCSC) = SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(f, S.nzval))

# Convert a (max,+) or (min,+) vector or matrix, dense or sparse into its
# equivalent in the classic (+,*).
plustimes(A::ArrTrop) = map(x -> x.λ, A)
plustimes(S::SpaTrop) = sparse_map(x -> x.λ, S)
full(S::ArrTrop) = Matrix(S)
full(S::SpaTrop) = Matrix(S)
dense(S::SpaTrop) = Matrix(S)
array(S::SpaTrop) = Matrix(map(x -> x.λ, S))

# Definition of neutral and absorbing elements for (max,+) and (min,+).
Base.zero(::Type{MP}) = MP(typemin(Float64))
Base.zero(::Type{MI}) = MI(typemax(Float64))
Base.zero(x::Trop) = zero(typeof(x))
Base.one(x::Trop) = one(typeof(x))
Base.one(::Type{Trop{T}}) where {T<:MM} = Trop{T}(zero(Float64))

# Shorter names for neutral and absorbing elements
const global mp0 = zero(MP)
const global ε = zero(MP)
const global mp1 = one(MP)
const global mpe = one(MP)
const global mptop = MP(Inf)
const global mi0 = zero(MI)
const global mi1 = one(MI)
const global mitop = MI(-Inf)

# Define the (max,+) and (min,+) alegbra operators. Since any Julia number shall
# be considered as number of the current tropical, we have to force the implicit
# conversion.
Base.:(+)(x::MP, y::MP) = MP(max(x.λ, y.λ))
Base.:(+)(x::MP, y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP) = MP(max(x, y.λ))
Base.:(+)(x::MI, y::MI) = MI(min(x.λ, y.λ))
Base.:(+)(x::MI, y::Real) = MI(min(x.λ, y))
Base.:(+)(x::Real, y::MI) = MI(min(x, y.λ))
Base.:(*)(x::MP, y::MP) = Trop(MP, x.λ + y.λ)
Base.:(*)(x::MP, y::Real) = Trop(MP, x.λ + y)
Base.:(*)(x::Real, y::MP) = Trop(MP, x + y.λ)
Base.:(*)(x::MI, y::MI) = Trop(MI, x.λ + y.λ)
Base.:(*)(x::MI, y::Real) = Trop(MI, x.λ + y)
Base.:(*)(x::Real, y::MI) = Trop(MI, x + y.λ)
Base.:(-)(x::MP, y::MP) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::MP, y::Real) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::Real, y::MP) = error("Minus operator does not exist in (max,+) algebra")
Base.:(-)(x::MI, y::MI) = error("Minus operator does not exist in (min,+) algebra")
Base.:(-)(x::MI, y::Real) = error("Minus operator does not exist in (min,+) algebra")
Base.:(-)(x::Real, y::MI) = error("Minus operator does not exist in (min,+) algebra")

# Comparaison operations. Since any Julia number shall be considered as number
# of the current tropical, we have to force the implicit conversion.
Base.:(==)(x::Trop{T}, y::Trop{T}) where {T<:MM} = (x.λ == y.λ)
Base.:(==)(x::Trop, y::Real) = (x.λ == y)
Base.:(==)(x::Real, y::Trop) = (x == y.λ)
Base.:(<=)(x::Trop, y::Trop) = (x.λ <= y.λ)
Base.:(<=)(x::Trop, y::Real) = (x.λ <= y)
Base.:(<=)(x::Real, y::Trop) = (x <= y.λ)
Base.:(<)(x::Trop, y::Trop) = (x.λ < y.λ)
Base.:(<)(x::Trop, y::Real) = (x.λ < y)
Base.:(<)(x::Real, y::Trop) = (x < y.λ)
Base.isless(x::Trop, y::Trop) = x.λ < y.λ
Base.isless(x::Trop, y::Real) = x.λ < y
Base.isless(x::Real, y::Trop) = x < y.λ

# Power, division, residuation
Base.literal_pow(::typeof(^), x::Trop, ::Val{0}) = one(x)
Base.:(^)(x::Trop{T}, y::Int) where {T<:MM} = Trop{T}(x.λ * y)
Base.:(^)(x::Trop{T}, y::Float64) where {T<:MM} = Trop{T}(x.λ * y)
Base.:(-)(x::Trop{T}) where {T<:MM} = (x == zero(Trop{T})) ? zero(Trop{T}) : Trop{T}(-x.λ)
Base.sign(x::Trop) = Base.sign(x.λ)
Base.inv(x::Trop{T}) where {T<:MM} = Trop{T}(-x.λ)
Base.inv(A::ArrMP) = Matrix(transpose(-A)) # TODO error if not(A * tr(A') == tr(A') * A == Id)
Base.inv(S::SpaMP) = SparseMatrixCSC(transpose(-S)) # TODO idem
Base.:(\)(x::Trop{T}, y::Trop{T}) where {T<:MM} = inv(inv(y) * x)
Base.:(\)(A::ArrTrop{T}, B::ArrTrop{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(\)(A::ArrTrop{T}, B::Trop{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(\)(A::Trop{T}, B::ArrTrop{T}) where {T<:MM} = inv(inv(B) * A)
Base.:(/)(x::Trop{T}, y::Trop{T}) where {T<:MM} = inv(y * inv(x))
Base.:(/)(A::ArrTrop{T}, B::ArrTrop{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::ArrTrop{T}, B::Trop{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::Trop{T}, B::ArrTrop{T}) where {T<:MM} = inv(B * inv(A))
Base.:(/)(A::SpaTrop{T,U}, B::ArrTrop{T}) where {T<:MM,U} = inv(B * inv(A)) #FIXME
Base.:min(x::MP, y::MP) = (x * y) / (x + y)
Base.:min(x::MP, y::Real) = min(x, MP(y))
Base.:min(x::Real, y::MP) = min(MP(x), y)
#FIXME Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
#FIXME Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

# Compatibility with Julia builtin operators.
Base.isnan(x::Trop) = isnan(x.λ)
Base.abs(x::Trop{T}) where {T<:MM} = Trop{T}(abs(x.λ))
Base.abs2(x::Trop{T}) where {T<:MM} = Trop{T}(x.λ + x.λ)
Base.round(x::Trop{T}, n::RoundingMode) where {T<:MM} = Trop{T}(round(x.λ, n))
Base.big(x::Trop) = Base.big(x.λ)

# Since Julia has remove the Matlab eye() function, making identity matrix, now
# replaced by the Julia I operator, we reintroduce it.
eye(::Type{Trop{T}}, n::Int64) where {T<:MM} = Matrix{Trop{T}}(I, n, n)
eye(::Type{Trop{T}}, m::Int64, n::Int64) where {T<:MM} = Matrix{Trop{T}}(I, m, n)
eye(A::Array{Trop{T},N}) where {T,N} = Array{Trop{T},N}(mpI, size(A,1), size(A,2))

# Julia builtin zero, zeros, spzero, and their equivalent functions for ones are
# defacto good for building scalar and matrices. We complete some functions.
Base.zeros(A::Array{Trop{T},N}) where {T,N} = zeros(Trop{T}, size(A,1), size(A,2))
SparseArrays.spzeros(A::Array{Trop{T},N}) where {T,N} = spzeros(Trop{T}, size(A,1), size(A,2))
Base.ones(A::Array{Trop{T},N}) where {T,N} = ones(Trop{T}, size(A,1), size(A,2))

# Equivalent to I but for (max,+) and (min,+)
const global mpI = UniformScaling(one(MP))
const global miI = UniformScaling(one(MI))

# Trace and norm. FIXME: I do not remeber why we have to redefine these
# functions and why they are not directly working !
LinearAlgebra.tr(A::ArrTrop) = isempty(A) ? mp0 : sum(diag(A))
LinearAlgebra.tr(S::SpaTrop) = isempty(S) ? mp0 : sum(diag(S))
LinearAlgebra.norm(S::ArrTrop) = MP(S[argmax(S)].λ - S[argmin(S)].λ) # FIXME MP + factorise
LinearAlgebra.norm(S::SpaTrop) = MP(S[argmax(S)].λ - S[argmin(S)].λ)

# Do the A*, A+, Ax+b operations
squared_size(A::AbstractArray) = (n = size(A, 1); (n != size(A, 2)) && error("Matrix shall be squared") || n)
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
