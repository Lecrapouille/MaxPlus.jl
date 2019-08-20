# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# In the way of the Scicoslab max-plus toolbox.
# ==============================================================================

module MaxPlus

using
    LinearAlgebra, SparseArrays, Printf

export
    MP, ⨁, ⨂,
    mpzero, mpone, mp0, mp1, mptop,
    mpI, mpeye, mpzeros, mpones,
    mpsparse, full, dense, plustimes, array, mparray, mptrace, merde

# ==============================================================================

"""
    MP{T}

Max-Plus immutable structure.
Promote a number of type T (such as ideally Float64) to a
max-plus number.

# Examples
```julia-repl
a=MP(3)
typeof(a)
```
"""
struct MP{T} <: Number λ::T end

# ==============================================================================

# Copy constructor
MP(x::MP) = MP(x.λ)

# ==============================================================================

"""
    show(io::IO,a::MP)

Display a max-plus number. -Inf are displayed with 'ε' symbols.
Note in Scicoslab -Inf are displayed with '.' symbols.

# Examples
```julia-repl
julia> mpeye(Float64, 5,5)
5×5 Array{MP{Float64},2}:
 0.000      ε      ε      ε      ε
     ε  0.000      ε      ε      ε
     ε      ε  0.000      ε      ε
     ε      ε      ε  0.000      ε
     ε      ε      ε      ε  0.000
```
"""
Base.show(io::IO,a::MP) =
    (a.λ == -Inf) ? (@printf io "ε") : show(io, a.λ)

#Base.show(io::IO,a::MP) =
#    (a.λ == -Inf) ? (@printf io "ε") : ((a.λ == 0) ? (@printf io "e") : show(io, a.λ))

# ==============================================================================

"""
    +(a::MP, b::MP)

Max operator. Return the maximum of `a` and `b`.
See also ⨁ operator.

# Examples
```julia-repl
julia> MP(1.0) + MP(3.0)
MP{Float64}(3.0)
```
"""
Base.:+(a::MP,   b::MP)   = MP(max(a.λ, b.λ))
Base.:+(a::MP,   b::Real) = MP(max(a.λ, b))
Base.:+(a::Real, b::MP)   = MP(max(a,   b.λ))

# ==============================================================================

"""
    ⨁(a::MP, b::MP)

Max operator. Return the maximum of `a` and `b`.
Unicode character "circled plus": U+2A01.

# Examples
```julia-repl
julia> ⨁(1.0, 3.0)
MP{Float64}(3.0)
```
"""
⨁(a::MP,   b::MP)   = MP(max(a.λ, b.λ))
⨁(a::MP,   b::Real) = MP(max(a.λ, b))
⨁(a::Real, b::MP)   = MP(max(a,   b.λ))
⨁(a::Real, b::Real) = MP(max(a,   b))

# ==============================================================================

"""
    *(a::MP, b::MP)

Addition operator. Return the summation of `a` and `b`.

# Examples
```julia-repl
julia> MP(1.0) * MP(3.0)
MP{Float64}(4.0)
```
"""
Base.:*(a::MP,   b::MP)   = MP(a.λ + b.λ)
Base.:*(a::MP,   b::Real) = MP(a.λ + b)
Base.:*(a::Real, b::MP)   = MP(a   + b.λ)

# ==============================================================================

"""
    ⨂(a::MP, b::MP)

Return the summation of `a` and `b`.
Unicode character "circled times": U+2A00.

# Examples
```julia-repl
julia> ⨂(1.0, 3.0)
MP{Float64}(4.0)
```
"""
⨂(a::MP,   b::MP)   = MP(a.λ + b.λ)
⨂(a::MP,   b::Real) = MP(a.λ + b)
⨂(a::Real, b::MP)   = MP(a   + b.λ)
⨂(a::Real, b::Real) = MP(a   + b)

# ==============================================================================

"""
    -(a::MP, b::MP)

Minus operator. Return the difference between `a` and `b`.

# Examples
```julia-repl
julia> MP(1.0) / MP(2.0)
MP{Float64}(-1.0)
```
"""
Base.:-(a::MP,   b::MP)   = MP(a.λ - b.λ)
Base.:-(a::MP,   b::Real) = MP(a.λ - b)
Base.:-(a::Real, b::MP)   = MP(a   - b.λ)

# ==============================================================================

"""
    /(a::MP, b::MP)

Divisor operator. Return the difference between `a` and `b`.

# Examples
```julia-repl
julia> MP(1.0) / MP(2.0)
MP{Float64}(-1.0)
```
"""
Base.:/(a::MP,   b::MP)   = MP(a.λ - b.λ)
Base.:/(a::MP,   b::Real) = MP(a.λ - b)
Base.:/(a::Real, b::MP)   = MP(a   - b.λ)

# ==============================================================================

"""
    min(a::MP,   b::MP)

Return the minimun of `a` and `b`.

# Examples
```julia-repl
julia> min(MP(1), -3)
MP{Int64}(-3)
```
"""
Base.:min(a::MP,   b::MP)   = MP(min(a.λ, b.λ))
Base.:min(a::MP,   b::Real) = MP(min(a.λ, b))
Base.:min(a::Real, b::MP)   = MP(min(a,   b.λ))

# ==============================================================================

Base.isless(a::MP{T}, b::MP{T}) where T = a.λ < b.λ
Base.promote_rule(::Type{MP{T}}, ::Type{U}) where T where U = MP{T}
Base.convert(::MP{T}, x) where T = MP(T(x))

# ==============================================================================

"""
    zero(::MP{T})

Create the max-plus constant zero equal to -Inf.
This value is the neutral for the ⨁ operator.
"""
Base.zero(::MP{T})       where T = MP(typemin(T))
Base.zero(::Type{MP{T}}) where T = MP(typemin(T))

"""
    mpzero(::Type{T})

Create the max-plus constant zero equal to -Inf.
This value is the neutral for the ⨁ operator.
"""
mpzero(::Type{T})        where T = MP(typemin(T))

"""
    one(::MP{T})

Create the max-plus constant one equal to 0.
This value is the neutral for operators ⨁ and ⨂.
"""
Base.one(::MP{T})        where T = MP(zero(T))
Base.one(::Type{MP{T}})  where T = MP(zero(T))

"""
    one(::Type{T})

Create the max-plus constant one equal to 0.
This value is the neutral for operators ⨁ and ⨂.
"""
mpone(::Type{T})         where T = MP(zero(T))

# ==============================================================================

"""
    mp0

Max-plus constant zero equal to -Inf.
This value is the neutral for the ⨁ operator.

Equivalent to Scicoslab code: `%0 = #(-%inf) = .`

# Examples
```julia-repl
julia> mp0 * 5
-Inf

julia> mp0 + 5
5.0
```
"""
mp0 = mpzero(Float64)

# ==============================================================================

"""
    mp1

Max-plus constant one equal to 0.
This value is the neutral for operators ⨁ and ⨂.

Equivalent to Scicoslab code: `%1 = #(1) = 0`

# Examples
```julia-repl
julia> mp1 * 5
5.0

julia> mp1 + 5
5.0
```
"""
mp1 = mpone(Float64)

# ==============================================================================

"""
    mptop

Max-plus constant "top" equal to +Inf.

Equivalent to Scicoslab code: `%top = +Inf`
"""
mptop = MP{Float64}(Inf)

# ==============================================================================

"""
    mpI

Fix an algebra conception in Julia official LinearAlgebra (uniformscaling.jl)
`I` representes an identity matrix of any size is defined with a booleen instead
of the function `one()`. As consequence, in Julia 0.4 the `eye(T,m,n)` could
created a max-plus identity matrix. Since Julia 0.7 `eye()` has been deprecated
and replaced by the buggy `Matrix{T}(I, m, n)`. This function uses `zero()` but
not `one()` and as consequence the max-plus identity matrix is not well formed.

This const allows to be more algebra compliant by calling `one()` and fixing the
fucntion `Matrix{T}(I, m, n)`.
"""
const mpI = UniformScaling(one(MP{Float64}).λ)

"""
    mpeye(::Type{T}, n::Int64)

Construct a max-plus identity n-by-n matrix.

# Examples
```julia-repl
julia> mpeye(Float64, 2)
2×2 Array{MP{Float64},2}:
  0.0  -Inf
 -Inf   0.0
```
"""
mpeye(::Type{T}, n::Int64) where T = Matrix{MP{T}}(mpI, n, n)

# ==============================================================================

"""
    mpeye(::Type{T}, m::Int64, n::Int64)

Construct a max-plus identity m-by-n matrix.

# Examples
```julia-repl
julia> mpeye(Float64, 2, 2)
2×2 Array{MP{Float64},2}:
  0.0  -Inf
 -Inf   0.0
```
"""
mpeye(::Type{T}, m::Int64, n::Int64) where T = Matrix{MP{T}}(mpI, m, n)

# ==============================================================================

"""
    mpzeros(::Type{T}, n::Int64)

Construct a max-plus zero n-by-m sparse matrix.

# Examples
```julia-repl
julia> mpzeros(Float64, 2)
2-element SparseVector{MP{Float64},Int64} with 0 stored entries
```
"""
mpzeros(::Type{T}, n::Int64) where T = spzeros(MP{T}, n)

# ==============================================================================

"""
    mpzeros(::Type{T}, n::Int64)

Construct a sparse max-plus zero m-by-n matrix.

# Examples
```julia-repl
julia> mpzeros(Float64, 2,5)
2×5 SparseMatrixCSC{MP{Float64},Int64} with 0 stored entries
```
"""
mpzeros(::Type{T}, m::Int64, n::Int64) where T = spzeros(MP{T}, m, n)

# ==============================================================================

"""
    mpones(::Type{T}, n::Int64)

Construct a max-plus one n-by-n matrix.

# Examples
```julia-repl
julia> mpones(Float64,2)
2-element Array{MP{Float64},1}:
 0.0
 0.0
```
"""
mpones(::Type{T}, n::Int64) where T = ones(MP{T}, n)

# ==============================================================================

"""
    mpones(::Type{T}, m::Int64, n::Int64)

Construct a max-plus one m-by-n matrix.

# Examples
```julia-repl
julia> mpones(Float64,2)
2×2 Array{MP{Float64},2}:
 0.0  0.0
 0.0  0.0
```
"""
mpones(::Type{T}, m::Int64, n::Int64) where T = ones(MP{T}, m, n)

# ==============================================================================

"""
    mpsparse(A::Array{T})

Transform a dense matrix to a max-plus sparse matrix.
"""
mpsparse(A::Array{T})     where T = mpsparse(mparray(A))
function mpsparse(M::Array{MP{T}}) where T
    R = Array{Int64, 1}[];
    C = Array{Int64, 1}[];
    V = Array{MP{T}, 1}[];

    for i = 1:size(M, 1), j=1:size(M, 2)
        M[i,j] != zero(MP{T}) && (R = [R; i]; C = [C; j]; V = [V; M[i,j]])
    end

    sparse(R, C, V)
end

# ==============================================================================

"""
    array(::Array{MP{T}})

Convert a max-plus array to a standard type array.

# Examples
```julia-repl
julia> array([MP(1.0) 2.0; 3.0 4.0])
2×2 Array{Float64,2}:
 1.0  2.0
 3.0  4.0
```
"""
array(A::Array{MP{T}}) where T = map(x -> x.λ, A)

# ==============================================================================

"""
    array(::SparseMatrixCSC{MP{T}})

Convert a sparse max-plus array to a dense non max-plus array.

# Examples
```julia-repl
julia> array(mpzeros(Float64, 2,2))
2×2 SparseMatrixCSC{Float64,Int64} with 4 stored entries:
  [1, 1]  =  -Inf
  [2, 1]  =  -Inf
  [1, 2]  =  -Inf
  [2, 2]  =  -Inf
```
"""
array(A::SparseMatrixCSC{MP{T},Int64}) where T = map(x -> x.λ, A)

# ==============================================================================

"""
    full(::SparseMatrixCSC{MP{T}})

Convert a sparse max-plus array to a dense max-plus array.
Alternative function name: dense.

# Examples
```julia-repl
julia> full(mpzeros(Float64, 2,5))
2×5 Array{MP{Float64},2}:
 -Inf  -Inf  -Inf  -Inf  -Inf
 -Inf  -Inf  -Inf  -Inf  -Inf
```
"""
full(A::SparseMatrixCSC{MP{T},Int64}) where T = Matrix(A)

# ==============================================================================

"""
    dense(::SparseMatrixCSC{MP{T}})

Convert a sparse max-plus array to a dense max-plus array.
Alternative function name: full.

# Examples
```julia-repl
julia> dense(mpzeros(Float64, 2,5))
2×5 Array{MP{Float64},2}:
 -Inf  -Inf  -Inf  -Inf  -Inf
 -Inf  -Inf  -Inf  -Inf  -Inf
```
"""
dense(A::SparseMatrixCSC{MP{T},Int64}) where T = Matrix(A)

# ==============================================================================

"""
    plustimes(a::MP{T})

Convert max-plus number to the standard type.
Note: the function name comes from Scicoslab.

# Examples
```julia-repl
julia> typeof(plustimes(MP(1)))
Int64

julia> plustimes([MP(1.0) 2.0; 3.0 4.0])
2×2 Array{Float64,2}:
 1.0  2.0
 3.0  4.0
```
"""
plustimes(a::MP{T})        where T = a.λ
plustimes(a::Array{MP{T}}) where T = array(a)

# ==============================================================================

"""
    mparray(A::Array)

Convert an array to a max-plus array.

# Examples
```julia-repl
julia> mparray([1.0 2; 3 4])
2×2 Array{MP{Float64},2}:
 1.0  2.0
 3.0  4.0
```
"""
mparray(A::Array)              = map(MP, A)
# FIXME missing mparray([1, 2, 3], [1, 2, 3], [0, 2, 0])

# ==============================================================================

"""
    mptrace(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> mptrace([MP(1) 2; 3 4])
MP{Float64}(4.0)
```
"""
mptrace(A::Array{MP{T}}) where T = sum(diag(A))
mptrace(A::Array{T}) where T = sum(diag(mparray(A)))
mptrace(A::SparseMatrixCSC{MP{T},Int64}) where T = sum(diag(A))
mptrace(A::SparseMatrixCSC{T,Int64}) where T = sum(diag(A))

end # MaxPlus module
