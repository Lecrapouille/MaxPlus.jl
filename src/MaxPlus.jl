# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# In the way of the ScicosLab max-plus toolbox.
# ==============================================================================

module MaxPlus

using
    LinearAlgebra, SparseArrays, Printf

# ==============================================================================
# Max-Plus core

export
    MP, SpaMP, SpvMP, ArrMP, VecMP,
    mpzero, mpone, mp0, mp1, mptop,
    mpI, mpeye, mpzeros, mpones,
    mpsparse, full, dense, array,
    plustimes, minplus,
    mptrace, mpnorm, mpstar,
    mp_change_display, LaTeX

# ==============================================================================
# Max-Plus Linear system

export
    MPSysLin,
    mpsyslin, mpsimul, mpexplicit

# ==============================================================================
# Max-Plus flowhop

export
    mpgraph
#    flowshop, flowshop_graph, flowshop_simu

# ==============================================================================

"""
    MP{T}

Max-Plus immutable structure. Promote a number of type T (ideally a Float64) to
a max-plus number.

# Examples
```julia-repl
julia> a=MP(3.0); typeof(a)
MP{Float64}

julia> a=MP(3); typeof(a)
MP{Int64}
```
"""
struct MP{T <: Real} <: Real
    λ::T
end

# ==============================================================================
# Type alias but with shorter number of characters

const Sparse{T,U} = SparseMatrixCSC{T,U}
const SpaMP{T,U}  = SparseMatrixCSC{MP{T},U}
const SpvMP{T,U}  = SparseVector{MP{T},U}
const ArrMP{T,N}  = Array{MP{T},N}
const VecMP{T,N}  = Vector{MP{T}}

# ==============================================================================
# Copy constructor

MP(x::MP) = MP(x.λ)

# ==============================================================================

"""
    zero(::MP{T})

Create the max-plus constant zero equal to -Inf.
This value is the neutral for the ⨁ operator.
"""
Base.zero(::Type{MP{T}}) where T = MP(typemin(T))
Base.zero(x::MP{T})      where T = zero(typeof(x))

"""
    mpzero(::Type{T})

Create the max-plus constant zero equal to -Inf.
This value is the neutral for the ⨁ operator.
"""
mpzero(::Type{T})        where T = MP(typemin(T))

# ==============================================================================

"""
    one(::MP{T})

Create the max-plus constant one equal to 0.
This value is the neutral for operators ⨁ and ⨂.
"""
Base.one(::Type{MP{T}})  where T = MP(zero(T))
Base.one(x::MP{T})       where T = one(typeof(x))

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

Equivalent to ScicosLab code: `%0 = #(-%inf) = .`

# Examples
```julia-repl
julia> mp0 * 5
-Inf

julia> mp0 + 5
5.0
```
"""
const global mp0 = mpzero(Float64)

# ==============================================================================

"""
    mp1

Max-plus constant one equal to 0.
This value is the neutral for operators ⨁ and ⨂.

Equivalent to ScicosLab code: `%1 = #(1) = 0`

# Examples
```julia-repl
julia> mp1 * 5
5.0

julia> mp1 + 5
5.0
```
"""
const global mp1 = mpone(Float64)

# ==============================================================================

"""
    mptop

Min-plus constant "top" equal to +Inf.

Equivalent to ScicosLab code: `%top = #(%inf)`
"""
const global mptop = MP{Float64}(Inf)

# ==============================================================================
# Conversions

Base.promote_rule(::Type{MP{T}}, ::Type{U}) where {T, U} = MP{T}
Base.convert(::MP{T}, x::Number)            where T = MP{T}(x)

# ==============================================================================

"""
    MP(x::UnitRange)

Create a max-plus dense vector.

# Examples
```julia-repl
julia> MP(1:3)
3-element Array{MP{Int64},1}:
 1
 2
 3
```
"""
MP(x::UnitRange{T}) where T = MP(Vector{T}(x))

# ==============================================================================

"""
    MP(A::Array)

Convert a dense array to a max-plus dense array.

# Examples
```julia-repl
julia> MP([1.0 -Inf; 0.0 4])
2×2 Array{MP{Float64},2}:
 1.0  -Inf
 0.0  4.0
```
"""
MP(A::Array) = map(MP, A)

# ==============================================================================

"""
    MP(A::SparseArray)

Convert a sparse array to a max-plus sparse array.
By default values are max-plus zeros values are preserved else if
parameter `preserve` is set then

# Examples
```julia-repl
julia> MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
  [1, 1]  =  -Inf
  [2, 2]  =  2.0
  [3, 3]  =  0.0

julia> MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]), preserve=false)
  [2, 2]  =  2.0
```
"""
function MP(S::Sparse{T,U}; preserve=true) where {T, U}
    if preserve
        convert(SpaMP{T,U}, S)
    else
        M = spzeros(MP{T}, size(S,1), size(S,2))
        for c = 1:size(S, 2), r = nzrange(S, c)
            if (S[r,c] != zero(T)) && (MP(S[r,c]) != mpzero(T))
                M[r,c] = convert(MP{T}, S[r,c])
            end
        end
        M
    end
end

# ==============================================================================

"""
    MP(A::SparseVector)

Convert a sparse vector to a max-plus sparse vector.
By default values are max-plus zeros values are preserved else if
parameter `preserve` is set then

# Examples
```julia-repl
julia> MP(sparse([1.0, 0.0, 1.0]))
3-element SparseVector{MP{Float64},Int64} with 2 stored entries:
  [1]  =  1.0
  [3]  =  1.0
```
"""
function MP(V::SparseVector{T,U}; preserve=true) where {T, U}
    if preserve
        convert(SparseVector{MP{T},U}, V)
    else
        S = spzeros(MP{T}, size(V,1), size(V,2))
        for c = V.n
            if (S[c] != zero(T)) && (MP(S[c]) != mpzero(T))
                M[c] = convert(MP{T}, S[c])
            end
        end
        S
    end
end

# ==============================================================================

global mpstyle = 1 # Style ScicosLab

"""
    mp_change_display(style::Int)

Change the style of display for the Base.show(io::IO, x::MP) function.
See help of show(io::IO,x::MP) for examples.
"""
mp_change_display(style::Int) = (global mpstyle = min(max(style, 0), 4))

"""
    show(io::IO,x::MP)

Display a max-plus number depending on the currently set style:
- -Inf are displayed either with 'ε' (style 2 or 3) or '.' symbols (style 1).
- 0 are displayed either with 'e' (style 3) or '0' symbols (style 1 or 2).
- else -Inf and 0 are displayed in Julia default sytle (style 0).
By default the ScicosLab style is used (style 1).

# Examples
```julia-repl
julia> mp_change_display(0); mpeye(Float64, 2,2)
2×2 Array{MP{Float64},2}:
  0.0  -Inf
 -Inf   0.0

julia> mp_change_display(0); mpeye(Int64, 2,2)
2×2 Array{MP{Int64},2}:
                    0  -9223372036854775808
 -9223372036854775808                     0

julia> mp_change_display(1); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Array{MP{Float64},2}:
 0  .
 .  0

julia> mp_change_display(2); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Array{MP{Float64},2}:
 e  .
 .  e

julia> mp_change_display(3); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Array{MP{Float64},2}:
 0  ε
 ε  0

julia> mp_change_display(4); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Array{MP{Float64},2}:
 e  ε
 ε  e
```
"""
function Base.show(io::IO, x::MP{T}) where T
    if mpstyle == 0
        show(io, x.λ)
    elseif x == mpzero(T)
        (mpstyle == 1 || mpstyle == 2) ? (@printf io ".") : (@printf io "ε")
    elseif x == mpone(T)
        (mpstyle == 1 || mpstyle == 3) ? (@printf io "0") : (@printf io "e")
    elseif x.λ == trunc(x.λ)
        (@printf io "%d" x.λ)
    else
        show(io, x.λ)
    end
end

# ==============================================================================
# FIXME: not a column-major traversal
"""
    LaTeX(io::IO, A::Array{MP})

Convert a max-plus dense matrix to a LaTeX formula. The display depends on
mp_change_display(style).

# Examples
```julia-repl
julia> LaTeX(stdout, MP([4 3; 7 -Inf]))
\\left[
\\begin{array}{*{20}c}
4 & 3 \\\\
7 & . \\\\
\\end{array}
\\right]
```
"""
function LaTeX(io::IO, A::ArrMP{T}) where T
    (@printf io "\\left[\n\\begin{array}{*{20}c}\n")
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == mpzero(T)
                if mpstyle == 0
                    (@printf io "-\\infty")
                elseif mpstyle == 3 || mpstyle == 4
                    (@printf io "\\varepsilon")
                else
                    (@printf io ".")
                end
            elseif A[i,j] == mpone(T)
                if mpstyle == 2 || mpstyle == 4
                    (@printf io "e")
                else
                    (@printf io "%d" A[i,j].λ)
                end
            elseif A[i,j].λ == trunc(A[i,j].λ)
                (@printf io "%d" A[i,j].λ)
            else
                show(io, A[i,j].λ)
            end
            if j < size(A, 2)
                (@printf io " & ")
            end
        end
        (@printf io " \\\\\n")
    end
    (@printf io "\\end{array}\n\\right]\n")
end

LaTeX(io::IO, S::SpaMP{T}) where {T,U} = LaTeX(io, full(S))

# ==============================================================================

"""
    +(x::MP, y::MP)

Max operator ⨁. Return the maximum of `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) + MP(3.0)
MP{Float64}(3.0)
```
"""
Base.:(+)(x::MP,   y::MP)   = MP(max(x.λ, y.λ))
Base.:(+)(x::MP,   y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP)   = MP(max(x,   y.λ))

# ==============================================================================

"""
    *(x::MP, y::MP)

Addition operator ⨂. Return the summation of `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) * MP(3.0)
MP{Float64}(4.0)
```
"""
Base.:(*)(x::MP,   y::MP)   = MP(x.λ + y.λ)
Base.:(*)(x::MP,   y::Real) = MP(x.λ + y)
Base.:(*)(x::Real, y::MP)   = MP(x   + y.λ)

# ==============================================================================
# These functions fix Julia bugs

@inline Base.literal_pow(::typeof(^), x::MP, ::Val{0}) = one(x)
@inline Base.literal_pow(::typeof(^), x::MP, ::Val{p}) where {p} = MP(x.λ * p)
@inline Base.:(^)(x::MP, y::Number) = MP(x.λ * y)
@inline Base.literal_pow(::typeof(^), A::ArrMP{T}, ::Val{p}) where {T, p} = A^p
@inline Base.abs2(x::MP) = x.λ + x.λ

# ==============================================================================

"""
    /(x::MP, y::MP)

Divisor operator. Return the difference between `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) / MP(2.0)
MP{Float64}(-1.0)
```
"""
Base.:(/)(x::MP,   y::MP)   = MP(x.λ - y.λ)
Base.:(/)(x::MP,   y::Real) = MP(x.λ - b)
Base.:(/)(x::Real, y::MP)   = MP(a   - y.λ)

# ==============================================================================

"""
    -(x::MP, y::MP)

Minus operator (not used for max-plus but for min-plus).
Return the difference between `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) - MP(2.0)
MP{Float64}(-1.0)
```
"""
Base.:(-)(x::MP,   y::MP)   = MP(x.λ - y.λ)
Base.:(-)(x::MP,   y::Real) = MP(x.λ - y)
Base.:(-)(x::Real, y::MP)   = MP(x   - y.λ)

# ==============================================================================

"""
    min(x::MP, y::MP)

Return the minimun of `x` and `y`.

# Examples
```julia-repl
julia> min(MP(1), -3)
MP{Int64}(-3)

julia> min(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 Array{MP{Int64},2}:
 4  1
 6  1
```
"""
Base.:min(x::MP,   y::MP)   = MP(min(x.λ, y.λ))
Base.:min(x::MP,   y::Real) = MP(min(x.λ, y))
Base.:min(x::Real, y::MP)   = MP(min(x,   y.λ))
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

# ==============================================================================

Base.:(==)(x::MP,   y::MP)   = (x.λ == y.λ)
Base.:(==)(x::MP,   y::Real) = (x.λ == y)
Base.:(==)(x::Real, y::MP)   = (x   == y.λ)

Base.:(<=)(x::MP,   y::MP)   = (x.λ <= y.λ)
Base.:(<=)(x::MP,   y::Real) = (x.λ <= y)
Base.:(<=)(x::Real, y::MP)   = (x   <= y.λ)

Base.:(<)(x::MP,   y::MP)   = (x.λ < y.λ)
Base.:(<)(x::MP,   y::Real) = (x.λ < y)
Base.:(<)(x::Real, y::MP)   = (x   < y.λ)

Base.isless(x::MP,   y::MP)   = x.λ < y.λ
Base.isless(x::MP,   y::Real) = x.λ < y
Base.isless(x::Real, y::MP)   = x < y.λ

# ==============================================================================

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

"""
    mpfzeros(::Type{T}, n::Int64)

Construct a dense max-plus zero m-by-n matrix.

# Examples
```julia-repl
julia> mpfzeros(Float64, 2,2)
2×2 Array{MP{Float64},2}:
 -Inf  -Inf
 -Inf  -Inf
```
"""
mpfzeros(::Type{T}, m::Int64, n::Int64) where T = full(mpzeros(T, m, n))

# ==============================================================================

"""
    mpzeros(x::UnitRange{T})

Construct a max-plus zero vector.
# Examples
```julia-repl
julia> mpzeros(1:3)
3-element Array{MP{Int64},1}:
 -Inf
 -Inf
 -Inf
```
"""
mpzeros(x::UnitRange{T}) where {T} = fill!(MP(similar(x)), mpzero(T))

# ==============================================================================

"""
    mpones(::Type{T}, n::Int64)

Construct a max-plus one n-by-n matrix.

# Examples
```julia-repl
julia> mpones(Float64, 2)
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
julia> mpones(Float64, 2,2)
2×2 Array{MP{Float64},2}:
 0.0  0.0
 0.0  0.0
```
"""
mpones(::Type{T}, m::Int64, n::Int64) where T = ones(MP{T}, m, n)

# ==============================================================================

"""
    mpones(x::UnitRange{T})

Construct a max-plus one vector.
# Examples
```julia-repl
julia> mpones(1:3)
3-element Array{MP{Int64},1}:
 0
 0
 0
```
"""
mpones(x::UnitRange{T}) where {T} = fill!(MP(similar(x)), mpone(T))

# ==============================================================================

"""
    mpsparse(A::Array{T}; preserve=false)

Transform a dense matrix to a max-plus sparse matrix.

# Examples
```julia-repl
julia> mpsparse([4 0; 7 -Inf])
2×2 SparseMatrixCSC{MP{Float64},Int64} with 3 stored entries:
  [1, 1]  =  4.0
  [2, 1]  =  7.0
  [1, 2]  =  0.0
```
"""
mpsparse(A::Array{T}) where T = mpsparse(MP(A))
# ; preserve=false) where T = mpsparse(MP(A), preserve)

"""
    mpsparse(A::ArrMP{T}; preserve=false)

Transform a dense matrix to a max-plus sparse matrix. By default
Remove -Inf values. To keep them set param preserve to true.

# Examples
```julia-repl
julia> mpsparse([4 0; 7 -Inf])
2×2 SparseMatrixCSC{MP{Float64},Int64} with 3 stored entries:
  [1, 1]  =  4.0
  [2, 1]  =  7.0
  [1, 2]  =  0.0

mpsparse([4 0; 7 -Inf], preserve=true)
2×2 SparseMatrixCSC{MP{Float64},Int64} with 3 stored entries:
  [1, 1]  =  4.0
  [2, 1]  =  7.0
  [1, 2]  =  0.0
  [2, 2]  = -Inf
```
"""
function mpsparse(M::ArrMP{T}) where T #; preserve=false) where T
    #if preserve
    #    sparse(M)
    #else
        R = Array{Int64, 1}[];
        C = Array{Int64, 1}[];
        V = Array{MP{T}, 1}[];

        for j = 1:size(M, 2), i=1:size(M, 1)
            M[i,j] != zero(MP{T}) && (R = [R; i]; C = [C; j]; V = [V; M[i,j]])
        end

        sparse(R, C, V)
    #end
end

# ==============================================================================

"""
    array(::ArrMP{T})

Convert a max-plus array to a standard type array.

# Examples
```julia-repl
julia> array([MP(1.0) 2.0; 3.0 4.0])
2×2 Array{Float64,2}:
 1.0  2.0
 3.0  4.0
```
"""
array(A::ArrMP{T}) where T = map(x -> x.λ, A)

# ==============================================================================

"""
    array(::SpaMP{T}})

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
array(S::SpaMP{T,U}) where {T, U} = map(x -> x.λ, S)

# ==============================================================================

SparseArrays.sparse(S::SpaMP{T,U}) where {T, U} = map(x -> x.λ, S)

# ==============================================================================

"""
    full(::SpaMP{T}})

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
full(S::SpaMP{T,U}) where {T, U} = Matrix(S)

# ==============================================================================

"""
    dense(::SpaMP{T}})

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
dense(S::SpaMP{T,U}) where {T, U} = Matrix(S)

# ==============================================================================

"""
    plustimes(x::MP{T})

Convert max-plus number to the standard type.
Note: the function name comes from ScicosLab.

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
plustimes(n::MP{T})      where T      = n.λ
plustimes(A::ArrMP{T})   where T      = array(A)
plustimes(S::SpaMP{T,U}) where {T, U} = array(S)

# ==============================================================================
# Conversion to min-plus algebra. Convert +Inf and -Inf

"""
    minplus(A::Array)

Conversion to min-plus algebra. Convert +Inf and -Inf

# Examples
```julia-repl
julia> A = MP([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
minplus(A);
```
"""
minplus(n::MP{T})      where T      = map(x -> (x.λ == mp0) ? mptop : ((x.λ == mptop) ? mp0 : x), n)
minplus(A::ArrMP{T})   where T      = map(x -> minplus(x), A)
minplus(S::SpaMP{T,U}) where {T, U} = map(x -> minplus(x), S)

# ==============================================================================

"""
    mptrace(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> mptrace([MP(1) 2; 3 4])
4
```
"""
mptrace(A::ArrMP{T})    where T = isempty(A) ? mp0 : sum(diag(A))
mptrace(A::Array{T})    where T = isempty(A) ? mp0 : sum(MP(diag(A)))
mptrace(S::SpaMP{T,U})  where {T, U} = isempty(S) ? mp0 : sum(diag(S))
mptrace(S::Sparse{T,U}) where {T, U} = sum(MP(diag(S)))

# ==============================================================================

"""
    mpnorm(A)

Compute the norm of the full or sparce matrix A.
Return the largest entry minus smallest entry of A.

# Examples
```
julia-repl
julia> A = MP([1 20 2;30 400 4;4 50 10]);
S = mpsparse(A);
[mpnorm(A) mpnorm(S)]
1×2 Array{MP{Int64},2}:
 399  399
```
"""
mpnorm(A::Array{T}) where T = A[argmax(A)].λ - A[argmin(A)].λ
mpnorm(S::SpaMP{T}) where T = S[argmax(S)].λ - S[argmin(S)].λ

# ==============================================================================
function hstar(B::Array{T}) where T
    n = size(B, 1)
    if n != size(B, 2) error("Matrix shall be squared") end
    C = B
    for k in 1:n
        t = (C[k,k] <= zero(T)) ? zero(T) : typemax(T);
        for j in 1:n, i in 1:n
            C[i,j] = max(C[i,j], C[i,k] + C[k,j] + t)
        end
    end

    for k in 1:n
        C[k,k] = max(C[k,k], zero(T));
    end
    C
end

"""
    mpstar(A::Array{T})

TODO
"""
mpstar(x::MP{T})    where T = mpstar(plustimes(x))
mpstar(x::T)        where T = MP(hstar([x]))[1,1]
mpstar(A::Array{T}) where T = MP(hstar(A))
mpstar(A::ArrMP{T}) where T = MP(hstar(plustimes(A)))
#mpstar(S::Sparse{T,U}) where {T, U} = MP(hstar(S))
#mpstar(S::SpaMP{T,U}) where {T, U} = MP(hstar(plustimes(S)))

# ==============================================================================
# TODO insertion out of bounds => Sparse:  d=MP([1.0 2; 3 4]); d[5,5] = MP(6.0)

# ==============================================================================
include("juliabugs.jl")
include("howard.jl")
include("syslin.jl")
include("flowshop.jl")

end # MaxPlus module
