# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# In the way of the ScicosLab Max-Plus toolbox.
# ==============================================================================

module MaxPlus

using
    LinearAlgebra, SparseArrays, PrettyTables, Printf

# ==============================================================================
# Max-Plus core

export
    MP, SpaMP, SpvMP, ArrMP, VecMP, mpzero, mpone, mp0, mp1, mptop, ϵ, mpe, mpI,
    mpeye, mpzeros, mpones, mpsparse, full, dense, array, plustimes, minplus,
    mpsparse_map, mptrace, mpnorm, mpstar, mpplus, howard, mp_change_display,
    mpshow, LaTeX

# ==============================================================================
# Max-Plus Linear system

export
    MPSysLin, mpsyslin, mpsimul, mpexplicit

# ==============================================================================
# Max-Plus flowhop

export
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

# ==============================================================================
# Max-Plus type

"""
    MP{T}

Max-Plus immutable structure. Promote a number of type T (ideally a Float64) to
a Max-Plus number.

MP(3.0) is the equivalent to ScicosLab code: `#(3.0)`

# Examples
```julia-repl
julia> a = MP(3.5)
MP{Float64}(3.5)

julia> b = MP(3)
MP{Int64}(3)

julia> typeof(a), typeof(b)
(MP{Float64}, MP{Int64})
```
"""
struct MP{T <: Real} <: Real
    λ::T
end

# ==============================================================================
# Type alias but with shorter number of characters

const Sparse{T,U} = SparseMatrixCSC{T,U}      # Sparse matrix (classic algebra)
const SpaMP{T,U}  = SparseMatrixCSC{MP{T},U}  # Sparse matrix (max-plus algebra)
const SpvMP{T,U}  = SparseVector{MP{T},U}     # Sparse vector (max-plus algebra)
const ArrMP{T,N}  = Array{MP{T},N}            # Dense matrix (max-plus algebra)
const VecMP{T,N}  = Vector{MP{T}}             # Dense vector (max-plus algebra)

# ==============================================================================
# Julia type promotions and conversions to/from Max-Plus number

Base.promote_rule(::Type{MP{T}}, ::Type{U}) where {T, U} = MP{T}
Base.convert(::MP{T}, x::Number)            where T = MP{T}(x)

# ==============================================================================
# Copy constructor

MP(x::MP) = MP(x.λ)

# ==============================================================================
# Constructor from non Max-Plus dense matrix

"""
    MP(A::Array)

Convert a dense array from classic algebra to a Max-Plus dense array.

# Examples
```julia-repl
julia> MP([1.0 -Inf; 0.0 4])
2×2 Matrix{MP{Float64}}:
  1.0   -Inf
  0.0    4.0
```
"""
MP(A::Array) = map(MP, A)

# ==============================================================================
# Constructor from non Max-Plus sparse matrix

"""
    MP(A::SparseArray; keepzeros=true)

Convert a sparse array from classic algebra to a Max-Plus sparse array. By
default, explicit Max-Plus zeros (`ϵ`, `mp0`, `MP(-Inf)`) are kept except if
the parameter `keepzeros` is set to `false`.

# Examples
```julia-repl
julia> using SparseArrays

julia> A = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
3×3 SparseMatrixCSC{MP{Float64}, Int64} with 3 stored entries:
  .   .   .
  .   2   .
  .   .   0

julia> A.nzval
3-element Vector{MP{Float64}}:
  -Inf
   2.0
   0.0

julia> B = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]), keepzeros=false)
3×3 SparseMatrixCSC{MP{Float64}, Int64} with 2 stored entries:
  .   .   .
  .   2   .
  .   .   0

julia> B.nzval
2-element Vector{MP{Float64}}:`
  2.0
  0.0

julia> A == B
true
```
"""
function MP(S::Sparse{T,U}; keepzeros::Bool=true) where {T, U}
    if keepzeros
        convert(SpaMP{T,U}, S)
    else
        M = spzeros(MP{T}, size(S,1), size(S,2))
        for c = 1:size(S, 2), r = nzrange(S, c)
            if (MP(S[r,c]) != mpzero(T))
                M[r,c] = convert(MP{T}, S[r,c])
            end
        end
        M
    end
end

# TODO keepzeros=true to remove mp0
MP(S::Sparse{MP{T},U}) where {T, U} = S

# ==============================================================================
# Constructor from non Max-Plus range

"""
    MP(x::UnitRange)

Create a Max-Plus dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1:3)
3-element Vector{MP{Int64}}:
 1
 2
 3
```
"""
MP(x::UnitRange{T}) where T = Vector{MP{T}}(x)

"""
    MP(x::StepRangeLen)

Create a Max-Plus dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1.0:0.5:3.0)
5-element Vector{MP{Float64}}:
  1.0
  1.5
  2.0
  2.5
  3.0
```
"""
MP(x::StepRangeLen{T, Base.TwicePrecision{T}, Base.TwicePrecision{T}}) where T = Vector{MP{T}}(x)

# ==============================================================================
# Algebra redefinition for Max-Plus: zero

"""
    zero(::MP{T})

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ϵ`.
"""
Base.zero(::Type{MP{T}}) where T = MP(typemin(T))
Base.zero(x::MP{T}) where T = zero(typeof(x))

"""
    mpzero(::Type{T})

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ϵ`.

# Examples
```julia-repl
julia> mpzero(Float64)
MP{Float64}(-Inf)

julia> mpzero(Int64)
MP{Int64}(-9223372036854775808)
```
"""
mpzero(::Type{T}) where T = MP(typemin(T))

"""
    mpzero()

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ϵ`.

# Examples
```julia-repl
julia> mpzero()
MP{Float64}(-Inf)
```
"""
mpzero() = MP(typemin(Float64))

# ==============================================================================
# Algebra redefinition for Max-Plus: one

"""
    one(::MP{T})

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
    the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.
"""
Base.one(::Type{MP{T}}) where T = MP(zero(T))
Base.one(x::MP{T}) where T = one(typeof(x))

"""
    one(::Type{T})

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.

# Examples
```julia-repl
julia> mpone(Float64)
MP{Float64}(0.0)

julia> mpone(Int64)
MP{Int64}(0)
```
"""
mpone(::Type{T}) where T = MP(zero(T))

"""
    mpone()

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.

# Examples
```julia-repl
julia> mpone()
MP{Float64}(0.0)
```
"""
mpone() = MP(zero(Float64))

# ==============================================================================
# Max-Plus core plus operator

"""
    +(x::MP, y::MP)

Max-Plus operator ⨁. Return the maximum of `x` and `y` as Max-Plus type. At
least one parameter shall be a Max-Plus number (its conversion to Max-Plus is
automatic).

# Examples
```julia-repl
julia> MP(1.0) + MP(3.0)
MP{Float64}(3.0)

julia> MP(1.0) + 3
MP{Float64}(3.0)

julia> 1 + MP(3.0)
MP{Float64}(3.0)

julia> MP(3.0) + -Inf
MP{Float64}(3.0)
```
"""
Base.:(+)(x::MP,   y::MP)   = MP(max(x.λ, y.λ))
Base.:(+)(x::MP,   y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP)   = MP(max(x,   y.λ))

# ==============================================================================
# Max-Plus core times operator

"""
    *(x::MP, y::MP)

Max-Plus operator ⨂. Return the sum of `x` and `y` as Max-Plus type. At least
one parameter shall be a Max-Plus number (its conversion to Max-Plus is
automatic).

# Examples
```julia-repl
julia> MP(1.0) * MP(3.0)
MP{Float64}(4.0)

julia> MP(1.0) * 3
MP{Float64}(4.0)

julia> 1 * MP(3.0)
MP{Float64}(4.0)

julia> MP(1.0) * -Inf
MP{Float64}(-Inf)
```
"""
Base.:(*)(x::MP,   y::MP)   = MP(x.λ + y.λ)
Base.:(*)(x::MP,   y::Real) = MP(x.λ + y)
Base.:(*)(x::Real, y::MP)   = MP(x   + y.λ)

# ==============================================================================
# Max-Plus and Min-Plus constants

"""
    mp0

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mpzero`.

Equivalent to ScicosLab code: `%0` sugar notation for `#(-%inf)`

# Examples
```julia-repl
julia> mp0
MP{Float64}(-Inf)

julia> mp0 * 5
MP{Float64}(-Inf)

julia> mp0 + 5
MP{Float64}(5.0)
```
"""
const global mp0 = mpzero()

"""
    ϵ (\\epsilon)

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mpzero`.

Equivalent to ScicosLab code: `%0` sugar notation for `#(-%inf)`

# Examples
```julia-repl
julia> ϵ
MP{Float64}(-Inf)

julia> ϵ * 5
MP{Float64}(-Inf)

julia> ϵ + 5
MP{Float64}(5.0)
```
"""
const global ϵ = mp0

"""
    mp1

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mpone`.

Equivalent to ScicosLab code: `%1` sugar notation for `#(1)`

# Examples
```julia-repl
julia> mp1
MP{Float64}(0.0)

julia> mp1 * 5
MP{Float64}(5.0)

julia> mp1 + 5
MP{Float64}(5.0)
```
"""
const global mp1 = mpone()

"""
    mpe

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mpone`.

Equivalent to ScicosLab code: `%1` sugar notation for `#(1)`

# Examples
```julia-repl
julia> mpe
MP{Float64}(0.0)

julia> mpe * 5
MP{Float64}(5.0)

julia> mpe + 5
MP{Float64}(5.0)
```
"""
const global mpe = mp1

"""
    mptop

Create the constant Min-Plus one (equals to `+∞`, infinity) in classic algebra which is
the neutral for operators ⨁ and ⨂.

Equivalent to ScicosLab code: `%top = #(%inf)`

# Examples
```julia-repl
julia> mptop
MP{Float64}(Inf)

julia> mptop * 5
MP{Float64}(Inf)

julia> mptop + 5
MP{Float64}(Inf)
```
"""
const global mptop = MP{Float64}(Inf)

# ==============================================================================
# Max-Plus sparse map function

"""
    mpsparse_map(f, M::SparseMatrixCSC{MP{T},U})

Map a function ot each element of the Max-Plus sparse matrix.

# Examples
```julia-repl
julia> mpsparse_map(x -> x.λ, mpsparse([1.0 0; 0 1.0]))
2×2 SparseMatrixCSC{Float64, Int64} with 4 stored entries:
 1.0  0.0
 0.0  1.0
```
"""
function mpsparse_map(f, M::SparseMatrixCSC{MP{T},U}) where {T, U}
    SparseMatrixCSC(M.m, M.n, M.colptr, M.rowval, map(f, M.nzval))
end

# ==============================================================================
# Algebra conversion: Max-Plus to Min-Plus or Max-Plus to classic algebra

"""
    plustimes(x::MP{T})

Convert a Max-Plus number to a number in standard algebra. An alternative way
could be `x.λ`.

# Examples
```julia-repl
julia> plustimes(MP(1.0))
1.0

julia> typeof(plustimes(MP(1)))
Int64
```
"""
plustimes(n::MP{T}) where T = n.λ

"""
    plustimes(A::ArrMP{T})

Convert a Max-Plus dense matrix to a dense matrix in standard algebra.

# Examples
```julia-repl
julia> A=[MP(1.0) 2.0; ϵ mpe]
2×2 Matrix{MP{Float64}}:
   1.0   2.0
  -Inf   0.0

julia> plustimes(A)
2×2 Matrix{Float64}:
   1.0  2.0
  -Inf  0.0
```
"""
plustimes(A::ArrMP{T}) where T = map(x -> x.λ, A)

"""
    plustimes(A::SpaMP{T})

Convert a Max-Plus sparse matrix to an sparse matrix in standard algebra.

# Examples
```julia-repl
julia> S = mpsparse(mpeye(2,2))
2×2 SparseMatrixCSC{MP{Float64}, Int64} with 2 stored entries:
 0.0   ⋅
  ⋅   0.0

julia> findnz(S)
([1, 2], [1, 2], MP{Float64}[0.0, 0.0])

julia> plustimes(S)
2×2 SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 0.0   ⋅
  ⋅   0.0

julia> findnz(plustimes(S))
([1, 2], [1, 2], [0.0, 0.0])
```
"""
function plustimes(S::SpaMP{T,U}) where {T, U}
    SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(x -> x.λ, S.nzval))
end

"""
    minplus(x::MP{T})

Conversion a Max-Plus number to Min-Plus number. This function convert `+∞` and
`-∞` to their opposite sign.

# Examples
```julia-repl
julia> minplus(mptop), minplus(mp0), minplus(MP(4.5)), minplus(MP(-4.5))
(-Inf, Inf, 4.5, -4.5)

julia> typeof(minplus(mptop))
MP{Float64}
```
"""
minplus(x::MP{T}) where T = (x == mp0) ? mptop : ((x == mptop) ? mp0 : x)

"""
    minplus(A::ArrMP)

Conversion a Max-Plus dense matrix to Min-Plus dense matrix. This function
convert `+∞` and `-∞` to their opposite sign.

# Examples
```julia-repl
julia> A = MP([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
3×4 Matrix{MP{Float64}}:
   0.0   3.0   Inf    1.0
   1.0   2.0   2.0   -Inf
  -Inf   Inf   1.0    0.0

julia> minplus(A)
3×4 Matrix{MP{Float64}}:
  0.0    3.0   -Inf   1.0
  1.0    2.0    2.0   Inf
  Inf   -Inf    1.0   0.0
```
"""
minplus(A::ArrMP{T})   where T = map(x -> minplus(x), A)

"""
    minplus(A::SpaMP)

Conversion a Max-Plus dense matrix to Min-Plus dense matrix. This function
convert `+∞` and `-∞` to their opposite sign.

# Examples
```julia-repl
julia> S = mpsparse([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
3×4 SparseMatrixCSC{MP{Float64}, Int64} with 10 stored entries:
  0     3   Inf   1
  1     2     2   .
  .   Inf     1   0

julia> minplus(S)
3×4 SparseMatrixCSC{MP{Float64}, Int64} with 10 stored entries:
    0   3   .     1
    1   2   2   Inf
  Inf   .   1     0
```
"""
minplus(S::SpaMP{T,U}) where {T, U} = dropzeros(map(x -> minplus(x), S))

"""
    sparse(S::SpaMP{T,U})

Convert a dense Max-Plus array to a sparse Max-Plus array.

# Examples
```julia-repl
julia> S = mpsparse([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
3×4 SparseMatrixCSC{MP{Float64}, Int64} with 10 stored entries:
  0     3   Inf   1
  1     2     2   .
  .   Inf     1   0

julia> minplus(S)
3×4 SparseMatrixCSC{MP{Float64}, Int64} with 12 stored entries:
    0   3   .     1
    1   2   2   Inf
  Inf   .   1     0
```
"""
SparseArrays.sparse(S::SpaMP{T,U}) where {T, U} = map(x -> x.λ, S)

"""
    full(::SpaMP{T}})

Convert a sparse Max-Plus array to a dense Max-Plus array.
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

"""
    dense(::SpaMP{T}})

Convert a sparse Max-Plus array to a dense Max-Plus array.
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

"""
    array(::SpaMP{T}})

Convert a sparse Max-Plus array to a dense non Max-Plus array.

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
# Display

"""
    mpstyle

Memory for saving the style of display for neutral and absorbing elements. By
default the display of ScicosLab will be used.
"""
global mpstyle = 1; # Style ScicosLab

"""
    mp_change_display(style::Int)

Change the style of behavior of functions `Base.show()`:
- `-Inf` are displayed either with `ε` (style 2 or 3) or `.` symbols (style 1).
- `0` are displayed either with `e `(style 3) or '0' symbols (style 1 or 2).
- else: `-Inf` and `0` are displayed in Julia default sytle (style 0).

If this function is not called, by default the ScicosLab style will be used
(style 1).
"""
mp_change_display(style::Int) = (global mpstyle = min(max(style, 0), 4))

# Base function for display a Max-Plus number.
function mpshow(io::IO, x::MP{T}) where T
    if (mpstyle == 0)
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

function mpshow(io::IO, A::ArrMP{T}) where T
    if (size(A,2) == 1)
        print(io, size(A,1), "-element Vector{MP{$T}}:\n")
    else
        print(io, size(A,1), '×', size(A,2), " Matrix{MP{$T}}:\n")
    end
    pretty_table(io, A, tf = tf_borderless, noheader = true)
end

function mpshow(io::IO, S::SpaMP{T,U}) where {T,U}
    print(io, size(S,1), '×', size(S,2), " SparseMatrixCSC{MP{$T}, $U} with ",
          size(S.nzval,1), " stored entries:\n")
    old_mpstyle = mpstyle # Force showing empty elements as .
    mp_change_display(1)
    pretty_table(io, S, tf = tf_borderless, noheader = true)
    mp_change_display(old_mpstyle)
end

"""
    LaTeX(io::IO, A::Array{MP})

Base function for convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
neutral and absorbing elements depends on mp_change_display(style).

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
# FIXME: not a column-major traversal
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

# Base function for displaying a Max-Plus sparse matrix.
LaTeX(io::IO, S::SpaMP{T}) where {T,U} = LaTeX(io, full(S))

"""
    show(io::IO, x::MP)

Display a Max-Plus number depending on the currently set style:

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
# Called by pretty_table() when REPL shows a MP matrix.
# julia> [MP(1) MP(2); MP(3) MP(4)]
Base.show(io::IO, x::MP{T}) where T = mpshow(io, x)

# Called by REPL when showing a MP scalar.
# julia> MP(1)
function Base.show(io::IO, ::MIME"text/plain", x::MP{T}) where T
    print(io, "MP{$T}(")
    mpshow(io, x)
    print(io, ")")
end

# Called by the REPL through the display() method. This function fixes
# misaliged columns made by the default show() Julia. We use the package
# PrettyTables.
# julia> [MP(1) MP(2); MP(3) MP(4)]
Base.show(io::IO, ::MIME"text/plain", A::ArrMP{T}) where T = mpshow(io, A)

Base.show(io::IO, ::MIME"text/plain", S::SpaMP{T,U}) where {T,U} = mpshow(io, S)

# Convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", A::ArrMP{T}) where T = LaTeX(io, A)

# Convert a Max-Plus sparse matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", x::SpaMP{T}) where T = LaTeX(io, A)

# ==============================================================================
# These functions fix Julia bugs

@inline Base.literal_pow(::typeof(^), x::MP, ::Val{0}) = one(x)
@inline Base.literal_pow(::typeof(^), x::MP, ::Val{p}) where {p} = MP(x.λ * p)
@inline Base.:(^)(x::MP, y::Number) = MP(x.λ * y)
@inline Base.literal_pow(::typeof(^), A::ArrMP{T}, ::Val{p}) where {T, p} = A^p
@inline Base.abs2(x::MP) = x.λ + x.λ

# ==============================================================================
# Hardly used operators

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

"""
    -(x::MP, y::MP)

Minus operator (not used for Max-Plus but for Min-Plus).
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

"""
    min(x::MP, y::MP)

Return the minimun of Max-Plus numbers `x` and `y`.

# Examples
```julia-repl
julia> min(MP(1), 3)
MP{Int64}:
  1

julia> min(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 Matrix{MP{Int64}}:
  4   1
  6   1
```
"""
Base.:min(x::MP,   y::MP)   = MP(min(x.λ, y.λ))
Base.:min(x::MP,   y::Real) = MP(min(x.λ, y))
Base.:min(x::Real, y::MP)   = MP(min(x,   y.λ))
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

"""
    max(x::MP, y::MP)

Return the max of `x` and `y`.

# Examples
```julia-repl
julia> max(MP(1), 3)
MP{Int64}:
  3

julia> max(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 Matrix{MP{Int64}}:
  10   5
  10   5
```
"""
Base.:max(x::MP,   y::MP)   = MP(max(x.λ, y.λ))
Base.:max(x::MP,   y::Real) = MP(max(x.λ, y))
Base.:max(x::Real, y::MP)   = MP(max(x,   y.λ))
Base.:max(A::ArrMP, B::ArrMP) = map(Base.:max, A, B)
Base.:max(A::SpaMP, B::SpaMP) = map(Base.:max, A, B)

# ==============================================================================
# Comparator

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
# Identity matrix

"""
    mpeye(::Type{T}, n::Int64)

Construct a Max-Plus identity dense n-by-n matrix.

# Examples
```julia-repl
julia> mpeye(Float64, 2)
2×2 Matrix{MP{Float64}}:
   0.0   -Inf
  -Inf    0.0
```
"""
mpeye(::Type{T}, n::Int64) where T = Matrix{MP{T}}(mpI, n, n)

"""
    mpeye(n::Int64)

Construct a Max-Plus identity dense n-by-n matrix.

# Examples
```julia-repl
julia> mpeye(2)
2×2 Matrix{MP{Float64}}:
   0.0   -Inf
  -Inf    0.0
```
"""
mpeye(n::Int64) = Matrix{MP{Float64}}(mpI, n, n)

"""
    mpeye(::Type{T}, m::Int64, n::Int64)

Construct a Max-Plus identity dense m-by-n matrix.

# Examples
```julia-repl
julia> mpeye(Float64, 2, 3)
2×3 Matrix{MP{Float64}}:
   0.0   -Inf   -Inf
  -Inf    0.0   -Inf
```
"""
mpeye(::Type{T}, m::Int64, n::Int64) where T = Matrix{MP{T}}(mpI, m, n)

"""
    mpeye(m::Int64, n::Int64)

Construct a Max-Plus identity dense m-by-n matrix.

# Examples
```julia-repl
julia> mpeye(2, 3)
2×3 Matrix{MP{Float64}}:
   0.0   -Inf   -Inf
  -Inf    0.0   -Inf
```
"""
mpeye(m::Int64, n::Int64) = Matrix{MP{Float64}}(mpI, m, n)

"""
    mpeye(A::Array{T})

Construct a Max-Plus identity dense matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> mpeye(A)
2×2 Matrix{MP{Float64}}:
  0.0  -Inf
 -Inf   0.0

julia> mpeye(MP(A))
2×2 Matrix{MP{Float64}}:
  0.0  -Inf
 -Inf   0.0
```
"""
mpeye(A::Array{T}) where T = Matrix{MP{T}}(mpI, size(A,1), size(A,2))
mpeye(A::ArrMP{T}) where T = Matrix{MP{T}}(mpI, size(A,1), size(A,2))

# ==============================================================================
# Zero matrix

"""
    mpzeros(::Type{T}, n::Int64)

Construct a Max-Plus zero n-by-m sparse matrix.

# Examples
```julia-repl
julia> mpzeros(Float64, 2)
2-element SparseVector{MP{Float64},Int64} with 0 stored entries
```
"""
mpzeros(::Type{T}, n::Int64) where T = spzeros(MP{T}, n)

"""
    mpzeros(n::Int64)

Construct a Max-Plus zero n-by-m sparse matrix.

# Examples
```julia-repl
julia> mpzeros(2)
2-element SparseVector{MP{Float64},Int64} with 0 stored entries
```
"""
mpzeros(n::Int64) = spzeros(MP{Float64}, n)

"""
    mpzeros(::Type{T}, n::Int64)

Construct a sparse Max-Plus zero m-by-n matrix.

# Examples
```julia-repl
julia> mpzeros(Float64, 2,5)
2×5 SparseMatrixCSC{MP{Float64},Int64} with 0 stored entries
```
"""
mpzeros(::Type{T}, m::Int64, n::Int64) where T = spzeros(MP{T}, m, n)

"""
    mpzeros(n::Int64)

Construct a sparse Max-Plus zero m-by-n matrix.

# Examples
```julia-repl
julia> mpzeros(2,5)
2×5 SparseMatrixCSC{MP{Float64},Int64} with 0 stored entries
```
"""
mpzeros(m::Int64, n::Int64) = spzeros(MP{Float64}, m, n)

"""
    mpzeros(A::Array{T})

Construct a sparse  Max-Plus zero matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> mpzeros(A)
2×2 SparseMatrixCSC{MP{Float64}, Int64} with 0 stored entries:
  -Inf   -Inf
  -Inf   -Inf

julia> mpzeros(MP(A))
2×2 SparseMatrixCSC{MP{Float64}, Int64} with 0 stored entries:
  -Inf   -Inf
  -Inf   -Inf
```
"""
mpzeros(A::Array{T}) where T = spzeros(MP{T}, size(A,1), size(A,2))
mpzeros(A::ArrMP{T}) where T = spzeros(MP{T}, size(A,1), size(A,2))

# ==============================================================================
# One matrix

"""
    mpones(::Type{T}, n::Int64)

Construct a Max-Plus one n-by-1 matrix.

# Examples
```julia-repl
julia> mpones(Float64, 2)
2×1 Matrix{MP{Float64}}:
  0.0
  0.0
```
"""
mpones(::Type{T}, n::Int64) where T = ones(MP{T}, n)

"""
    mpones(n::Int64)

Construct a Max-Plus one n-by-1 matrix.

# Examples
```julia-repl
julia> mpones(2)
2×1 Matrix{MP{Float64}}:
  0.0
  0.0
```
"""
mpones(n::Int64) where T = ones(MP{Float64}, n)

"""
    mpones(::Type{T}, m::Int64, n::Int64)

Construct a Max-Plus one m-by-n matrix.

# Examples
```julia-repl
julia> mpones(Float64, 3,2)
3×2 Matrix{MP{Float64}}:
  0.0   0.0
  0.0   0.0
  0.0   0.0
```
"""
mpones(::Type{T}, m::Int64, n::Int64) where T = ones(MP{T}, m, n)

"""
    mpones(m::Int64, n::Int64)

Construct a Max-Plus one m-by-n matrix.

# Examples
```julia-repl
julia> mpones(3,2)
3×2 Matrix{MP{Float64}}:
  0.0   0.0
  0.0   0.0
  0.0   0.0
```
"""
mpones(m::Int64, n::Int64) where T = ones(MP{Float64}, m, n)

"""
    mpones(A::Array{T})

Construct a sparse Max-Plus ones matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> mpones(A)
2×2 Matrix{MP{Float64}}:
  0.0   0.0
  0.0   0.0
```
"""
mpones(A::Array{T}) where T = mpones(T, size(A,1), size(A,2))
mpones(A::ArrMP{T}) where T = mpones(T, size(A,1), size(A,2))

# ==============================================================================
# Dense and Sparse matrices

"""
    mpsparse(A::Array{T}; keepzeros::Bool=false)

Transform a dense matrix from classic algebra to a Max-Plus sparse matrix.
Max-Plus zeros (`ϵ`, `mp0`, `MP(-Inf)`) are removed and classic algebra
zeros are removed if the argument `keepzeros` is set to `true`.

# Arguments
- keepzeros: if true then 0 values from classic algebra are considered as values
and are not removed.
- keepzeros: if false then 0 values from classic algebra are removed.

# Examples
```julia-repl
julia> S = mpsparse([-Inf 0; 0 -Inf])
2×2 SparseMatrixCSC{MP{Float64}, Int64} with 2 stored entries:
  .   0
  0   .

julia> findnz(S)
([2, 1], [1, 2], MP{Float64}[0.0, 0.0])

julia> S = mpsparse([-Inf 0; 0 -Inf], keepzeros=false)
2×2 SparseMatrixCSC{MP{Float64}, Int64} with 0 stored entries:
  .   .
  .   .

julia> findnz(S)
(Int64[], Int64[], MP{Float64}[])
```
"""
mpsparse(A::Array{T}; keepzeros::Bool=true) where T =
    keepzeros ? sparse(MP(A)) : MP(sparse(A), keepzeros=false)

"""
    mpsparse(A::ArrMP{T})

Transform a dense Max-Plus matrix to a sparse Max-Plus matrix.
Zero values from Max-Plus algebra are removed (therefore -Inf
from classic algebra are removed).

# Examples
```julia-repl
julia> S = mpsparse(MP([4 0; 7 -Inf]))
2×2 SparseMatrixCSC{MP{Float64}, Int64} with 3 stored entries:
  4   0
  7   .

julia> findnz(S)
([1, 2, 1], [1, 1, 2], MP{Float64}[4, 7, 0])
```
"""
mpsparse(A::ArrMP{T}) where T = sparse(A)

# ==============================================================================
# Max-Plus matrices operations

"""
    mptrace(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> mptrace([MP(1) 2; 3 4])
MP{Int64}(4)
```
"""
mptrace(A::ArrMP{T})    where T = isempty(A) ? mp0 : sum(diag(A))
mptrace(A::Array{T})    where T = isempty(A) ? mp0 : sum(MP(diag(A)))
mptrace(S::SpaMP{T,U})  where {T, U} = isempty(S) ? mp0 : sum(diag(S))
mptrace(S::Sparse{T,U}) where {T, U} = sum(MP(diag(S)))

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
# Max-Plus star FIXME not working with int

"""
    B = mpstar(A::ArrMP{T})

Solve `x = Ax + I` in the Max-Plus algebra. When there is no circuits with
positive weight in G(A) (the incidence graph of A) `B = I + A + ... + A^(n-1)`
where n denotes the order of the square matrix A.

See also mpplus.

# Arguments
- A : Max-Plus full square matrix.
- B : Id + A + A^2 + ...

# Examples
```julia-repl
julia> mpstar(MP(1.0))
MP{Float64}:
  Inf

julia> mpstar(MP(-1.0))
MP{Float64}:
  0

julia> mpstar(MP([1.0 2; 3 4]))
2×2 Matrix{MP{Float64}}:
  Inf   Inf
  Inf   Inf

julia> A = mpstar(MP([-3.0 -2;-1 0])); B = mpstar(A)
2×2 Matrix{MP{Float64}}:
   0   -2
  -1    0

julia> B == (mpeye(2,2) + A)
true

julia> B == (B + A * A)
true
```
"""
function mpstar(A::ArrMP{T}) where T
    n = size(A, 1)
    if n != size(A, 2)
        error("Matrix shall be squared")
    end
    C = A
    for k in 1:n
        t::T = (C[k,k].λ <= zero(T)) ? zero(T) : typemax(T);
        for j in 1:n, i in 1:n
            C[i,j] = MP{T}(max(C[i,j].λ, C[i,k].λ + C[k,j].λ + t))
        end
    end

    for k in 1:n
        C[k,k] = MP{T}(max(C[k,k].λ, zero(T)));
    end
    C
end

"""
    mpstar(x::MP{T})

Make x a 1x1 matrix then call mpstar(A::ArrMP{T}).
See mpstar(A::ArrMP{T}) for more information.
"""
mpstar(x::MP{T}) where T = mpstar([x])[1,1]

"""
    mpastarb(A::ArrMP{T}, b::ArrMP{T})

Max-Plus linear system solution.

Solve `x = Ax + b` in the Max-Plus algebra when there is no circuits with
positive weight in `G(A')` (the incidence graph of `A'`, that is it exists an
arc from `j` to `i` if `A_ij` is nonzero).

TODO It is much more efficient in time and memory than `mpstar(A) * b`.
"""
mpastarb(A::ArrMP{T}, b::ArrMP{T}) where T = mpstar(A) * b
# TODO: optimal code

"""
    B = mpplus(A::ArrMP{T})

Compute `A * A^* = A + A^2 + ...` of a maxplus matrix A.
See also mpstar.

# Arguments
- A : Max-Plus full square matrix.
- B : Id + A + A^2 + ...

# Examples
```julia-repl
julia> mpplus(MP(1.0))
MP{Float64}:
  Inf

julia> mpplus(MP(-1.0))
MP{Float64}:
  -1

julia> A = MP([-3.0 -2; -1 0]); B = mpplus(A)
2×2 Matrix{MP{Float64}}:
  -3   -2
  -1    0

julia> B == (A * mpstar(A))
true
```
"""
function mpplus(A::ArrMP{T}) where T
    n = size(A, 1)
    if n != size(A, 2)
        error("Matrix shall be squared")
    end
    C = A
    for k in 1:n
        t::T = (C[k,k].λ <= zero(T)) ? zero(T) : typemax(T);
        for j in 1:n, i in 1:n
            C[i,j] = MP{T}(max(C[i,j].λ, C[i,k].λ + C[k,j].λ + t))
        end
    end
    C
end

"""
    mpplus(x::MP{T})

Make x a 1x1 matrix then call mpplus(A::ArrMP{T}).
See mpplus(A::ArrMP{T}) for more information.
"""
mpplus(x::MP{T}) where T = mpplus([x])[1,1]

# ==============================================================================
# TODO insertion out of bounds => Sparse:  d=MP([1.0 2; 3 4]); d[5,5] = MP(6.0)

# ==============================================================================
include("fallbacks.jl")
include("howard.jl")
include("syslin.jl")
include("flowshop.jl")

end # MaxPlus module
