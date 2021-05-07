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
    MP, SpaMP, SpvMP, ArrMP, VecMP, mpzero, mpone, mp0, mp1, mptop, ε, mpe, mpI,
    mpeye, mpzeros, mpones, mpsparse, full, dense, array, plustimes, minplus,
    mpsparse_map, mptrace, mpnorm, mpstar, mpplus, howard, mp_change_display,
    mpshow, LaTeX

# ==============================================================================
# Max-Plus Linear system

export
    MPSysLin, mpsimul, mpexplicit

# ==============================================================================
# Max-Plus flowhop

export
    mpgraph, flowshop, LaTeXG, flowshop_graph, flowshop_simu

# ==============================================================================
# Max-Plus type

"""
    MP

Immutable Julia structure for Max-Plus scalar. Promote a number of type Float64
or Int64 to a number in the tropical semi-ring max,+ (ℝ ∪ {-∞}, ⊕, ⊙) where ℝ is
the domain of reals, ⊕ is the usual multiplication and ⊙ is the usual maximum.

`MP(3)` is the equivalent to ScicosLab code: `#(3)`

# Examples
```julia-repl
julia> a = MP(3.5)
Max-Plus 3.5)

julia> b = MP(3)
Max-Plus 3)

julia> typeof(a), typeof(b)
(MP, MP)
```
"""
struct MP <: Real
    λ::Float64
end

# ==============================================================================
# Type alias but with shorter number of characters

const Sparse{T,U} = SparseMatrixCSC{T,U}  # Sparse matrix (classic algebra)
const SpaMP{U}  = SparseMatrixCSC{MP,U}   # Sparse matrix (max-plus algebra)
const SpvMP{U}  = SparseVector{MP,U}      # Sparse vector (max-plus algebra)
const ArrMP{N}  = Array{MP,N}             # Dense matrix (max-plus algebra)
const VecMP{N}  = Vector{MP}              # Dense vector (max-plus algebra)

# ==============================================================================
# Julia type promotions and conversions to/from Max-Plus number

Base.promote_rule(::Type{MP}, ::Type{T}) where T = MP
Base.convert(::MP, x::Number) = MP(x)

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
2×2 Max-Plus dense matrix:
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
default, explicit Max-Plus zeros (`ε`, `mp0`, `MP(-Inf)`) are kept except if
the parameter `keepzeros` is set to `false`.

# Examples
```julia-repl
julia> using SparseArrays

julia> A = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
3×3 Max-Plus sparse matrix with 3 stored entries:
  .   .   .
  .   2   .
  .   .   0

julia> A.nzval
3-element Vector{max+}:
  -Inf
   2.0
   0.0

julia> B = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]), keepzeros=false)
3×3 Max-Plus sparse matrix with 2 stored entries:
  .   .   .
  .   2   .
  .   .   0

julia> B.nzval
2-element Vector{max+}:`
  2.0
  0.0

julia> A == B
true
```
"""
function MP(S::Sparse{T,U}; keepzeros::Bool=true) where {T,U}
    if keepzeros
        convert(SpaMP{U}, S)
    else
        M = spzeros(MP, size(S,1), size(S,2))
        for c = 1:size(S, 2), r = nzrange(S, c)
            if (MP(S[r,c]) != mpzero())
                M[r,c] = convert(MP, S[r,c])
            end
        end
        M
    end
end

# ==============================================================================
# Constructor from non Max-Plus range

"""
    MP(x::UnitRange)

Create a Max-Plus dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1:3)
3-element Vector{max+}:
 1.0
 2.0
 3.0
```
"""
MP(x::UnitRange) = Vector{MP}(x)

"""
    MP(x::StepRangeLen)

Create a Max-Plus dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1.0:0.5:3.0)
5-element Vector{max+}:
  1.0
  1.5
  2.0
  2.5
  3.0
```
"""
MP(x::StepRangeLen) = Vector{MP}(x)

# ==============================================================================
# Algebra redefinition for Max-Plus: zero

"""
    zero(::MP)

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.
"""
Base.zero(::Type{MP}) = MP(typemin(Float64))
Base.zero(x::MP) = zero(typeof(x))

"""
    mpzero()

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.

# Examples
```julia-repl
julia> mpzero()
Max-Plus -Inf
```
"""
mpzero() = MP(typemin(Float64))

# ==============================================================================
# Algebra redefinition for Max-Plus: one

"""
    one(::MP)

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
    the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.
"""
Base.one(::Type{MP}) = MP(zero(Float64))
Base.one(x::MP) = one(typeof(x))

"""
    mpone()

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.

# Examples
```julia-repl
julia> mpone()
MP(0.0)
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
MP(3.0)

julia> MP(1.0) + 3
MP(3.0)

julia> 1 + MP(3.0)
MP(3.0)

julia> MP(3.0) + -Inf
MP(3.0)
```
"""
Base.:(+)(x::MP, y::MP) = MP(max(x.λ, y.λ))
Base.:(+)(x::MP, y::Real) = MP(max(x.λ, y))
Base.:(+)(x::Real, y::MP) = MP(max(x, y.λ))

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
MP(4.0)

julia> MP(1.0) * 3
MP(4.0)

julia> 1 * MP(3.0)
MP(4.0)

julia> MP(1.0) * -Inf
Max-Plus -Inf)
```
"""
Base.:(*)(x::MP, y::MP) = MP(x.λ + y.λ)
Base.:(*)(x::MP, y::Real) = MP(x.λ + y)
Base.:(*)(x::Real, y::MP) = MP(x + y.λ)

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
Max-Plus -Inf

julia> mp0 * 5
Max-Plus -Inf

julia> mp0 + 5
MP(5.0)
```
"""
const global mp0 = mpzero()

"""
    ε (\\varepsilon)

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mpzero`.

Equivalent to ScicosLab code: `%0` sugar notation for `#(-%inf)`

# Examples
```julia-repl
julia> ε
Max-Plus -Inf

julia> ε * 5
Max-Plus -Inf

julia> ε + 5
MP(5.0)
```
"""
const global ε = mp0

"""
    mp1

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mpone`.

Equivalent to ScicosLab code: `%1` sugar notation for `#(1)`

# Examples
```julia-repl
julia> mp1
MP(0.0)

julia> mp1 * 5
MP(5.0)

julia> mp1 + 5
MP(5.0)
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
MP(0.0)

julia> mpe * 5
MP(5.0)

julia> mpe + 5
MP(5.0)
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
MP(Inf)

julia> mptop * 5
MP(Inf)

julia> mptop + 5
MP(Inf)
```
"""
const global mptop = MP(typemax(Float64))

# ==============================================================================
# Max-Plus sparse map function

"""
    mpsparse_map(f, M::SparseMatrixCSC{MP,U})

Map a function ot each element of the Max-Plus sparse matrix.

# Examples
```julia-repl
julia> mpsparse_map(x -> x.λ, mpsparse([1.0 0; 0 1.0]))
2×2 SparseMatrixCSC{Float64, Int64} with 4 stored entries:
 1.0  0.0
 0.0  1.0
```
"""
function mpsparse_map(f, M::SpaMP)
    SparseMatrixCSC(M.m, M.n, M.colptr, M.rowval, map(f, M.nzval))
end

# ==============================================================================
# Algebra conversion: Max-Plus to Min-Plus or Max-Plus to classic algebra

"""
    plustimes(x::MP)

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
plustimes(n::MP) = n.λ

"""
    plustimes(A::ArrMP)

Convert a Max-Plus dense matrix to a dense matrix in standard algebra.

# Examples
```julia-repl
julia> A=[MP(1.0) 2.0; ε mpe]
2×2 Max-Plus dense matrix:
   1.0   2.0
  -Inf   0.0

julia> plustimes(A)
2×2 Matrix{Float64}:
   1.0  2.0
  -Inf  0.0
```
"""
plustimes(A::ArrMP) = map(x -> x.λ, A)

"""
    plustimes(A::SpaMP)

Convert a Max-Plus sparse matrix to an sparse matrix in standard algebra.

# Examples
```julia-repl
julia> S = mpsparse(mpeye(2,2))
2×2 Max-Plus sparse matrix with 2 stored entries:
 0.0   ⋅
  ⋅   0.0

julia> findnz(S)
([1, 2], [1, 2], MP[0.0, 0.0])

julia> plustimes(S)
2×2 SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 0.0   ⋅
  ⋅   0.0

julia> findnz(plustimes(S))
([1, 2], [1, 2], [0.0, 0.0])
```
"""
function plustimes(S::SpaMP)
    SparseMatrixCSC(S.m, S.n, S.colptr, S.rowval, map(x -> x.λ, S.nzval))
end

"""
    minplus(x::MP)

Conversion a Max-Plus number to Min-Plus number. This function convert `+∞` and
`-∞` to their opposite sign.

# Examples
```julia-repl
julia> minplus(mptop), minplus(mp0), minplus(MP(4.5)), minplus(MP(-4.5))
(-Inf, Inf, 4.5, -4.5)

julia> typeof(minplus(mptop))
MP
```
"""
minplus(x::MP) = (x == mp0) ? mptop : ((x == mptop) ? mp0 : x)

"""
    minplus(A::ArrMP)

Conversion a Max-Plus dense matrix to Min-Plus dense matrix. This function
convert `+∞` and `-∞` to their opposite sign.

# Examples
```julia-repl
julia> A = MP([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
3×4 Max-Plus dense matrix:
   0.0   3.0   Inf    1.0
   1.0   2.0   2.0   -Inf
  -Inf   Inf   1.0    0.0

julia> minplus(A)
3×4 Max-Plus dense matrix:
  0.0    3.0   -Inf   1.0
  1.0    2.0    2.0   Inf
  Inf   -Inf    1.0   0.0
```
"""
minplus(A::ArrMP) = map(x -> minplus(x), A)

"""
    minplus(A::SpaMP)

Conversion a Max-Plus dense matrix to Min-Plus dense matrix. This function
convert `+∞` and `-∞` to their opposite sign.

# Examples
```julia-repl
julia> S = mpsparse([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
3×4 Max-Plus sparse matrix with 10 stored entries:
  0     3   Inf   1
  1     2     2   .
  .   Inf     1   0

julia> minplus(S)
3×4 Max-Plus sparse matrix with 10 stored entries:
    0   3   .     1
    1   2   2   Inf
  Inf   .   1     0
```
"""
minplus(S::SpaMP) = dropzeros(map(x -> minplus(x), S))

"""
    sparse(S::SpaMP)

Convert a dense Max-Plus array to a sparse Max-Plus array.

# Examples
```julia-repl
julia> S = mpsparse([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
3×4 Max-Plus sparse matrix with 10 stored entries:
  0     3   Inf   1
  1     2     2   .
  .   Inf     1   0

julia> minplus(S)
3×4 Max-Plus sparse matrix with 12 stored entries:
    0   3   .     1
    1   2   2   Inf
  Inf   .   1     0
```
"""
SparseArrays.sparse(S::SpaMP) = map(x -> x.λ, S)

"""
    full(::SpaMP})

Convert a sparse Max-Plus array to a dense Max-Plus array.
Alternative function name: dense.

# Examples
```julia-repl
julia> full(mpzeros(Float64, 2,5))
2×5 Max-Plus dense array:
 -Inf  -Inf  -Inf  -Inf  -Inf
 -Inf  -Inf  -Inf  -Inf  -Inf
```
"""
full(S::SpaMP) = Matrix(S)

"""
    dense(::SpaMP})

Convert a sparse Max-Plus array to a dense Max-Plus array.
Alternative function name: full.

# Examples
```julia-repl
julia> dense(mpzeros(Float64, 2,5))
2×5 Max-Plus dense array:
 -Inf  -Inf  -Inf  -Inf  -Inf
 -Inf  -Inf  -Inf  -Inf  -Inf
```
"""
dense(S::SpaMP) = Matrix(S)

"""
    array(::SpaMP})

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
array(S::SpaMP) = map(x -> x.λ, S)

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
function mpshow(io::IO, x::MP)
    if (mpstyle == 0)
        show(io, x.λ)
    elseif x == mpzero()
        (mpstyle == 1 || mpstyle == 2) ? (@printf io ".") : (@printf io "ε")
    elseif x == mpone()
        (mpstyle == 1 || mpstyle == 3) ? (@printf io "0") : (@printf io "e")
    elseif x.λ == trunc(x.λ)
        (@printf io "%d" x.λ)
    else
        show(io, x.λ)
    end
end

function mpshow(io::IO, A::ArrMP)
    if (size(A,2) == 1)
        print(io, size(A,1), "-element Max-Plus vector:\n")
    else
        print(io, size(A,1), '×', size(A,2), " Max-Plus dense matrix:\n")
    end
    pretty_table(io, A, tf = tf_borderless, noheader = true)
end

function mpshow(io::IO, S::SpaMP)
    print(io, size(S,1), '×', size(S,2), " Max-Plus sparse matrix with ",
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
function LaTeX(io::IO, A::ArrMP)
    (@printf io "\\left[\n\\begin{array}{*{20}c}\n")
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == mpzero()
                if mpstyle == 0
                    (@printf io "-\\infty")
                elseif mpstyle == 3 || mpstyle == 4
                    (@printf io "\\varepsilon")
                else
                    (@printf io ".")
                end
            elseif A[i,j] == mpone()
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
LaTeX(io::IO, S::SpaMP) = LaTeX(io, full(S))

"""
    show(io::IO, x::MP)

Display a Max-Plus number depending on the currently set style:

# Examples
```julia-repl
julia> mp_change_display(0); mpeye(Float64, 2,2)
2×2 Max-Plus dense array:
  0.0  -Inf
 -Inf   0.0

julia> mp_change_display(0); mpeye(Int64, 2,2)
2×2 Max-Plus dense array:
                    0  -9223372036854775808
 -9223372036854775808                     0

julia> mp_change_display(1); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Max-Plus dense array:
 0  .
 .  0

julia> mp_change_display(2); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Max-Plus dense array:
 e  .
 .  e

julia> mp_change_display(3); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Max-Plus dense array:
 0  ε
 ε  0

julia> mp_change_display(4); mpeye(Float64, 2,2) # or mpeye(Int64, 2,2)
2×2 Max-Plus dense array:
 e  ε
 ε  e
```
"""
# Called by pretty_table() when REPL shows a MP matrix.
# julia> [MP(1) MP(2); MP(3) MP(4)]
Base.show(io::IO, x::MP) = mpshow(io, x)

# Called by REPL when showing a MP scalar.
# julia> MP(1)
function Base.show(io::IO, ::MIME"text/plain", x::MP)
    print(io, "Max-Plus ")
    mpshow(io, x)
end

# Called by the REPL through the display() method. This function fixes
# misaliged columns made by the default show() Julia. We use the package
# PrettyTables.
# julia> [MP(1) MP(2); MP(3) MP(4)]
Base.show(io::IO, ::MIME"text/plain", A::ArrMP) = mpshow(io, A)

Base.show(io::IO, ::MIME"text/plain", S::SpaMP) = mpshow(io, S)

# Convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", A::ArrMP) = LaTeX(io, A)

# Convert a Max-Plus sparse matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", x::SpaMP) = LaTeX(io, A)

# ==============================================================================
# These functions fix Julia bugs

@inline Base.literal_pow(::typeof(^), x::MP, ::Val{0}) = one(x)
@inline Base.literal_pow(::typeof(^), x::MP, ::Val{p}) where p = MP(x.λ * p)
@inline Base.:(^)(x::MP, y::Number) = MP(x.λ * y)
@inline Base.literal_pow(::typeof(^), A::ArrMP, ::Val{p}) where p = A^p
@inline Base.abs2(x::MP) = x.λ + x.λ

# ==============================================================================
# Hardly used operators

"""
    /(x::MP, y::MP)

Divisor operator. Return the difference between `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) / MP(2.0)
MP(-1.0)
```
"""
Base.:(/)(x::MP, y::MP) = MP(x.λ - y.λ)
Base.:(/)(x::MP, y::Real) = MP(x.λ - b)
Base.:(/)(x::Real, y::MP) = MP(a - y.λ)

"""
    -(x::MP, y::MP)

Minus operator (not used for Max-Plus but for Min-Plus).
Return the difference between `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) - MP(2.0)
MP(-1.0)
```
"""
Base.:(-)(x::MP, y::MP) = MP(x.λ - y.λ)
Base.:(-)(x::MP, y::Real) = MP(x.λ - y)
Base.:(-)(x::Real, y::MP) = MP(x - y.λ)

"""
    min(x::MP, y::MP)

Return the minimun of Max-Plus numbers `x` and `y`.

# Examples
```julia-repl
julia> min(MP(1), 3)
MP:
  1

julia> min(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 Max-Plus dense matrix:
  4   1
  6   1
```
"""
Base.:min(x::MP, y::MP) = MP(min(x.λ, y.λ))
Base.:min(x::MP, y::Real) = MP(min(x.λ, y))
Base.:min(x::Real, y::MP) = MP(min(x, y.λ))
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

"""
    max(x::MP, y::MP)

Return the max of `x` and `y`.

# Examples
```julia-repl
julia> max(MP(1), 3)
MP:
  3

julia> max(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 Max-Plus dense matrix:
  10   5
  10   5
```
"""
Base.:max(x::MP, y::MP) = MP(max(x.λ, y.λ))
Base.:max(x::MP, y::Real) = MP(max(x.λ, y))
Base.:max(x::Real, y::MP) = MP(max(x, y.λ))
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
    mpeye(n::Int64)

Construct a Max-Plus identity dense n-by-n matrix.

# Examples
```julia-repl
julia> mpeye(2)
2×2 Max-Plus dense matrix:
   0.0   -Inf
  -Inf    0.0
```
"""
mpeye(n::Int64) = Matrix{MP}(mpI, n, n)

"""
    mpeye(m::Int64, n::Int64)

Construct a Max-Plus identity dense m-by-n matrix.

# Examples
```julia-repl
julia> mpeye(2, 3)
2×3 Max-Plus dense matrix:
   0.0   -Inf   -Inf
  -Inf    0.0   -Inf
```
"""
mpeye(m::Int64, n::Int64) = Matrix{MP}(mpI, m, n)

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
2×2 Max-Plus dense matrix:
  0.0  -Inf
 -Inf   0.0

julia> mpeye(MP(A))
2×2 Max-Plus dense matrix:
  0.0  -Inf
 -Inf   0.0
```
"""
mpeye(A::Array) = Matrix{MP}(mpI, size(A,1), size(A,2))
mpeye(A::ArrMP) = Matrix{MP}(mpI, size(A,1), size(A,2))

# ==============================================================================
# Zero matrix

"""
    mpzeros(n::Int64)

Construct a Max-Plus zero n-by-m sparse matrix.

# Examples
```julia-repl
julia> mpzeros(2)
2-element SparseVector{MP,Int64} with 0 stored entries
```
"""
mpzeros(n::Int64) = spzeros(MP, n)

"""
    mpzeros(n::Int64)

Construct a sparse Max-Plus zero m-by-n matrix.

# Examples
```julia-repl
julia> mpzeros(2,5)
2×5 SparseMatrixCSC{MP,Int64} with 0 stored entries
```
"""
mpzeros(m::Int64, n::Int64) = spzeros(MP, m, n)

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
2×2 Max-Plus sparse matrix with 0 stored entries:
  -Inf   -Inf
  -Inf   -Inf

julia> mpzeros(MP(A))
2×2 Max-Plus sparse matrix with 0 stored entries:
  -Inf   -Inf
  -Inf   -Inf
```
"""
mpzeros(A::Array) = spzeros(MP, size(A,1), size(A,2))
mpzeros(A::ArrMP) = spzeros(MP, size(A,1), size(A,2))

# ==============================================================================
# One matrix

"""
    mpones(n::Int64)

Construct a Max-Plus one n-by-1 matrix.

# Examples
```julia-repl
julia> mpones(2)
2×1 Max-Plus dense matrix:
  0.0
  0.0
```
"""
mpones(n::Int64) = ones(MP, n)

"""
    mpones(m::Int64, n::Int64)

Construct a Max-Plus one m-by-n matrix.

# Examples
```julia-repl
julia> mpones(3,2)
3×2 Max-Plus dense matrix:
  0.0   0.0
  0.0   0.0
  0.0   0.0
```
"""
mpones(m::Int64, n::Int64) = ones(MP, m, n)

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
2×2 Max-Plus dense matrix:
  0.0   0.0
  0.0   0.0
```
"""
mpones(A::Array) = mpones(Float64, size(A,1), size(A,2))
mpones(A::ArrMP) = mpones(Float64, size(A,1), size(A,2))

# ==============================================================================
# Dense and Sparse matrices

"""
    mpsparse(A::Array{T}; keepzeros::Bool=false)

Transform a dense matrix from classic algebra to a Max-Plus sparse matrix.
Max-Plus zeros (`ε`, `mp0`, `Max-Plus -Inf)`) are removed and classic algebra
zeros are removed if the argument `keepzeros` is set to `true`.

# Arguments
- keepzeros: if true then 0 values from classic algebra are considered as values
and are not removed.
- keepzeros: if false then 0 values from classic algebra are removed.

# Examples
```julia-repl
julia> S = mpsparse([-Inf 0; 0 -Inf])
2×2 Max-Plus sparse matrix with 2 stored entries:
  .   0
  0   .

julia> findnz(S)
([2, 1], [1, 2], MP[0.0, 0.0])

julia> S = mpsparse([-Inf 0; 0 -Inf], keepzeros=false)
2×2 Max-Plus sparse matrix with 0 stored entries:
  .   .
  .   .

julia> findnz(S)
(Int64[], Int64[], MP[])
```
"""
mpsparse(A::Array; keepzeros::Bool=true) =
    keepzeros ? sparse(MP(A)) : MP(sparse(A), keepzeros=false)

"""
    mpsparse(A::ArrMP)

Transform a dense Max-Plus matrix to a sparse Max-Plus matrix.
Zero values from Max-Plus algebra are removed (therefore -Inf
from classic algebra are removed).

# Examples
```julia-repl
julia> S = mpsparse(MP([4 0; 7 -Inf]))
2×2 Max-Plus sparse matrix with 3 stored entries:
  4   0
  7   .

julia> findnz(S)
([1, 2, 1], [1, 1, 2], MP[4, 7, 0])
```
"""
mpsparse(A::ArrMP) = sparse(A)

# ==============================================================================
# Max-Plus matrices operations

"""
    mptrace(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> mptrace([MP(1) 2; 3 4])
MP(4)
```
"""
mptrace(A::ArrMP) = isempty(A) ? mp0 : sum(diag(A))
mptrace(A::Array) = isempty(A) ? mp0 : sum(MP(diag(A)))
mptrace(S::SpaMP) = isempty(S) ? mp0 : sum(diag(S))
mptrace(S::Sparse) = sum(MP(diag(S)))

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
1×2 Max-Plus dense array:
 399  399
```
"""
mpnorm(A::Array) = A[argmax(A)].λ - A[argmin(A)].λ
mpnorm(S::SpaMP) = S[argmax(S)].λ - S[argmin(S)].λ

# ==============================================================================
# Max-Plus star FIXME not working with int

"""
    B = mpstar(A::ArrMP)

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
MP:
  Inf

julia> mpstar(MP(-1.0))
MP:
  0

julia> mpstar(MP([1.0 2; 3 4]))
2×2 Max-Plus dense matrix:
  Inf   Inf
  Inf   Inf

julia> A = mpstar(MP([-3.0 -2;-1 0])); B = mpstar(A)
2×2 Max-Plus dense matrix:
   0   -2
  -1    0

julia> B == (mpeye(2,2) + A)
true

julia> B == (B + A * A)
true
```
"""
function mpstar(A::ArrMP)
    n = size(A, 1)
    if n != size(A, 2)
        error("Matrix shall be squared")
    end
    C = A
    for k in 1:n
        t = (C[k,k].λ <= zero(Float64)) ? zero(Float64) : typemax(Float64);
        for j in 1:n, i in 1:n
            C[i,j] = MP(max(C[i,j].λ, C[i,k].λ + C[k,j].λ + t))
        end
    end

    for k in 1:n
        C[k,k] = MP(max(C[k,k].λ, zero(Float64)));
    end
    C
end

"""
    mpstar(x::MP)

Make x a 1x1 matrix then call mpstar(A::ArrMP).
See mpstar(A::ArrMP) for more information.
"""
mpstar(x::MP) = mpstar([x])[1,1]

"""
    mpastarb(A::ArrMP, b::ArrMP)

Max-Plus linear system solution.

Solve `x = Ax + b` in the Max-Plus algebra when there is no circuits with
positive weight in `G(A')` (the incidence graph of `A'`, that is it exists an
arc from `j` to `i` if `A_ij` is nonzero).

TODO It is much more efficient in time and memory than `mpstar(A) * b`.
"""
mpastarb(A::ArrMP, b::ArrMP) = mpstar(A) * b
# TODO: optimal code

"""
    B = mpplus(A::ArrMP)

Compute `A * A^* = A + A^2 + ...` of a maxplus matrix A.
See also mpstar.

# Arguments
- A : Max-Plus full square matrix.
- B : Id + A + A^2 + ...

# Examples
```julia-repl
julia> mpplus(MP(1.0))
MP:
  Inf

julia> mpplus(MP(-1.0))
MP:
  -1

julia> A = MP([-3.0 -2; -1 0]); B = mpplus(A)
2×2 Max-Plus dense matrix:
  -3   -2
  -1    0

julia> B == (A * mpstar(A))
true
```
"""
function mpplus(A::ArrMP)
    n = size(A, 1)
    if n != size(A, 2)
        error("Matrix shall be squared")
    end
    C = A
    for k in 1:n
        t = (C[k,k].λ <= zero(Float64)) ? zero(Float64) : typemax(Float64);
        for j in 1:n, i in 1:n
            C[i,j] = MP(max(C[i,j].λ, C[i,k].λ + C[k,j].λ + t))
        end
    end
    C
end

"""
    mpplus(x::MP)

Make x a 1x1 matrix then call mpplus(A::ArrMP).
See mpplus(A::ArrMP) for more information.
"""
mpplus(x::MP) = mpplus([x])[1,1]

# ==============================================================================
# TODO insertion out of bounds => Sparse:  d=MP([1.0 2; 3 4]); d[5,5] = MP(6.0)

# ==============================================================================
include("fallbacks.jl")
include("howard.jl")
include("syslin.jl")
#include("flowshop.jl")

end # MaxPlus module
