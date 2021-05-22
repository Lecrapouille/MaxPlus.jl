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
    mpeye, mpzeros, mpones, full, dense, array, plustimes, inv, mpsparse_map,
    mptrace, mpnorm, mpastarb, mpstar, mpplus, howard, mp_change_display,
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
or Int64 to a number in the tropical semi-ring max,+ (ℝ ∪ {-∞}, ⨁, ⨂) where ℝ is
the domain of reals, ⨁ is the usual multiplication and ⨂ is the usual maximum.

`MP(3)` is the equivalent to ScicosLab code: `#(3)`

# Examples
```julia-repl
julia> a = MP(3.5)
Max-Plus 3.5

julia> b = MP(3)
Max-Plus 3

julia> typeof(a), typeof(b)
(MP, MP)
```
"""
struct MP <: Real
    λ::Float64

    # Avoid pathological cases: MP(Nan) made with mp0 * mptop
    MP(x::Real) = isnan(x) ? new(typemin(Float64)) : new(Float64(x))
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
    MP(A::SparseArray)

Convert a sparse array from classic algebra to a Max-Plus sparse array. By
default, explicit Max-Plus zeros (`ε`, `mp0`, `MP(-Inf)`) are removed except if
the parameter `keepzeros` is set to `true`.

# Examples
```julia-repl
julia> using SparseArrays

julia> A = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
3×3 Max-Plus sparse matrix with 2 stored entries:
  .   .   .
  .   2   .
  .   .   0

julia> A.nzval
2-element Max-Plus vector:
  2.0
  0.0
```
"""
MP(S::Sparse{T,U}) where {T,U} = convert(SpaMP{U}, S)

"""
    MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)

Construct a sparse Max-Plus matrix such as S[I[k], J[k]] = V[k]

# Examples
```julia-repl
julia> MP([1; 2; 3], [1; 2; 3], [42; 43; 44])
3×3 Max-Plus sparse matrix with 3 stored entries:
  42    .    .
   .   43    .
   .    .   44
```
"""
function MP(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {Tv,Ti<:Integer}
    sparse(I, J, MP(V))
end

# ==============================================================================
# Constructor from non Max-Plus sparse vector

"""
    MP(A::SparseVector)

Convert a sparse vector from classic algebra to a Max-Plus sparse vector. By
default, explicit Max-Plus zeros (`ε`, `mp0`, `MP(-Inf)`) are removed except if
the parameter `keepzeros` is set to `true`.


# Examples
```julia-repl
julia> MP(sparse([1.0, 0.0, 1.0]))
3-element SparseVector{MP{Float64},Int64} with 2 stored entries:
  [1]  =  1.0
  [3]  =  1.0

julia> A = MP(sparse([1; -Inf; 3]))
3-element SparseVector{MP, Int64} with 2 stored entries:
  [1]  =  1.0
  [3]  =  3.0
```
"""
MP(V::SparseVector{T,U}) where {T, U} = convert(SparseVector{MP,U}, V)

# ==============================================================================
# Constructor from non Max-Plus range

"""
    MP(x::UnitRange)

Create a Max-Plus dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1:3)
3-element Max-Plus Vector:
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
5-element Max-Plus Vector:
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
Max-Plus 0.0
```
"""
mpone() = MP(zero(Float64))

# ==============================================================================
# Max-Plus conversion from Bool. Look weird but needed for Matrix{MP}(I, 2, 2)

"""
    MP(x::Bool)

Force conversion from Bool to Max-Plus needed for example by the identity matrix
operator `LinearAlgebra.I` else identity matrices are not well created.

# Examples
```julia-repl
julia> MP(true)
Max-Plus 0.0

julia> MP(false)
Max-Plus -Inf

julia> using LinearAlgebra

julia> Matrix{MP}(I, 2, 2)
2×2 Max-Plus dense matrix:
   0.0   -Inf
  -Inf    0.0
```
"""
MP(x::Bool) = x ? mpone() : mpzero()

# ==============================================================================

"""
    mpI

Is the equivalent of `LinearAlgebra.I` but for Max-Plus type. Allow to create
identity Max-Plus matrices.

# Examples
```julia-repl
julia> typeof(mpI)
LinearAlgebra.UniformScaling{MP}

julia> Matrix(mpI, 2, 2)
2×2 Max-Plus dense matrix:
   0.0   -Inf
  -Inf    0.0

julia> using LinearAlgebra

julia> Matrix(mpI, 2, 2) == Matrix{MP}(I, 2, 2)
true
```
"""
const global mpI = UniformScaling(mpone())

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
Max-Plus 3.0

julia> MP(1.0) + 3
Max-Plus 3.0

julia> 1 + MP(3.0)
Max-Plus 3.0

julia> MP(3.0) + -Inf
Max-Plus 3.0
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
Max-Plus 4.0

julia> MP(1.0) * 3
Max-Plus 4.0

julia> 1 * MP(3.0)
Max-Plus 4.0

julia> MP(1.0) * -Inf
Max-Plus -Inf
```
"""
Base.:(*)(x::MP, y::MP) = MP(x.λ + y.λ)
Base.:(*)(x::MP, y::Real) = MP(x.λ + y)
Base.:(*)(x::Real, y::MP) = MP(x + y.λ)

# ==============================================================================
# Max-Plus divisor operator

"""
    /(x::MP, y::MP)

Divisor operator. Return the difference between `x` and `y` in classic algebra.

x ⨸ y = x ⨂ y^-1

# Examples
```julia-repl
julia> MP(1.0) / MP(2.0)
Max-Plus -1.0
```
"""
Base.:(/)(x::MP, y::MP) = MP(x.λ - y.λ)
Base.:(/)(x::MP, y::Real) = MP(x.λ - b)
Base.:(/)(x::Real, y::MP) = MP(a - y.λ)

# ==============================================================================
# Max-Plus substraction operator

"""
    -(x::MP, y::MP)

Minus operator (not used for Max-Plus but for Min-Plus).
Return the difference between `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) - MP(2.0)
Max-Plus -1.0

julia> MP(1.0) / MP(0.0)
Max-Plus 1.0

julia> MP(1.0) / mp0
Max-Plus Inf
```
"""
Base.:(-)(x::MP, y::MP) = MP(x.λ - y.λ)
Base.:(-)(x::MP, y::Real) = ((x == mp0) && (y == typemin(Real))) ? mp0 : MP(x.λ - y)
Base.:(-)(x::Real, y::MP) = ((x == typemin(Real)) && (y == mp0)) ? mp0 : MP(x - y.λ)
Base.:(-)(x::MP) = (x == mpzero()) ? mpzero() : MP(-x.λ)

# ==============================================================================
# Max-Plus residual operator

"""
    inv(A::ArrMP)

Return the inverse of a scalar or a Max-Plus matrix which is `transpose(-A)` if
the matrix can be inversed. Throw an error if the matrix cannot be inversed.

# Example 1: Scalar
```julia-repl
julia> inv(MP(5))
Max-Plus -5.0

julia> MP(5)^-1
Max-Plus -5.0
```

# Example 2: Square matrix inversible
```julia-repl
julia> A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
3×3 Max-Plus dense matrix:
  .   1   .
  2   .   .
  .   .   3

julia> A^-1
3×3 Max-Plus transpose(dense matrix):
   .   -2    .
  -1    .    .
   .    .   -3

julia> inv(A)
3×3 Max-Plus transpose(dense matrix):
   .   -2    .
  -1    .    .
   .    .   -3

julia> A * inv(A)
3×3 Max-Plus dense matrix:
  0   .   .
  .   0   .
  .   .   0

julia> A * inv(A) == inv(A) * A == mpeye(A)
true
```

# Example 3: Square matrix inversible
```julia-repl
julia> A = [mp0 1 mp0; 2 mp0 mp0]
2×3 Max-Plus dense matrix:
  .   1   .
  2   .   .

julia> A^-1
3×2 Max-Plus transposed dense matrix:
   .   -2
  -1    .
   .    .

julia> inv(A)
3×2 Max-Plus transposed dense matrix:
   .   -2
  -1    .
   .    .

julia> A * inv(A)
2×2 Max-Plus dense matrix:
  0   .
  .   0

julia> inv(A) * A
3×3 Max-Plus dense matrix:
  0   .   .
  .   0   .
  .   .   .

julia> (A * inv(A)) != (inv(A) * A)
true
```

# Example 4: Non inversible
```julia-repl
julia> A = MP([1 2; 3 4])
2×2 Max-Plus dense matrix:
  1.0   2.0
  3.0   4.0

julia> inv(A)
ERROR: The matrix cannot be inversed
```
"""
function Base.inv(A::ArrMP)
    isempty(A) && return A

    Z = transpose(-A)

    # Checks
    (A * Z != mpeye(size(A,1), size(Z,2))) && error("The matrix cannot be inversed")
    if (size(A,1) == size(A,2))
        (Z * A != mpeye(size(Z,1), size(A,2))) && error("The matrix cannot be inversed")
    end

    return Z
end

"""
    \\(A::ArrMP, b::ArrMP)

`x = A \\ b` is a solution to `A ⨂ x = b` and its simply computed as
`inv(A) * b`.

# Examples
```julia-repl
julia> A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
3×3 Max-Plus dense matrix:
  .   1   .
  2   .   .
  .   .   3

julia> B = [3 mp0 mp0; mp0 mp0 4; mp0 5 mp0]
3×3 Max-Plus dense matrix:
  3   .   .
  .   .   4
  .   5   .

julia> x = A \\ B
3×3 Max-Plus dense matrix:
  .   .   2
  2   .   .
  .   2   .

julia> A * x == B
true

julia> A \\ A
3×3 Max-Plus dense matrix:
  0   .   .
  .   0   .
  .   .   0
```
"""
Base.:(\)(A::ArrMP, b::ArrMP) = inv(A) * b

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
Max-Plus 5.0
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
Max-Plus 5.0
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
Max-Plus 0.0

julia> mp1 * 5
Max-Plus 5.0

julia> mp1 + 5
Max-Plus 5.0
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
Max-Plus 0.0

julia> mpe * 5
Max-Plus 5.0

julia> mpe + 5
Max-Plus 5.0
```
"""
const global mpe = mp1

# ==============================================================================
# Min-Plus constante

"""
    mptop

Create the constant Min-Plus one (equals to `+∞`, infinity) in classic algebra which is
the neutral for operators ⨁ and ⨂.

Equivalent to ScicosLab code: `%top = #(%inf)`

# Examples
```julia-repl
julia> mptop
Max-Plus Inf

julia> mptop * 5
Max-Plus Inf

julia> mptop + 5
Max-Plus Inf
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
julia> mpsparse_map(x -> x.λ, sparse(MP([1.0 0; 0 1.0])))
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
Float64
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
julia> S = sparse(mpeye(2,2))
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
    full(::SpaMP})

Convert a sparse Max-Plus array to a dense Max-Plus array.
Alternative function name: dense.

# Examples
```julia-repl
julia> full(mpzeros(2,5))
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
julia> dense(mpzeros(2,5))
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
julia> array(mpzeros(2,2))
2×2 Matrix{Float64}:
 -Inf  -Inf
 -Inf  -Inf
```
"""
array(S::SpaMP) = Matrix(map(x -> x.λ, S))

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
function mp_change_display(style::Int)
    global mpstyle = min(max(style, 0), 4)
end

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

function mpshow(io::IO, A::LinearAlgebra.Transpose{MP, Matrix{MP}})
    print(io, size(A,1), '×', size(A,2), " Max-Plus transposed dense matrix:\n")
    pretty_table(io, A, tf = tf_borderless, noheader = true)
end

function mpshow(io::IO, V::LinearAlgebra.Transpose{MP, Vector{MP}})
    print(io, size(A,1), "-element Max-Plus transposed vector:\n")
    pretty_table(io, V, tf = tf_borderless, noheader = true)
end

# Note: we force the sparse display like done in Julia 1.5 because since Julia
# 1.6 sparse matrices are displayed like dense matrices with dots for zeros.
# This sounds weird since displayed huge sparse matrices take the same space
# than dense matrix.
mpshow(io::IO, S::SpaMP) = show(io, S)

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

function LaTeX(A::ArrMP)
    s = "\\left[\n\\begin{array}{*{20}c}\n"
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == mpzero()
                if mpstyle == 0
                    s = s * "-\\infty"
                elseif mpstyle == 3 || mpstyle == 4
                    s = s * "\\varepsilon"
                else
                    s = s * "."
                end
            elseif A[i,j] == mpone()
                if mpstyle == 2 || mpstyle == 4
                    s = s * "e"
                else
                    s = s * string(Int64(A[i,j].λ))
                end
            elseif A[i,j].λ == trunc(A[i,j].λ)
                s = s * string(Int64(A[i,j].λ))
            else
                s = s * string(A[i,j].λ)
            end
            if j < size(A, 2)
                s = s * " & "
            end
        end
        s = s * " \\\\\n"
    end
    return s * "\\end{array}\n\\right]\n"
end

"""
    show(io::IO, x::MP)

Display a Max-Plus number depending on the currently set style:

# Examples
```julia-repl
julia> mp_change_display(0); mpeye(2,2)
2×2 Max-Plus dense array:
  0.0  -Inf
 -Inf   0.0

julia> mp_change_display(1); mpeye(2,2)
2×2 Max-Plus dense array:
 0  .
 .  0

julia> mp_change_display(2); mpeye(2,2)
2×2 Max-Plus dense array:
 e  .
 .  e

julia> mp_change_display(3); mpeye(2,2)
2×2 Max-Plus dense array:
 0  ε
 ε  0

julia> mp_change_display(4); mpeye(2,2)
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

# Convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", A::ArrMP) = LaTeX(io, A)

# Convert a Max-Plus sparse matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", x::SpaMP) = LaTeX(io, A)

Base.show(io::IO, ::MIME"text/plain", A::LinearAlgebra.Transpose{MP, Matrix{MP}}) = mpshow(io, A)
Base.show(io::IO, ::MIME"text/plain", V::LinearAlgebra.Transpose{MP, Vector{MP}}) = mpshow(io, V)

# ==============================================================================
# These functions fix Julia bugs

@inline Base.literal_pow(::typeof(^), x::MP, ::Val{0}) = one(x)
@inline Base.literal_pow(::typeof(^), x::MP, ::Val{p}) where p = MP(x.λ * p)
@inline Base.:(^)(x::MP, y::Number) = MP(x.λ * y)
@inline Base.literal_pow(::typeof(^), A::ArrMP, ::Val{p}) where p = A^p
@inline Base.literal_pow(::typeof(^), x::ArrMP, ::Val{-1}) = transpose(-x)
@inline Base.abs(x::MP) = MP(abs(x.λ))
@inline Base.abs2(x::MP) = x.λ + x.λ
@inline Base.float(x::MP) = x.λ
@inline Base.round(x::MP, n::RoundingMode) = MP(round(x.λ, n))
@inline Base.big(x::MP) = Base.big(x.λ)
@inline Base.sign(x::MP) = Base.sign(x.λ)

# ==============================================================================
# Hardly used operators

"""
    min(x::MP, y::MP)

Return the minimun of Max-Plus numbers `x` and `y`.

min(x, y) = (x ⨂ y) ⨸ (x ⨁ y)

# Examples
```julia-repl
julia> min(MP(1), 3)
Max-Plus 1

julia> min(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 Max-Plus dense matrix:
  4   1
  6   1
```
"""
Base.:min(x::MP, y::MP) = (x * y) / (x + y)
Base.:min(x::MP, y::Real) = min(x, MP(y))
Base.:min(x::Real, y::MP) = min(MP(x), y)
Base.:min(A::ArrMP, B::ArrMP) = map(Base.:min, A, B)
Base.:min(A::SpaMP, B::SpaMP) = map(Base.:min, A, B)

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
  .   .   .   .   .
  .   .   .   .   .

julia> full(mpzeros(2,5))
2×5 Max-Plus dense matrix:
  -Inf   -Inf   -Inf   -Inf   -Inf
  -Inf   -Inf   -Inf   -Inf   -Inf
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
  .   .
  .   .

julia> mpzeros(MP(A))
2×2 Max-Plus sparse matrix with 0 stored entries:
  .   .
  .   .
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
2-element Max-Plus vector:
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
mpones(A::Array) = mpones(size(A,1), size(A,2))
mpones(A::ArrMP) = mpones(size(A,1), size(A,2))

# ==============================================================================
# Max-Plus matrices operations

"""
    mptrace(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> mptrace([MP(1) 2; 3 4])
Max-Plus 4.0
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
julia> A = MP([1 20 2;30 400 4;4 50 10])
3×3 Max-Plus dense matrix:
   1.0    20.0    2.0
  30.0   400.0    4.0
   4.0    50.0   10.0

julia> S = MP(sparse(A))
3×3 Max-Plus sparse matrix with 9 stored entries:
   1    20    2
  30   400    4
   4    50   10

julia> mpnorm(A)
Max-Plus 399

julia> mpnorm(S)
Max-Plus 399
```
"""
mpnorm(A::Array) = MP(A[argmax(A)].λ - A[argmin(A)].λ)
mpnorm(S::SpaMP) = MP(S[argmax(S)].λ - S[argmin(S)].λ)

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
Max-Plus Inf

julia> mpstar(MP(-1.0))
Max-Plus 0

julia> mpstar(MP([1.0 2; 3 4]))
2×2 Max-Plus dense matrix:
  Inf   Inf
  Inf   Inf

julia> A = mpstar(MP([-3.0 -2;-1 0]))
2×2 Max-Plus dense matrix:
   0   -2
  -1    0

julia> B = mpstar(A)
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
Max-Plus Inf

julia> mpplus(MP(-1.0))
Max-Plus -1

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
