# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# A portage of the ScicosLab Max-Plus toolbox http://www.scicoslab.org/
# License: public domain
#
# Note: the documentation of functions for the REPL are placed in docstrings.jl
# ==============================================================================
# Docstring for functions
# ==============================================================================

"""
    MP

Immutable Julia structure for Max-Plus scalar. Promote a number of type Float64
or Int64 to a number in the tropical semi-ring (max, +) (ℝ ∪ {-∞}, ⨁, ⨂) where ℝ
is the domain of reals, ⨁ is the usual multiplication and ⨂ is the usual
maximum.

# Notes
`MP(3)` is the equivalent to ScicosLab code: `#(3)`

# Examples
```julia-repl
julia> a = MP(3.5)
Max-Plus 3.5

julia> typeof(a)
MP

julia> a = MP(3)
Max-Plus 3

julia> typeof(a)
MP
```
"""
MP

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
MP(x::Bool)

"""
    MP(A::Array)

Convert a dense array from classic algebra to a Max-Plus dense array.
Note: if the array alaready contains Max-Plus number then the `MP()`
is not necessary since Max-Plus number contaminate other numbers.

# Examples (constructor)
```julia-repl
julia> A = MP([1.0 -Inf; 0.0 4])
2×2 Max-Plus dense matrix:
  1.0   -Inf
  0.0    4.0

julia> typeof(A)
Matrix{MP} (alias for Array{MP, 2})
```

# Examples (constructor not needed)
```julia-repl
julia> A = [MP(1.0) -Inf; 0.0 4]
2×2 Max-Plus dense matrix:
  1.0   -Inf
  0.0    4.0

julia> typeof(A)
Matrix{MP} (alias for Array{MP, 2})
```
"""
MP(A::Array)

"""
    MP(A::SparseMatrixCSC)

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
MP(S::SparseMatrixCSC)

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
MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)

"""
    MP(V::SparseVector)

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
MP(V::SparseVector)

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
MP(x::UnitRange)

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
MP(x::StepRangeLen)

"""
    zero(::Type{MP})

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.
"""
Base.zero(::Type{MP})

"""
    zero(::MP)

Create the constant Max-Plus zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.
"""
Base.zero(x::MP)

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
mpzero()

"""
    one(::MP)

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
    the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.
"""
Base.one(::Type{MP})

"""
    one(::MP)

Create the constant Max-Plus one (equals to `0`, zero) in classic algebra which is
    the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.
"""
Base.one(x::MP)

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
mpone()

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
mpI

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
Base.:(+)(x::MP, y::MP)

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
Base.:(*)(x::MP, y::MP)

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
Base.:(/)(x::MP, y::MP)

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
Base.:(-)(x::MP, y::MP)

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
Base.inv(A::ArrMP)

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
Base.:(\)(A::ArrMP, b::ArrMP)

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
mp0

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
ε

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
mp1

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
mpe


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
mptop

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
mpsparse_map(f, M::SpaMP)

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
plustimes(n::MP)

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
plustimes(A::ArrMP)

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
plustimes(S::SpaMP)

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
full(S::SpaMP)

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
dense(S::SpaMP)

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
array(S::SpaMP)

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
Base.:min(x::MP, y::MP)

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
mpeye(n::Int64)

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
mpeye(m::Int64, n::Int64)

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
mpeye(A::Array)

"""
    mpzeros(n::Int64)

Construct a Max-Plus zero n-by-m sparse matrix.

# Examples
```julia-repl
julia> mpzeros(2)
2-element SparseVector{MP,Int64} with 0 stored entries
```
"""
mpzeros(n::Int64)

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
mpzeros(m::Int64, n::Int64)

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
mpzeros(A::Array)

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
mpones(n::Int64)

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
mpones(m::Int64, n::Int64)

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
mpones(A::Array)

"""
    mptrace(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> mptrace([MP(1) 2; 3 4])
Max-Plus 4.0
```
"""
mptrace(A::ArrMP)

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
mpnorm(A::Array)

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
mpstar(A::ArrMP)

"""
    mpstar(x::MP)

Make x a 1x1 matrix then call mpstar(A::ArrMP).
See mpstar(A::ArrMP) for more information.
"""
mpstar(x::MP)

"""
    mpastarb(A::ArrMP, b::ArrMP)

Max-Plus linear system solution.

Solve `x = Ax + b` in the Max-Plus algebra when there is no circuits with
positive weight in `G(A')` (the incidence graph of `A'`, that is it exists an
arc from `j` to `i` if `A_ij` is nonzero).

TODO It is much more efficient in time and memory than `mpstar(A) * b`.
"""
mpastarb(A::ArrMP, b::ArrMP)

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
mpplus(A::ArrMP)

"""
    mpplus(x::MP)

Make x a 1x1 matrix then call mpplus(A::ArrMP).
See mpplus(A::ArrMP) for more information.
"""
mpplus(x::MP)
