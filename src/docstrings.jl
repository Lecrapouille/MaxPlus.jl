# ==============================================================================
# Docstring for functions
# ==============================================================================

"""
    MP

Immutable Julia structure for (max,+) scalar. Promote a number (i.e. type
Float64 or Int64) to a number in the tropical semi-ring (max, +) (ℝ ∪ {-∞}, ⨁,
⨂) where ℝ is the domain of reals, ⨁ is the usual multiplication and ⨂ is the
usual maximum.

# Notes
`MP(3)` is the equivalent to ScicosLab code: `#(3)`

# Examples
```julia-repl
julia> a = MP(3.5)
(max,+) 3.5

julia> typeof(a)
MP (alias for Trop{MaxPlus.Max})

julia> a = MP(3)
(max,+) 3

julia> typeof(a)
MP (alias for Trop{MaxPlus.Max})
```
"""
MP

"""
    MP(x::Bool)

Force conversion from Bool to (max,+) needed for example by the identity matrix
operator `LinearAlgebra.I` else identity matrices are not well created.

# Examples
```julia-repl
julia> mp_change_display(0)

julia> MP(true)
(max,+) 0.0

julia> MP(false)
(max,+) -Inf

julia> using LinearAlgebra

julia> I
UniformScaling{Bool}
true*I

julia> Matrix{MP}(I, 2, 2)
2×2 (max,+) dense matrix:
   0.0   -Inf
  -Inf    0.0

julia> mp_change_display(1)

julia> MP(true)
(max,+) 0

julia> MP(false)
(max,+) .
```
"""
MP(x::Bool)

"""
    MP(A::Array)

Convert a dense array from classic algebra to a (max,+) dense array.
Note: if the array alaready contains (max,+) number then the `MP()`
is not necessary since (max,+) number contaminate other numbers.

# Examples (constructor)
```julia-repl
julia> A = MP([1.0 -Inf; 0.0 4])
2×2 (max,+) dense matrix:
  1   .
  0   4

julia> typeof(A)
Matrix{MP} (alias for Array{Trop{MaxPlus.Max}, 2}
```

# Examples (constructor not needed)
```julia-repl
julia> A = [MP(1.0) -Inf; 0.0 4]
2×2 (max,+) dense matrix:
  1   .
  0   4

julia> typeof(A)
Matrix{MP} (alias for Array{Trop{MaxPlus.Max}, 2})
```
"""
MP(A::Array)

"""
    MP(A::SparseMatrixCSC)

Convert a sparse array from classic algebra to a (max,+) sparse array. By
default, and following Julia rules, explicit (max,+) zeros (`ε`, `mp0`,
`MP(-Inf)`) are not removed. Call dropzeros() to remove them.

# Examples
```julia-repl
julia> using SparseArrays

julia> A = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
3×3 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  .
  [2, 2]  =  2
  [3, 3]  =  0

julia> dropzeros(A)
3×3 (max,+) sparse matrix with 2 stored entries:
  [2, 2]  =  2
  [3, 3]  =  0

julia> A = sparse([1, 2, 3], [1, 2, 3], MP([-Inf, 2, 0]))
3×3 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  .
  [2, 2]  =  2
  [3, 3]  =  0

julia> dropzeros(A)
3×3 (max,+) sparse matrix with 2 stored entries:
  [2, 2]  =  2
  [3, 3]  =  0
```
"""
MP(S::SparseMatrixCSC)

"""
    MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)

Construct a sparse (max,+) matrix such as S[I[k], J[k]] = V[k]. By default, and
following Julia rules, explicit (max,+) zeros (`ε`, `mp0`, `MP(-Inf)`) are not
removed. Call dropzeros() to remove them.

# Examples
```julia-repl
julia> A = MP([1; 2; 3], [1; 2; 3], [42.5; -Inf; 44])
3×3 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  42.5
  [2, 2]  =  .
  [3, 3]  =  44

julia> dropzeros(A)
3×3 (max,+) sparse matrix with 2 stored entries:
  [1, 1]  =  42.5
  [3, 3]  =  44
```
"""
MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)

"""
    MP(V::SparseVector)

Convert a sparse vector from classic algebra to a (max,+) sparse vector. By
default, explicit (max,+) zeros (`ε`, `mp0`, `MP(-Inf)`) are removed.

# Examples
```julia-repl
julia> MP(sparse([1.0, 0.0, 1.0]))
3-element SparseVector{MP, Int64} with 2 stored entries:
  [1]  =  1
  [3]  =  1

julia> A = sparse(MP([1; -Inf; 3]))
3-element SparseVector{MP, Int64} with 2 stored entries:
  [1]  =  1
  [3]  =  3

julia> A = MP(sparse([1; -Inf; 3]))
3-element SparseVector{MP, Int64} with 3 stored entries:
  [1]  =  1
  [2]  =  .
  [3]  =  3
```
"""
MP(V::SparseVector)

"""
    MP(x::UnitRange)

Create a (max,+) dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1:3)
3-element (max,+) vector:
  1
  2
  3
```
"""
MP(x::UnitRange)

"""
    MP(x::StepRangeLen)

Create a (max,+) dense column vector from a given range.

# Examples
```julia-repl
julia> MP(1.0:0.5:3.0)
5-element (max,+) vector:
    1
  1.5
    2
  2.5
    3
```
"""
MP(x::StepRangeLen)

"""
    zero(::Type{MP})

Create the constant (max,+) zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.
"""
Base.zero(::Type{MP})

"""
    zero(::MP)

Create the constant (max,+) zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.

# Examples
```julia-repl
julia> zero(MP)
(max,+) -Inf
```
"""
Base.zero(x::MP)

"""
    one(::Type{MP})

Create the constant (max,+) one (equals to `0`, zero) in classic algebra which
is the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.

# Examples
```julia-repl
julia> one(MP)
(max,+) 0.0
```
"""
Base.one(::Type{MP})

"""
    one(::MP)

Create the constant (max,+) one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.
"""
Base.one(x::MP)

"""
    mpI

Is the equivalent of `LinearAlgebra.I` but for (max,+) type. Allow to create
identity (max,+) matrices.

# Examples
```julia-repl
julia> typeof(mpI)
LinearAlgebra.UniformScaling{MP}

julia> mp_change_display(0)

julia> Matrix(mpI, 2, 2)
2×2 (max,+) dense matrix:
   0.0   -Inf
  -Inf    0.0

julia> mp_change_display(1)

julia> Matrix(mpI, 2, 2)
2×2 (max,+) dense matrix:
  0   .
  .   0

julia> using LinearAlgebra

julia> Matrix(mpI, 2, 2) == Matrix{MP}(I, 2, 2)
true
```
"""
mpI

"""
    +(x::MP, y::MP)

(max,+) operator ⨁. Return the maximum of `x` and `y` as (max,+) type. At
least one parameter shall be a (max,+) number (its conversion to (max,+) is
automatic).

# Examples
```julia-repl
julia> MP(1.0) + MP(3.0)
(max,+) 3

julia> MP(1.0) + 3
(max,+) 3

julia> 1 + MP(3.0)
(max,+) 3

julia> MP(3.0) + -Inf
(max,+) 3
```
"""
Base.:(+)(x::MP, y::MP)

"""
    *(x::MP, y::MP)

(max,+) operator ⨂. Return the sum of `x` and `y` as (max,+) type. At least
one parameter shall be a (max,+) number (its conversion to (max,+) is
automatic).

# Examples
```julia-repl
julia> MP(1.0) * MP(3.0)
(max,+) 4

julia> MP(1.0) * 3
(max,+) 4

julia> 1 * MP(3.0)
(max,+) 4

julia> mp_change_display(0)

julia> MP(1.0) * -Inf
(max,+) -Inf

julia> mp_change_display(1)

julia> MP(1.0) * -Inf
(max,+) .
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
(max,+) -1
```
"""
Base.:(/)(x::MP, y::MP)

"""
    -(x::MP, y::MP)

Minus operator (not used for (max,+) or for (min,+)).
Use the operator / to return the difference between `x` and `y`.

# Examples
```julia-repl
julia> MP(1.0) - MP(2.0)
ERROR: Minus operator does not exist in (max,+) algebra

julia> MP(1.0) / MP(0.0)
(max,+) 1.0

julia> MP(1.0) / mp0
(max,+) Inf
```
"""
Base.:(-)(x::MP, y::MP)

"""
    inv(A::ArrMP)

Return the inverse of a scalar or a (max,+) matrix which is `transpose(-A)` if
the matrix can be inversed. Throw an error if the matrix cannot be inversed.

# Example 1: Scalar
```julia-repl
julia> inv(MP(5))
(max,+) -5

julia> MP(5)^-1
(max,+) -5

julia> MP(5)^0
(max,+) 0
```

# Example 2: Square matrix inversible
```julia-repl
julia> A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
3×3 (max,+) dense matrix:
  .   1   .
  2   .   .
  .   .   3

julia> A^-1
3×3 (max,+) transpose(dense matrix):
   .   -2    .
  -1    .    .
   .    .   -3

julia> inv(A)
3×3 (max,+) dense matrix:
   .   -2    .
  -1    .    .
   .    .   -3

julia> A * inv(A)
3×3 (max,+) dense matrix:
  0   .   .
  .   0   .
  .   .   0

julia> A * inv(A) == inv(A) * A == eye(A)
true
```

# Example 3: Square matrix inversible
```julia-repl
julia> A = [mp0 1 mp0; 2 mp0 mp0]
2×3 (max,+) dense matrix:
  .   1   .
  2   .   .

julia> A^-1
3×2 (max,+) dense matrix:
   .   -2
  -1    .
   .    .

julia> inv(A)
3×2 (max,+) dense matrix:
   .   -2
  -1    .
   .    .

julia> A * inv(A)
2×2 (max,+) dense matrix:
  0   .
  .   0

julia> inv(A) * A
3×3 (max,+) dense matrix:
  0   .   .
  .   0   .
  .   .   .

julia> (A * inv(A)) != (inv(A) * A)
true
```

# Example 4: Non inversible
```julia-repl
julia> A = MP([1 2; 3 4])
2×2 (max,+) dense matrix:
  1   2
  3   4

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
3×3 (max,+) dense matrix:
  .   1   .
  2   .   .
  .   .   3

julia> B = [3 mp0 mp0; mp0 mp0 4; mp0 5 mp0]
3×3 (max,+) dense matrix:
  3   .   .
  .   .   4
  .   5   .

julia> x = A \\ B
3×3 (max,+) dense matrix:
  .   .   2
  2   .   .
  .   2   .

julia> A * x == B
true

julia> A \\ A
3×3 (max,+) dense matrix:
  0   .   .
  .   0   .
  .   .   0
```
"""
Base.:(\)(A::ArrMP, b::ArrMP)

"""
    mp0

Create the constant (max,+) zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mpzero`.

Equivalent to ScicosLab code: `%0` sugar notation for `#(-%inf)`

# Examples
```julia-repl
julia> mp0
(max,+) -Inf

julia> mp0 * 5
(max,+) -Inf

julia> mp0 + 5
(max,+) 5
```
"""
mp0

"""
    ε (\\varepsilon)

Create the constant (max,+) zero (equals to `-∞`, minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mpzero`.

Equivalent to ScicosLab code: `%0` sugar notation for `#(-%inf)`

# Examples
```julia-repl
julia> ε
(max,+) -Inf

julia> ε * 5
(max,+) -Inf

julia> ε + 5
(max,+) 5
```
"""
ε

"""
    mp1

Create the constant (max,+) one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mpone`.

Equivalent to ScicosLab code: `%1` sugar notation for `#(1)`

# Examples
```julia-repl
julia> mp1
(max,+) 0

julia> mp1 * 5
(max,+) 5

julia> mp1 + 5
(max,+) 5
```
"""
mp1

"""
    mpe

Create the constant (max,+) one (equals to `0`, zero) in classic algebra which is
the neutral for operators ⨁ and ⨂. See also `mpone`.

Equivalent to ScicosLab code: `%1` sugar notation for `#(1)`

# Examples
```julia-repl
julia> mpe
(max,+) 0

julia> mpe * 5
(max,+) 5

julia> mpe + 5
(max,+) 5
```
"""
mpe

"""
    mptop

Create the constant (min,+) one (equals to `+∞`, infinity) in classic algebra which is
the neutral for operators ⨁ and ⨂.

Equivalent to ScicosLab code: `%top = #(%inf)`

# Examples
```julia-repl
julia> mptop
(max,+) Inf

julia> mptop * 5
(max,+) Inf

julia> mptop + 5
(max,+) Inf
```
"""
mptop

"""
    mpI

An object of type UniformScaling, representing a (max,+) identity matrix of any size.

# Examples
```julia-repl
julia> mpI
LinearAlgebra.UniformScaling{MP}
e*I

julia> Matrix(mpI, 2, 2)
2×2 (max,+) dense matrix:
  0   .
  .   0
```
"""
mpI

"""
    sparse_map(f, M::SparseMatrixCSC{MP,U})

Map a function ot each element of the (max,+) sparse matrix.

# Examples
```julia-repl
julia> using SparseArrays

julia> sparse_map(x -> x.λ, sparse(MP([1.0 0; 0 1.0])))
2×2 SparseMatrixCSC{Float64, Int64} with 4 stored entries:
 1.0  0.0
 0.0  1.0

julia> sparse_map(x -> MP(x), sparse([1.0 0; 0 1.0]))
2×2 Max-Plus sparse matrix with 2 stored entries:
  [1, 1]  =  1
  [2, 2]  =  1
```
"""
mpsparse_map(f, M::SpaMP)

"""
    plustimes(x::MP)

Convert a (max,+) number to a number in standard algebra. An alternative way
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

Convert a (max,+) dense matrix to a dense matrix in standard algebra.

# Examples
```julia-repl
julia> A = [MP(1.0) 2.0; ε mpe]
2×2 (max,+) dense matrix:
  1   2
  .   0

julia> plustimes(A)
2×2 Matrix{Float64}:
   1.0  2.0
  -Inf  0.0
```
"""
plustimes(A::ArrMP)

"""
    plustimes(A::SpaMP)

Convert a (max,+) sparse matrix to an sparse matrix in standard algebra.

# Examples
```julia-repl
julia> using SparseArrays

julia> S = sparse(eye(MP, 2,2))
2×2 Max-Plus sparse matrix with 2 stored entries:
  [1, 1]  =  0
  [2, 2]  =  0

julia> findnz(S)
([1, 2], [1, 2], MP[0, 0])

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

Convert a sparse (max,+) array to a dense (max,+) array.
Alternative function name: dense.

# Examples
```julia-repl
julia> using SparseArrays

julia> mp_change_display(0)

julia> full(spzeros(MP, 2,5))
2×5 (max,+) dense matrix:
  -Inf   -Inf   -Inf   -Inf   -Inf
  -Inf   -Inf   -Inf   -Inf   -Inf

julia> mp_change_display(1)

julia> full(spzeros(MP, 2,5))
2×5 (max,+) dense array:
  .   .   .   .   .
  .   .   .   .   .
```
"""
full(S::SpaMP)

"""
    dense(::SpaMP})

Convert a sparse (max,+) array to a dense (max,+) array.
Alternative function name: full.

# Examples
```julia-repl
julia> using SparseArrays

julia> mp_change_display(0)

julia> dense(spzeros(MP, 2,5))
2×5 (max,+) dense array:
 -Inf  -Inf  -Inf  -Inf  -Inf
 -Inf  -Inf  -Inf  -Inf  -Inf

julia> mp_change_display(1)

julia> dense(spzeros(MP, 2,5))
2×5 (max,+) dense matrix:
  .   .   .   .   .
  .   .   .   .   .
```
"""
dense(S::SpaMP)

"""
    array(::SpaMP})

Convert a sparse (max,+) array to a dense non (max,+) array.

# Examples
```julia-repl
julia> using SparseArrays

julia> array(spzeros(MP, 2,2))
2×2 Matrix{Float64}:
 -Inf  -Inf
 -Inf  -Inf
```
"""
array(S::SpaMP)

"""
    min(x::MP, y::MP)

Return the minimun of (max,+) numbers `x` and `y`.

min(x, y) = (x ⨂ y) ⨸ (x ⨁ y)

# Examples
```julia-repl
julia> min(MP(1), 3)
(max,+) 1

julia> min(MP([10 1; 10 1]), MP([4 5; 6 5]))
2×2 (max,+) dense matrix:
  4   1
  6   1
```
"""
Base.:min(x::MP, y::MP)

"""
    eye(MP, n::Int64)

Construct a (max,+) identity dense n-by-n matrix.

# Examples
```julia-repl
julia> eye(MP, 2)
2×2 (max,+) dense matrix:
   0.0   -Inf
  -Inf    0.0
```
"""
eye(MP, n::Int64)

"""
    eye(MP, m::Int64, n::Int64)

Construct a (max,+) identity dense m-by-n matrix.

# Examples
```julia-repl
julia> eye(MP, 2, 3)
2×3 (max,+) dense matrix:
   0.0   -Inf   -Inf
  -Inf    0.0   -Inf
```
"""
eye(MP, m::Int64, n::Int64)

"""
    eye(MP, A::Array{T})

Construct a (max,+) identity dense matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> eye(MP, A)
2×2 (max,+) dense matrix:
  0.0  -Inf
 -Inf   0.0

julia> eye(MP, MP(A))
2×2 (max,+) dense matrix:
  0.0  -Inf
 -Inf   0.0
```
"""
eye(MP, A::Array)

"""
    spzeros(MP, n::Int64)

Construct a (max,+) zero n-by-m sparse matrix.

# Examples
```julia-repl
julia> using SparseArrays

julia> spzeros(MP, 2)
2-element SparseVector{MP,Int64} with 0 stored entries
```
"""
spzeros(MP, n::Int64)

"""
    spzeros(MP, n::Int64)

Construct a sparse (max,+) zero m-by-n matrix.

# Examples
```julia-repl
julia> using SparseArrays

julia> spzeros(MP, 2,5)
2×5 SparseMatrixCSC{MP,Int64} with 0 stored entries
  .   .   .   .   .
  .   .   .   .   .

julia> full(spzeros(MP, 2,5))
2×5 (max,+) dense matrix:
  -Inf   -Inf   -Inf   -Inf   -Inf
  -Inf   -Inf   -Inf   -Inf   -Inf
```
"""
spzeros(MP, m::Int64, n::Int64)

"""
    spzeros(MP, A::Array{T})

Construct a sparse  (max,+) zero matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> using SparseArrays

julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> spzeros(MP, A)
2×2 (max,+) sparse matrix with 0 stored entries:
  .   .
  .   .

julia> spzeros(MP, MP(A))
2×2 (max,+) sparse matrix with 0 stored entries:
  .   .
  .   .
```
"""
spzeros(MP, A::Array)

"""
    ones(MP, n::Int64)

Construct a (max,+) one n-by-1 matrix.

# Examples
```julia-repl
julia> ones(MP, 2)
2-element (max,+) vector:
  0.0
  0.0
```
"""
ones(MP, n::Int64)

"""
    ones(MP, m::Int64, n::Int64)

Construct a (max,+) one m-by-n matrix.

# Examples
```julia-repl
julia> ones(MP, 3,2)
3×2 (max,+) dense matrix:
  0.0   0.0
  0.0   0.0
  0.0   0.0
```
"""
ones(MP, m::Int64, n::Int64)

"""
    ones(MP, A::Array{T})

Construct a sparse (max,+) ones matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> ones(MP, A)
2×2 (max,+) dense matrix:
  0.0   0.0
  0.0   0.0
```
"""
ones(MP, A::Array)

"""
    tr(A::Array)

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> tr([MP(1) 2; 3 4])
(max,+) 4.0
```
"""
tr(A::ArrMP)

"""
    norm(A)

Compute the norm of the full or sparce matrix A.
Return the largest entry minus smallest entry of A.

# Examples
```
julia-repl
julia> using SparseArrays

julia> A = MP([1 20 2;30 400 4;4 50 10])
3×3 (max,+) dense matrix:
   1.0    20.0    2.0
  30.0   400.0    4.0
   4.0    50.0   10.0

julia> S = MP(sparse(A))
3×3 (max,+) sparse matrix with 9 stored entries:
   1    20    2
  30   400    4
   4    50   10

julia> mpnorm(A)
(max,+) 399

julia> mpnorm(S)
(max,+) 399
```
"""
norm(A::Array)

"""
    B = star(A::ArrMP)

Solve `x = Ax + I` in the (max,+) algebra. When there is no circuits with
positive weight in G(A) (the incidence graph of A) `B = I + A + ... + A^(n-1)`
where n denotes the order of the square matrix A.

See also plus.

# Arguments
- A : (max,+) full square matrix.
- B : Id + A + A^2 + ...

# Examples
```julia-repl
julia> star(MP(1.0))
(max,+) Inf

julia> star(MP(-1.0))
(max,+) 0

julia> star(MP([1.0 2; 3 4]))
2×2 (max,+) dense matrix:
  Inf   Inf
  Inf   Inf

julia> A = star(MP([-3.0 -2;-1 0]))
2×2 (max,+) dense matrix:
   0   -2
  -1    0

julia> B = star(A)
2×2 (max,+) dense matrix:
   0   -2
  -1    0

julia> B == (eye(MP, 2,2) + A)
true

julia> B == (B + A * A)
true
```
"""
star(A::ArrMP)

"""
    star(x::MP)

Make x a 1x1 matrix then call star(A::ArrMP).
See star(A::ArrMP) for more information.
"""
star(x::MP)

"""
    astarb(A::ArrMP, b::ArrMP)

(max,+) linear system solution.

Solve `x = Ax + b` in the (max,+) algebra when there is no circuits with
positive weight in `G(A')` (the incidence graph of `A'`, that is it exists an
arc from `j` to `i` if `A_ij` is nonzero).

TODO It is much more efficient in time and memory than `star(A) * b`.
"""
astarb(A::ArrMP, b::ArrMP)

"""
    B = plus(A::ArrMP)

Compute `A * A^* = A + A^2 + ...` of a maxplus matrix A.
See also star.

# Arguments
- A : (max,+) full square matrix.
- B : Id + A + A^2 + ...

# Examples
```julia-repl
julia> plus(MP(1.0))
(max,+) Inf

FIXME
julia> A = MP([-3.0 -2; -1 0])
2×2 (max,+) dense matrix:
  -3   -2
  -1    0

julia> B = plus(A)
2×2 (max,+) dense matrix:
   0   -2
  -1    0

julia> B == (A * star(A))
true
```
"""
plus(A::ArrMP)

"""
    plus(x::MP)

Make x a 1x1 matrix then call plus(A::ArrMP).
See plus(A::ArrMP) for more information.

# Examples
```julia-repl
FIXME
julia> plus(MP(-1.0))
(max,+) -1
```
"""
plus(x::MP)

"""
    mp_change_display(style::Int)

Change the style of behavior of functions `Base.show()`:
- `-Inf` are displayed either with `ε` (style 2 or 3) or `.` symbols (style 1).
- `0` are displayed either with `e` (style 3) or '0' symbols (style 1 or 2).
- else: `-Inf` and `0` are displayed in Julia default sytle (style 0).

If this function is not called, by default the ScicosLab style will be used
(style 1).

# Examples
```julia-repl
julia> mp_change_display(0)
I will show -Inf and 0.0

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
   0.0   -Inf   -Inf
  -Inf    0.0   -Inf
  -Inf   -Inf    0.0

julia> mp_change_display(1)
I will show -Inf as .

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  0   .   .
  .   0   .
  .   .   0

julia> mp_change_display(2)
I will show -Inf as . and 0.0 as e

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  e   .   .
  .   e   .
  .   .   e

julia> mp_change_display(3)
I will show -Inf as ε

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  0   ε   ε
  ε   0   ε
  ε   ε   0

julia> mp_change_display(4)
I will show -Inf as ε and 0.0 as e

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  e   ε   ε
  ε   e   ε
  ε   ε   e
```
"""
mp_change_display(style::Int)

"""
    Base.show(io::IO, ::MIME"text/plain", A::ArrTrop)

Display a tropical array on the desired output (i.e. console).
Controled by mp_change_display().

# Examples
```julia-repl

julia> mp_change_display(4)
I will show -Inf as ε and 0.0 as e

julia> show(stdout, "text/plain", MP([1 0; -Inf 6]))
2×2 (max,+) dense matrix:
  1   e
  ε   6

julia> mp_change_display(0)
I will show -Inf and 0.0

julia> show(stdout, "text/plain", MP([1 0; -Inf 6]))
2×2 (max,+) dense matrix:
  1   0
  .   6
```
"""
Base.show(io::IO, ::MIME"text/plain", A::ArrTrop)

"""
    Base.show(io::IO, ::MIME"text/latex", A::ArrTrop)

Generate the \\LaTeX formula from the given tropical array.
Controled by mp_change_display().

# Examples
```julia-repl
julia> mp_change_display(4)
I will show -Inf as ε and 0.0 as e

julia> show(stdout, "text/latex", MP([1 0; -Inf 6]))
\\left[
\\begin{array}{*{20}c}
1 & e \\\\
\\varepsilon & 6 \\\\
\\end{array}
\\right]

julia> mp_change_display(0)
I will show -Inf and 0.0

julia> show(stdout, "text/latex", MP([1 0; -Inf 6]))
\\left[
\\begin{array}{*{20}c}
1 & 0 \\\\
-\\infty & 6 \\\\
\\end{array}
\\right]
```
"""
Base.show(io::IO, ::MIME"text/latex", A::ArrTrop)
