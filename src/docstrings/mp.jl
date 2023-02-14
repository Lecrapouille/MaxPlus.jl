# ==============================================================================
# Docstring for Max-Plus functions
# ==============================================================================

################################################################################
###
### Max-Plus constructors
###
################################################################################

# ==============================================================================
"""
    MP(x::Float64)
    MP(x::Int64)
    MP(x::MP)
    MP(x::MI)

Immutable Julia structure for creating a (max,+) scalar. This promotes the given
number (Float64 or Int64 or MaxPlus or MinPlus) to a number in the tropical
semi-ring (max, +) (ℝ ∪ {-∞}, ⨁, ⨂) where ℝ is the domain of reals, ⨁ is the
usual multiplication and ⨂ is the usual maximum.

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

julia> MP(MI(3.5))
(max,+) 3.5
```
"""
MP(::Float64)

# ==============================================================================
"""
    MP(b::Bool)

Immutable Julia structure for promoting the given Boolean value to an neutral
number in the tropical semi-ring (max, +). This constructor does not make sense
mathematically speaking but it is needed by Julia for defining the identity
operator (max, +) `LinearAlgebra.I`.  Without this method, Julia will not be
able to create correctly (max,+) identity matrices.

# Examples
```julia-repl
julia> set_tropical_display(0)

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
```
"""
MP(::Bool)

# ==============================================================================
"""
    MP(A::Array{Float64})
    MP(A::Array{Int64})
    MP(A::Array{MP})
    MP(A::Array{MI})

Promote the given dense array of elements either from classic algebra (Float64
or Int64) or array of Max-Plus or Min-Plus elements to a dense array of tropical
semi-ring (max, +) numbers.

Note: if the array already contains at least one (max,+) element then the `MP()`
is not necessary since (max,+) numbers contaminate other Float64 or Int64
numbers.

# Example where the MP() constructor is needed
```julia-repl
julia> A = MP([1.0 -Inf; 0.0 4])
2×2 (max,+) dense matrix:
  1   .
  0   4

julia> typeof(A)
Matrix{MP} (alias for Array{Trop{MaxPlus.Max}, 2}
```

# Example where the MP() constructor is not needed
```julia-repl
julia> A = [MP(1.0) -Inf; 0.0 4]
2×2 (max,+) dense matrix:
  1   .
  0   4

julia> typeof(A)
Matrix{MP} (alias for Array{Trop{MaxPlus.Max}, 2})

julia> MP(MI([1.0 -Inf; 0.0 4]))
2×2 (max,+) dense matrix:
  1.0   -Inf
  0.0    4.0
```
"""
MP(::Array)

# ==============================================================================
"""
    MP(S::SparseMatrixCSC{Float64})
    MP(S::SparseMatrixCSC{Int64})
    MP(S::SparseMatrixCSC{MP})
    MP(S::SparseMatrixCSC{MI})

Promote the given sparse matrix of elements from classic algebra (Float64 or
Int64) to a sparse array of tropical semi-ring (max, +) numbers. Also accept to
promote elements from tropical semi-ring (max, +). By default, and following
Julia rules, explicit (max,+) zeros (`ε`, `mp0`, `MP(-Inf)`) are not
removed. You shall call dropzeros() to remove them.

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
MP(::SparseMatrixCSC)

# ==============================================================================
"""
    MP(I::AbstractVector, J::AbstractVector, V::AbstractVector{Float64})
    MP(I::AbstractVector, J::AbstractVector, V::AbstractVector{Int64})
    MP(I::AbstractVector, J::AbstractVector, V::AbstractVector{MP})
    MP(I::AbstractVector, J::AbstractVector, V::AbstractVector{MI})

Construct a sparse CSC matrix of tropical semi-ring (max, +) elements from
classic algebra (Float64 or Int64) to a sparse array of tropical semi-ring (max,
+) numbers. Also accept to promote elements from tropical semi-ring (max,
+). Note that zeros elements (`ε`, `mp0`, `MP(-Inf)`) are removed.

Elements are stored as: S[I[k], J[k]] = V[k] where elements of V are typed of
Float64 or Int64 or MP.

# Examples
```julia-repl
julia> A = MP([1; 2; 3], [1; 2; 3], [42.5; -Inf; 44])
3×3 (max,+) sparse matrix with 2 stored entries:
  [1, 1]  =  42.5
  [3, 3]  =  44

julia> A = MP([1; 2; 3], [1; 2; 3], MP([42.5; -Inf; 44]))
3×3 (max,+) sparse matrix with 2 stored entries:
  [1, 1]  =  42.5
  [3, 3]  =  44
```
"""
MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)

# ==============================================================================
"""
    MP(V::SparseVector{Float64})
    MP(V::SparseVector{Int64})
    MP(V::SparseVector{MP})
    MP(V::SparseVector{MI})

Convert a sparse vector of tropical semi-ring (max, +) elements from given
classic algebra to a (max,+) sparse vector. By default, explicit (max,+) zeros
(`ε`, `mp0`, `MP(-Inf)`) are removed.

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
MP(::SparseVector)

# ==============================================================================
"""
    MP(I::AbstractVector, V::AbstractVector{Float64})
    MP(I::AbstractVector, V::AbstractVector{Int64})
    MP(I::AbstractVector, V::AbstractVector{MP})
    MP(I::AbstractVector, V::AbstractVector{MI})

Construct a sparse (max,+) vector such as S[I[k]] = V[k]. Explicit (max,+)
zeros (`ε`, `mp0`, `MP(-Inf)`) are removed.

# Examples
```julia-repl
julia> A = MP([1; 2; 3], [42.5; -Inf; 44])
3-element (max,+) sparse vector with 2 stored entries:
  [1]  =  42.5
  [3]  =  44
```
"""
MP(I::AbstractVector, V::AbstractVector)

# ==============================================================================
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
MP(::UnitRange)

# ==============================================================================
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

################################################################################
###
### Max-Plus scalars
###
################################################################################

# ==============================================================================
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

julia> x = MP(1)
(max,+) 1.0

julia> x.λ
1.0

julia> typeof(x.λ)
Float64
```
"""
plustimes(n::MP)

# ==============================================================================
"""
    star(x::MP)

Make x a 1x1 matrix then call star(A::Array{MP}).
See star(A::Array{MP}) for more information.
"""
star(x::MP)

# ==============================================================================
"""
    plus(x::MP)

Make x a 1x1 matrix then call plus(A::Array{MP}).
See plus(A::Array{MP}) for more information.

# Examples
```julia-repl
FIXME
julia> plus(MP(-1.0))
(max,+) -1
```
"""
plus(x::MP)

################################################################################
###
### Max-Plus Algebra
###
################################################################################

# ==============================================================================
"""
    zero(::Type{MP})

Create the constant (max,+) zero equals to `-∞` (minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.
"""
Base.zero(::Type{MP})

# ==============================================================================
"""
    zero(::MP)

Create the constant (max,+) zero equals to `-∞` (minus infinity) in classic
algebra which is the neutral for the ⨁ operator. See also `mp0` and `ε`.

# Examples
```julia-repl
julia> zero(MP)
(max,+) -Inf
```
"""
Base.zero(x::MP)

# ==============================================================================
"""
    one(::Type{MP})

Create the constant (max,+) one equals to `0.0` (zero) in classic algebra which
is the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.

# Examples
```julia-repl
julia> one(MP)
(max,+) 0.0
```
"""
Base.one(::Type{MP})

# ==============================================================================
"""
    one(::MP)

Create the constant (max,+) one equals to `0.0` (zero) in classic algebra which
is the neutral for operators ⨁ and ⨂. See also `mp1` and `mpe`.
"""
Base.one(x::MP)

# ==============================================================================
"""
    +(x::MP, y::MP)

(max,+) operator ⨁. Return the maximum of `x` and `y` as (max,+) type. At least
one parameter shall be a (max,+) number (the other can be typed of Float64 or
Int64 and its conversion to a (max,+) number is automatic).

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

# ==============================================================================
"""
    *(x::MP, y::MP)

(max,+) operator ⨂. Return the sum of `x` and `y` as (max,+) type. At least one
parameter shall be a (max,+) number (the other can be typed of Float64 or Int64
and its conversion to a (max,+) number is automatic).

# Examples
```julia-repl
julia> MP(1.0) * MP(3.0)
(max,+) 4

julia> MP(1.0) * 3
(max,+) 4

julia> 1 * MP(3.0)
(max,+) 4

julia> set_tropical_display(0)

julia> MP(1.0) * -Inf
(max,+) -Inf

julia> set_tropical_display(1)

julia> MP(1.0) * -Inf
(max,+) .
```
"""
Base.:(*)(x::MP, y::MP)

# ==============================================================================
"""
    Base.:(^)(::MP, ::Number)

In (max,+) algebra the power operator behaves like a multiplication in classical
algebra.

# Examples
```julia-repl
julia> MP(2)^5
(max,+) 10

julia> MP(2)^-1
(max,+) -2
```
"""
Base.:(^)(::MP, ::Number)

# ==============================================================================
"""
    /(x::MP, y::MP)

Divisor operator. Return the difference between `x` and `y` in classic algebra.

x / y ≜ (x ⨂ y^-1)^-1

# Examples
```julia-repl
julia> MP(1.0) / MP(2.0)
(max,+) -1
```
"""
Base.:(/)(x::MP, y::MP)

# ==============================================================================
"""
    Base.:(\\)(x::MP, y::MP)

Divisor operator. Return the difference between `y` and `x` in classic algebra.

x \\ y ≜ (y^-1 ⨂ x)^-1

# Examples
```julia-repl

julia> MP(2) \\ MP(3)
(max,+) 1
```
"""
Base.:(\)(x::MP, y::MP)

# ==============================================================================
"""
    min(x::MP, y::MP)

Return the minimun of (max,+) numbers `x` and `y`.

min(x, y) ≜ (x ⨂ y) ⨸ (x ⨁ y)

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

# ==============================================================================
"""
    -(x::MP, y::MP)

The minus operator is not used in (max,+) algebra. Calling this operator will
throw an error. Use the operator / to return the difference between `x` and `y`.

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

################################################################################
###
### Max-Plus constants
###
################################################################################

# ==============================================================================
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

# ==============================================================================
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
MaxPlus.ε

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
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

################################################################################
###
### Max-Plus Matrix (TODO speye)
###
################################################################################

# ==============================================================================
"""
    sparse_map(f, S::SparseMatrixCSC{MP})

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
sparse_map(f, S::SparseMatrixCSC{MP})

# ==============================================================================
"""
    mpI

An object of type UniformScaling, representing a (max,+) identity matrix of any size.
Is the equivalent of `LinearAlgebra.I` but for (max,+) type. Allow to create
identity (max,+) matrices.

# Examples
```julia-repl
julia> typeof(mpI)
LinearAlgebra.UniformScaling{MP}

julia> set_tropical_display(0)

julia> Matrix(mpI, 2, 2)
2×2 (max,+) dense matrix:
   0.0   -Inf
  -Inf    0.0

julia> set_tropical_display(1)

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

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
"""
    speye(MP, n::Int64)

Construct a (max,+) identity sparse n-by-n matrix.

# Examples
```julia-repl
julia> speye(MP, 2)
2×2 (max,+) sparse matrix with 2 stored entries:
  [1, 1]  =  0.0
  [2, 2]  =  0.0
```
"""
speye(MP, n::Int64)

# ==============================================================================
"""
    speye(MP, m::Int64, n::Int64)

Construct a (max,+) identity sparse m-by-n matrix.

# Examples
```julia-repl
julia> speye(MP, 2, 3)
2×3 (max,+) sparse matrix with 2 stored entries:
  [1, 1]  =  0.0
  [2, 2]  =  0.0
```
"""
speye(MP, m::Int64, n::Int64)

# ==============================================================================
"""
    zeros(MP, n::Int64)

Construct a (max,+) zero n-by-1 (max,+) dense matrix.

# Examples
```julia-repl
julia> zeros(MP, 2)
2-element (max,+) vector:
  0.0
  0.0
```
"""
zeros(MP, n::Int64)

# ==============================================================================
"""
    zeros(MP, m::Int64, n::Int64)

Construct a (max,+) zero m-by-n (max,+) dense matrix.

# Examples
```julia-repl
julia> zeros(MP, 3,2)
3×2 (max,+) dense matrix:
  0.0   0.0
  0.0   0.0
  0.0   0.0
```
"""
zeros(MP, m::Int64, n::Int64)

# ==============================================================================
"""
    zeros(MP, A::Array{T})

Construct a (max,+) zero dense matrix of same dimension
and of the same type that the matrix given as parameter.

# Examples
```julia-repl
julia> A=[1.0 2; 3 4]
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> zeros(MP, A)
2×2 (max,+) dense matrix:
  0.0   0.0
  0.0   0.0
```
"""
zeros(MP, A::Array)

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
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

# ==============================================================================
"""
    full(::SparseMatrixCSC{MP}})

Convert a sparse (max,+) array to a dense (max,+) array.
Alternative function name: dense.

# Examples
```julia-repl
julia> using SparseArrays

julia> set_tropical_display(0)

julia> full(spzeros(MP, 2,5))
2×5 (max,+) dense matrix:
  -Inf   -Inf   -Inf   -Inf   -Inf
  -Inf   -Inf   -Inf   -Inf   -Inf

julia> set_tropical_display(1)

julia> full(spzeros(MP, 2,5))
2×5 (max,+) dense array:
  .   .   .   .   .
  .   .   .   .   .
```
"""
full(S::SparseMatrixCSC{MP})

# ==============================================================================
"""
    dense(::SparseMatrixCSC{MP}})

Convert a sparse (max,+) array to a dense (max,+) array.
Alternative function name: full.

# Examples
```julia-repl
julia> using SparseArrays

julia> set_tropical_display(0)

julia> dense(spzeros(MP, 2,5))
2×5 (max,+) dense array:
 -Inf  -Inf  -Inf  -Inf  -Inf
 -Inf  -Inf  -Inf  -Inf  -Inf

julia> set_tropical_display(1)

julia> dense(spzeros(MP, 2,5))
2×5 (max,+) dense matrix:
  .   .   .   .   .
  .   .   .   .   .
```
"""
dense(S::SparseMatrixCSC{MP})

# ==============================================================================
"""
    plustimes(A::Array{MP})

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
plustimes(A::Array{MP})

# ==============================================================================
"""
    plustimes(A::SparseMatrixCSC{MP})

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
plustimes(S::SparseMatrixCSC{MP})

# ==============================================================================
"""
    inv(A::Array{MP})

Return the inverse of a scalar or a (max,+) matrix which is `transpose(-A)`.

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
2×2 (max,+) dense matrix:
  -1   -3
  -2   -4

julia> (A * inv(A)) != eye(A)
true

julia> (inv(A) * A) != eye(A)
true
```
"""
Base.inv(A::Array{MP})

# ==============================================================================
"""
    \\(A::AbstractMatrix{MP}, b::AbstractMatrix{MP})

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
Base.:(\)(A::AbstractMatrix{MP}, b::AbstractMatrix{MP})

# ==============================================================================
"""
    Base.:(\\)(A::AbstractMatrix{MP}, b::MP)

TODO

# Examples
```julia-repl
julia>
```
"""
Base.:(\)(A::AbstractMatrix{MP}, b::MP)

# ==============================================================================
"""
    Base.:(\\)(A::MP, b::AbstractMatrix{MP})

TODO

# Examples
```julia-repl
julia>
```
"""
Base.:(\)(A::MP, b::AbstractMatrix{MP})

# ==============================================================================
"""
    Base.:(/)(A::AbstractMatrix{MP}, B::AbstractMatrix{MP})

TODO

# Examples
```julia-repl
julia>
```
"""
Base.:(/)(A::AbstractMatrix{MP}, B::AbstractMatrix{MP})

# ==============================================================================
"""
    Base.:(/)(A::AbstractMatrix{MP}, b::MP)

TODO

# Examples
```julia-repl
julia>
```
"""
Base.:(/)(A::AbstractMatrix{MP}, b::MP)

# ==============================================================================
"""
    Base.:(/)(a::MP, b::AbstractMatrix{MP})

TODO

# Examples
```julia-repl
julia>
```
"""
Base.:(/)(a::MP, b::AbstractMatrix{MP})

# ==============================================================================
"""
    tr(A::Array{MP})

Compute the trace of the matrix (summation of diagonal elements).

# Examples
```julia-repl
julia> tr([MP(1) 2; 3 4])
(max,+) 4.0
```
"""
tr(A::Array{MP})

# ==============================================================================
"""
    norm(A::Array{MP})

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

# ==============================================================================
"""
    astarb(A::Array{MP}, b::Array{MP})

(max,+) linear system solution.

Solve `x = Ax + b` in the (max,+) algebra when there is no circuits with
positive weight in `G(A')` (the incidence graph of `A'`, that is it exists an
arc from `j` to `i` if `A_ij` is nonzero).

TODO It is much more efficient in time and memory than `star(A) * b`.
"""
astarb(A::Array{MP}, b::Array{MP})

# ==============================================================================
"""
    B = star(A::Array{MP})

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
star(A::Array{MP})

# ==============================================================================
"""
    B = plus(A::Array{MP})

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
plus(A::Array{MP})

################################################################################
###
### Max-Plus Spectral (TODO: semi_howard)
###
################################################################################

# ==============================================================================
"""
    [λ,v,p,c,n] = howard(S::SparseMatrixCSC{MP})

(max,+) mpeigenvalues mpeigenvectors (Howard algorithm).

Maxplus right mpeigenvalues and mpeigenvectors of a full or sparse maxplus matrix by
Howard algorithm. The mpeigenvalues are considered as the average cost per unit of
time for the corresponding dynamic programming problem.

The values taken by the entries of `λ` are the mpeigenvalues. If `S` is
irreducible, `λ` is constant, it is the mpeigenvalue and `v` is a corresponding
mpeigenvector (in this case, there exits only one mpeigenvalue but more than one
mpeigenvectors may exist).

Otherwise, `S` can be decomposed into irreducible components (in a certain
numbering of rows and columns, it becomes block-triangular with diagonal
irreducible blocks), `λ` is constant over each component and this constant is
the mpeigenvalue, the corresponding entries of `v`, completed by -inf for the
other blocks, provide a corresponding mpeigenvector.

`p` gives an optimal policy which satisfies \$S\\_{i,p(i)} v\\_{p(i)} = λ + v\\_i\$

# Remark:

- For the block triangular case, take a look at the examples to see what happen
  precisely on the transient block. All the mpeigen values are not found and the
  support of the mpeigenvectors depends of the mpeigenvalues of the blocks.

- For the block diagonal case all the mpeigen values are found and the support of
  the mpeigenvectors are clear.

# Outputs:

- `λ`: mpeigenvalues
- `v`: mpeigenvectors
- `p`: optimal policy (indices of the saturating entries of S)
- `c`: number of connected components of the optimal policy
- `n`: number of iterations of Howard algorithm

# Example:
```julia-repl
julia> using SparseArrays

julia> A = MP([1 2; 3 4])
2×2 (max,+) dense matrix:
  1   2
  3   4

julia> r = howard(sparse(A))
MaxPlus.HowardResult(MP[4, 4], MP[2, 4], [2, 2], 1, 1)
```
"""
howard(S::SparseMatrixCSC{MP}, max_iterations::Int64 = 1000)

# ==============================================================================
"""
    λ,v = mpeigen(S::SparseMatrixCSC{MP})
    λ,v = mpeigen(S::Matrix{MP})

(max,+) mpeigenvalues mpeigenvectors (interface the Howard algorithm):
- `λ`: mpeigenvalues
- `v`: mpeigenvectors

# Example:
```julia-repl
julia> using SparseArrays, LinearAlgebra

julia> A = MP([1 2; 3 4])
2×2 (max,+) dense matrix:
  1   2
  3   4

julia> λ,v = mpeigen(A)
(MP[4, 4], MP[2, 4])

# λ is constant
julia> (A * v) == (λ[1] * v)
true
```

# Example 2: Two blocks diagonal matrix.
```julia-repl
julia> using SparseArrays, LinearAlgebra

julia> S = sparse([mp0 2 mp0; mp1 mp0 mp0; mp0 mp0 2])
3×3 (max,+) sparse matrix with 3 stored entries:
  [2, 1]  =  0
  [1, 2]  =  2
  [3, 3]  =  2

julia> λ,v = mpeigen(S)
(MP[1, 1, 2], MP[1, 0, 2])

# The entries of λ take two values
julia> (S / diagm(λ)) * v == v
true
```

# Example 3: Block triangular matrix with 2 mpeigenvalues.
```julia-repl
julia> S = sparse([1 1; mp0 2])
2×2 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  1
  [1, 2]  =  1
  [2, 2]  =  2

julia> λ,v = mpeigen(S)
(MP[2, 2], MP[1, 2])

julia> (S * v) == (λ[1] * v)
true

# But MP(1) is also mpeigen value
julia> S * [0; mp0] == MP(1) * [0; mp0]
true
```

# Example 4: Block triangular matrix with 1 mpeigenvalue
```julia-repl
julia> using SparseArrays, LinearAlgebra

julia> S = sparse([2 1; mp0 mp1])
2×2 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  2
  [1, 2]  =  1
  [2, 2]  =  0

julia> λ,v = mpeigen(S)
(MP[2, 0], MP[2, 0])

# λ is not constant λ[1] is mpeigen value
# with mpeigen vector [v(1);0]
julia> (S / diagm(λ)) * v == v
true
```
"""
mpeigen(S::SparseMatrixCSC{MP})

################################################################################
###
### Max-Plus Display
###
################################################################################

# ==============================================================================
"""
    set_tropical_display(style::Int)

Change the style of behavior of functions `Base.show()`:
- `-Inf` are displayed either with `ε` (style 2 or 3) or `.` symbols (style 1).
- `0` are displayed either with `e` (style 3) or '0' symbols (style 1 or 2).
- else: `-Inf` and `0` are displayed in Julia default sytle (style 0).

If this function is not called, by default the ScicosLab style will be used
(style 1).

# Examples
```julia-repl
julia> set_tropical_display(0)
I will show -Inf and 0.0

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
   0.0   -Inf   -Inf
  -Inf    0.0   -Inf
  -Inf   -Inf    0.0

julia> set_tropical_display(1)
I will show -Inf as .

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  0   .   .
  .   0   .
  .   .   0

julia> set_tropical_display(2)
I will show -Inf as . and 0.0 as e

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  e   .   .
  .   e   .
  .   .   e

julia> set_tropical_display(3)
I will show -Inf as ε

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  0   ε   ε
  ε   0   ε
  ε   ε   0

julia> set_tropical_display(4)
I will show -Inf as ε and 0.0 as e

julia> eye(MP, 3, 3)
3×3 (max,+) dense matrix:
  e   ε   ε
  ε   e   ε
  ε   ε   e
```
"""
set_tropical_display(style::Int)

# ==============================================================================
"""
    Base.show(io::IO, ::MIME"text/plain", A::Matrix{MP})

Display a tropical array on the desired output (i.e. console).
Controled by set_tropical_display().

# Examples
```julia-repl

julia> set_tropical_display(4)
I will show -Inf as ε and 0.0 as e

julia> show(stdout, "text/plain", MP([1 0; -Inf 6]))
2×2 (max,+) dense matrix:
  1   e
  ε   6

julia> set_tropical_display(0)
I will show -Inf and 0.0

julia> show(stdout, "text/plain", MP([1 0; -Inf 6]))
2×2 (max,+) dense matrix:
  1   0
  .   6
```
"""
Base.show(io::IO, ::MIME"text/plain", A::Matrix{MP})

# ==============================================================================
"""
    Base.show(io::IO, ::MIME"text/latex", A::Matrix{MP})

Generate the \\LaTeX formula from the given tropical array.
Controled by set_tropical_display().

# Examples
```julia-repl
julia> set_tropical_display(4)
I will show -Inf as ε and 0.0 as e

julia> show(stdout, "text/latex", MP([1 0; -Inf 6]))
\\left[
\\begin{array}{*{20}c}
1 & e \\\\
\\varepsilon & 6 \\\\
\\end{array}
\\right]

julia> set_tropical_display(0)
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
Base.show(io::IO, ::MIME"text/latex", A::Matrix{MP})

# ==============================================================================
"""
    show(io::IO, x::MP)

Display a Max-Plus number depending on the currently set style:

# Examples
```julia-repl
julia> set_tropical_display(0); eye(MP, 2,2)
2×2 (max,+) dense matrix:
   0.0   -Inf
  -Inf    0.0

julia> set_tropical_display(1); eye(MP, 2,2)
2×2 (max,+) dense matrix:
  0   .
  .   0

julia> set_tropical_display(2); eye(MP, 2,2)
2×2 (max,+) dense matrix:
  e   .
  .   e

julia> set_tropical_display(3); eye(MP, 2,2)
2×2 (max,+) dense matrix:
  0   ε
  ε   0

julia> set_tropical_display(4); eye(MP, 2,2)
2×2 (max,+) dense matrix:
  e   ε
  ε   e
```
"""

# ==============================================================================
"""
    LaTeX(io::IO, A::Array{MP})

Base function for convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
neutral and absorbing elements depends on set_tropical_display(style).

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
LaTeX(io::IO, A::Matrix{MP})
