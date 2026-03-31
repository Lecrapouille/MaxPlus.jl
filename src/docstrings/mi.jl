# ==============================================================================
# Docstring for Min-Plus functions
# ==============================================================================

################################################################################
###
### Min-Plus constructors
###
################################################################################

# ==============================================================================
"""
    MI(x::Float64)
    MI(x::Int64)
    MI(x::MP)
    MI(x::MI)

Immutable Julia structure for creating a (min,+) scalar. This promotes the given
number (Float64 or Int64 or MaxPlus or MinPlus) to a number in the tropical
semi-ring (min, +) (ℝ ∪ {+∞}, ⨁, ⨂) where ℝ is the domain of reals, ⨁ is the
usual minimum and ⨂ is the usual addition.

# Notes
Similar to (max,+) but with minimum instead of maximum.

# Examples
```julia-repl
julia> a = MI(3.5)
(min,+) 3.5

julia> typeof(a)
MI (alias for Tropical{MaxPlus.Min})

julia> a = MI(3)
(min,+) 3

julia> typeof(a)
MI (alias for Tropical{MaxPlus.Min})

julia> MI(MP(3.5))
(min,+) 3.5
```
"""
MI(::Float64)

# ==============================================================================
"""
    MI(b::Bool)

Immutable Julia structure for promoting the given Boolean value to an neutral
number in the tropical semi-ring (min, +). This constructor is needed by Julia
for defining the identity operator (min, +) `LinearAlgebra.I`.

# Examples
```julia-repl
julia> set_tropical_display(0)

julia> MI(true)
(min,+) 0.0

julia> MI(false)
(min,+) +Inf

julia> using LinearAlgebra

julia> Matrix{MI}(I, 2, 2)
2×2 (min,+) dense matrix:
   0.0   +Inf
  +Inf    0.0
```
"""
MI(::Bool)

# ==============================================================================
"""
    MI(A::Array{Float64})
    MI(A::Array{Int64})
    MI(A::Array{MP})
    MI(A::Array{MI})

Promote the given dense array of elements either from classic algebra (Float64
or Int64) or array of Max-Plus or Min-Plus elements to a dense array of tropical
semi-ring (min, +) numbers.

# Example
```julia-repl
julia> A = MI([1.0 Inf; 0.0 4])
2×2 (min,+) dense matrix:
  1   .
  0   4
```
"""
MI(::Array)

################################################################################
###
### Min-Plus constants
###
################################################################################

# ==============================================================================
"""
    mi0
    zero(MI)

Return the absorbing (min,+) element also named zero: +∞.
For (min,+), zero is +∞ since min(+∞, x) = x.

# Examples
```julia-repl
julia> mi0
(min,+) .

julia> zero(MI)
(min,+) .

julia> mi0 + MI(5)
(min,+) 5
```
"""
mi0

# ==============================================================================
"""
    mi1
    mie
    one(MI)

Return the neutral (min,+) element also named one: 0.
For (min,+), one is 0 since 0 + x = x in classic algebra (the ⨂ operation).

# Examples
```julia-repl
julia> mi1
(min,+) 0

julia> one(MI)
(min,+) 0

julia> mi1 * MI(5)
(min,+) 5
```
"""
mi1

# ==============================================================================
"""
    mitop

Return the top element for (min,+): -∞.
This is the absorbing element for the max operation in (min,+) context.

# Examples
```julia-repl
julia> mitop
(min,+) -Inf
```
"""
mitop

################################################################################
###
### Min-Plus identity matrices
###
################################################################################

# ==============================================================================
"""
    miI

Return the (min,+) identity scaling operator for creating identity matrices.
Diagonal elements are mi1 (= 0), off-diagonal elements are mi0 (= +∞).

# Examples
```julia-repl
julia> Matrix{MI}(miI, 2, 2)
2×2 (min,+) dense matrix:
  0   .
  .   0

julia> eye(MI, 3)
3×3 (min,+) dense matrix:
  0   .   .
  .   0   .
  .   .   0
```
"""
miI

################################################################################
###
### Min-Plus operations
###
################################################################################

# ==============================================================================
"""
    +(x::MI, y::MI)

(min,+) addition: returns the minimum of x and y.

# Examples
```julia-repl
julia> MI(3) + MI(5)
(min,+) 3

julia> MI(3) + mi0
(min,+) 3
```
"""
Base.:(+)(::MI, ::MI)

# ==============================================================================
"""
    *(x::MI, y::MI)

(min,+) multiplication: returns the sum of x and y (in classic algebra).

# Examples
```julia-repl
julia> MI(3) * MI(5)
(min,+) 8

julia> MI(3) * mi1
(min,+) 3
```
"""
Base.:(*)(::MI, ::MI)

# ==============================================================================
"""
    -(x::MI, y::MI)

Minus operator does not exist in (min,+) algebra. Throws an error.

# Examples
```julia-repl
julia> MI(3) - MI(5)
ERROR: Minus operator does not exist in (min,+) algebra
```
"""
Base.:(-)(::MI, ::MI)

# ==============================================================================
"""
    inv(x::MI)

(min,+) inverse: returns -x in classic algebra.

# Examples
```julia-repl
julia> inv(MI(3))
(min,+) -3
```
"""
Base.inv(::MI)

# ==============================================================================
"""
    /(x::MI, y::MI)

(min,+) division (residuation).

# Examples
```julia-repl
julia> MI(8) / MI(3)
(min,+) 5
```
"""
Base.:(/)(::MI, ::MI)
