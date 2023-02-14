# ==============================================================================
# Docstring for Max-Plus linear systems
# ==============================================================================

################################################################################
###
### Max-Plus Linear system
###
################################################################################

# ==============================================================================
"""
    MPSysLin(A::MPAbstractVecOrMat, B::MPAbstractVecOrMat, C::MPAbstractVecOrMat [, D::MPAbstractVecOrMat, [x0::MPAbstractVecOrMat]])

Structure for state space representation of Max-Plus linear systems.
Creation of max-plus dynamical linear systems in implicit state form:
```
    X(n) = D ⨂ X(n) ⨁ A ⨂ X(n-1) ⨁ B ⨂ U(n),
    Y(n) = C ⨂ X(n)
```

`MPSysLin` is the equivalent to ScicosLab function `mpsyslin`. It accepts dense
or sparse Max-Plus arrays. `D` and `x0` can be omited: they will be filled
with Max-Plus zeros `ε`.

# Arguments
- `A` shall be squared.
- `B` shall has its row dimension in accordance with dimensions of `A`.
- `C` shall has its column dimension in accordance with dimensions of `A`.
- `D` shall has its dimension in accordance with dimensions of `A`.
- `x0` shall has its dimensions in accordance with dimensions of `A`.

Sizes are checked and an error is thrown during the compilation.

Elements of the system can access directly ie. `S.A`, `S.B`, `S.C`, `S.D`,
`S.x0`.

# Examples
```julia-repl
julia> S = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]), # Max-Plus matrix A
                    MP([    0.0;      0;      0]), # Max-Plus column vector B
                    MP([    0.0       0       0]), # Max-Plus row vector C
                    mpeye(3, 3),                   # Max-Plus matrix D
                    full(zeros(MP, 3, 1)))         # Max-Plus column vector x0

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 3×3 Max-Plus dense matrix:
  0   .   .
  .   0   .
  .   .   0

A = 3×3 Max-Plus dense matrix:
  1   2   3
  4   5   6
  7   8   9

B = 3-element Max-Plus vector:
  0
  0
  0

C = 1×3 Max-Plus dense matrix:
  0   0   0

x0 = 3-element Max-Plus vector:
  .
  .
  .

julia> typeof(S)
MPSysLin

julia> S.A
3×3 Max-Plus dense matrix:
  1   2   3
  4   5   6
  7   8   9
```
"""
MPSysLin

# ==============================================================================
"""
    mpexplicit(S::MPSysLin)

Convert to an explicit linear system.

# Examples
```julia-repl
julia> S = MPSysLin(MP([1.0 2; 3 4]),
                    MP([0.0; 0]),
                    MP([0.0 0]),
                    mpeye(Float64, 2,2))

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 2×2 Matrix{MP{Float64}}:
  0   .
  .   0

A = 2×2 Matrix{MP{Float64}}:
  1   2
  3   4

B = 2-element Vector{MP{Float64}}:
  0
  0

C = 1×2 Matrix{MP{Float64}}:
  0   0

x0 = 2-element Vector{MP{Float64}}:
  .
  .

julia> mpexplicit(S)

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 2×2 Matrix{MP{Float64}}:
  0   .
  .   0

A = 2×2 Matrix{MP{Float64}}:
  1   2
  3   4

B = 2-element Vector{MP{Float64}}:
  0
  0

C = 1×2 Matrix{MP{Float64}}:
  0   0

x0 = 2-element Vector{MP{Float64}}:
  .
  .
```
"""
mpexplicit(S::MPSysLin)

# ==============================================================================
"""
    mpsimul(S::MPSysLin, u::MPAbstractVecOrMat, history::Bool)

Compute states X of an autonomous linear max-plus system:
```
    x(n+1) = A ⨂ x(n)`  for n = 0 .. k
```

Where:
- `A` is a system matrix
- `x0` is a initial state vector,
- `k` is the number of iterations
- when history is set to true save all computed states, else return the last one.

# Arguments
- u: inputs to inject in the system
- history: if true then return a vector holding all results, else return the
  last result.

# Examples
```julia-repl
julia> S = MPSysLin(MP([1.0 2; 3 4]),
                    MP([0.0; 0]),
                    MP([0.0 0]),
                    mpeye(Float64, 2,2))

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 2×2 Matrix{MP{Float64}}:
  0   .
  .   0

A = 2×2 Matrix{MP{Float64}}:
  1   2
  3   4

B = 2-element Vector{MP{Float64}}:
  0
  0

C = 1×2 Matrix{MP{Float64}}:
  0   0

x0 = 2-element Vector{MP{Float64}}:
  .
  .

julia> mpsimul(S, MP(1:10))
1×10 Array{MP{Float64},2}:
 1  5  9  13  17  21  25  29  33  37

julia> mpsimul(S1, MP(1:10), history=false)
1×1 Array{MP{Float64},2}:
 37
```
"""
mpsimul(S::MPSysLin, u::MPAbstractVecOrMat, history::Bool)

# ==============================================================================
"""
    Base.:(+)(y::MPSysLin, x::MPSysLin)

Parallel composition of two Max-Plus linear systems.

# Examples
```julia-repl
julia> S1 = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                     MP([    1.0;      2;      3]),
                     MP([    4.0       5       6]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S1 + S2
Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 6×6 Max-Plus dense matrix:
  0   .   .   .   .   .
  .   0   .   .   .   .
  .   .   0   .   .   .
  .   .   .   0   .   .
  .   .   .   .   0   .
  .   .   .   .   .   0

A = 6×6 Max-Plus dense matrix:
  1   2   3    .    .    .
  4   5   6    .    .    .
  7   8   9    .    .    .
  .   .   .   10   11   12
  .   .   .   13   14   15
  .   .   .   16   17   18

B = 6-element Max-Plus vector:
  1
  2
  3
  0
  0
  0

C = 1×6 Max-Plus dense matrix:
  4   5   6   0   0   0

x0 = 6-element Max-Plus vector:
  .
  .
  .
  .
  .
  .
```
"""
Base.:(+)(x::MPSysLin, y::MPSysLin)

# ==============================================================================
"""
    Base.:(*)(y::MPSysLin, x::MPSysLin)

Series composition of two Max-Plus linear systems.

# Examples
```julia-repl
julia> S1 = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                     MP([    1.0;      2;      3]),
                     MP([    4.0       5       6]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S1 * S2
TODO
```
"""
Base.:(*)(y::MPSysLin, x::MPSysLin)

# ==============================================================================
"""
    Base.:(|)(y::MPSysLin, x::MPSysLin)

Diagonal composition of two Max-Plus linear systems.

# Examples
```julia-repl
julia> S1 = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                     MP([    1.0;      2;      3]),
                     MP([    4.0       5       6]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S1 | S2

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 6×6 Matrix{MP{Float64}}:
  0   .   .   .   .   .
  .   0   .   .   .   .
  .   .   0   .   .   .
  .   .   .   0   .   .
  .   .   .   .   0   .
  .   .   .   .   .   0

A = 6×6 Matrix{MP{Float64}}:
  1   2   3    .    .    .
  4   5   6    .    .    .
  7   8   9    .    .    .
  .   .   .   10   11   12
  .   .   .   13   14   15
  .   .   .   16   17   18

B = 6×2 Matrix{MP{Float64}}:
  1   .
  2   .
  3   .
  .   0
  .   0
  .   0

C = 2×6 Matrix{MP{Float64}}:
  4   5   6   .   .   .
  .   .   .   0   0   0

x0 = 6-element Vector{MP{Float64}}:
  .
  .
  .
  .
  .
  .
```
"""
Base.:(|)(x::MPSysLin, y::MPSysLin)

# ==============================================================================
"""
    Base.:vcat(x::MPSysLin, y::MPSysLin)

Composition of two Max-Plus linear systems: inputs in common, concatenation of
outputs.

# Examples
```julia-repl
julia> S1 = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                     MP([    1.0;      2;      3]),
                     MP([    4.0       5       6]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> [S1 ; S2]
Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 6×6 Matrix{MP{Float64}}:
  0   .   .   .   .   .
  .   0   .   .   .   .
  .   .   0   .   .   .
  .   .   .   0   .   .
  .   .   .   .   0   .
  .   .   .   .   .   0

A = 6×6 Matrix{MP{Float64}}:
  1   2   3    .    .    .
  4   5   6    .    .    .
  7   8   9    .    .    .
  .   .   .   10   11   12
  .   .   .   13   14   15
  .   .   .   16   17   18

B = 6-element Vector{MP{Float64}}:
  1
  2
  3
  0
  0
  0

C = 2×6 Matrix{MP{Float64}}:
  4   5   6   .   .   .
  .   .   .   0   0   0

x0 = 6-element Vector{MP{Float64}}:
  .
  .
  .
  .
  .
  .
```
"""
Base.:vcat(x::MPSysLin, y::MPSysLin)

# ==============================================================================
"""
    Base.:hcat(x::MPSysLin, y::MPSysLin)

Composition of two Max-Plus linear systems: concatenation of inputs, addition of
outputs.

# Examples
```julia-repl
julia> S1 = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                     MP([    1.0;      2;      3]),
                     MP([    4.0       5       6]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> [S1 S2]

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 6×6 Matrix{MP{Float64}}:
  0   .   .   .   .   .
  .   0   .   .   .   .
  .   .   0   .   .   .
  .   .   .   0   .   .
  .   .   .   .   0   .
  .   .   .   .   .   0

A = 6×6 Matrix{MP{Float64}}:
  1   2   3    .    .    .
  4   5   6    .    .    .
  7   8   9    .    .    .
  .   .   .   10   11   12
  .   .   .   13   14   15
  .   .   .   16   17   18

B = 6×2 Matrix{MP{Float64}}:
  1   .
  2   .
  3   .
  .   0
  .   0
  .   0

C = 1×6 Matrix{MP{Float64}}:
  4   5   6   0   0   0

x0 = 6-element Vector{MP{Float64}}:
  .
  .
  .
  .
  .
  .
```
"""
Base.:hcat(x::MPSysLin, y::MPSysLin)

# ==============================================================================
"""
    Base.:(/)(y::MPSysLin, x::MPSysLin)

Feedback composition: computes `mpstar(S1*S2)*S1` in state-space form.

# Examples
```julia-repl
julia> S1 = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                     MP([    1.0;      2;      3]),
                     MP([    4.0       5       6]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(zeros(MP, 3, 1)));

julia> S1 / S2

Implicit dynamic linear Max-Plus system:
Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 6×6 Matrix{MP{Float64}}:
  0   .   .   1   1   1
  .   0   .   2   2   2
  .   .   0   3   3   3
  4   5   6   0   .   .
  4   5   6   .   0   .
  4   5   6   .   .   0

A = 6×6 Matrix{MP{Float64}}:
  1   2   3    .    .    .
  4   5   6    .    .    .
  7   8   9    .    .    .
  .   .   .   10   11   12
  .   .   .   13   14   15
  .   .   .   16   17   18

B = 6-element Vector{MP{Float64}}:
  1
  2
  3
  .
  .
  .

C = 1×6 Matrix{MP{Float64}}:
  4   5   6   .   .   .

x0 = 6-element Vector{MP{Float64}}:
  .
  .
  .
  .
  .
  .
```
"""
Base.:(/)(x::MPSysLin, y::MPSysLin)

# ==============================================================================
"""
    Base.:(==)(x::MPSysLin, y::MPSysLin)

Test whether two Max-Plus linear system are equal.
```
"""
Base.:(==)(x::MPSysLin, y::MPSysLin)

# ==============================================================================
"""
    mpshow(io::IO, S::MPSysLin)

Base function for displaying a Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
mpshow(io::IO, S::MPSysLin)

# ==============================================================================
"""
    show(io::IO, S::MPSysLin)

Display a Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
Base.show(io::IO, S::MPSysLin)

# ==============================================================================
"""
    show(io::IO, ::MIME"text/plain", S::MPSysLin)

Display a Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
Base.show(io::IO, ::MIME"text/plain", S::MPSysLin)

# ==============================================================================
"""
    LaTeX(S::MPSysLin)

Return the LaTeX code as string from the given Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
LaTeX(S::MPSysLin)

# ==============================================================================
"""
    LaTeX(io::IO, S::MPSysLin)

Display the LaTeX code from the given Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
LaTeX(io::IO, S::MPSysLin)
