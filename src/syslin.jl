# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# A portage of the ScicosLab Max-Plus toolbox http://www.scicoslab.org/
# License: public domain
#
# Note: the documentation of functions for the REPL are placed in docstrings.jl
# ==============================================================================

export # Max-Plus Linear system
    MPSysLin, mpsimul, mpexplicit

# ##############################################################################
#
# State space representation of Max-Plus linear systems.
# Creation of Max-Plus dynamical linear systems in implicit state form:
#    X(n) = D.X(n) ⨁ A.X(n-1) ⨁ B.U(n),
#    Y(n) = C.X(n)
#
# ##############################################################################

# Utility function: return an unified size for a matrix or a vector
size2(A::SparseMatrixCSC) = (size(A, 1), size(A, 2))
size2(A::Array) = (size(A, 1), size(A, 2))
size2(A::Vector) = (size(A, 1), size(A, 2))

# ==============================================================================
# Implicit dynamic linear Max-Plus system.

"""
    MPSysLin(A::AbsArrMP, B::AbsArrMP, C::AbsArrMP [, D::AbsArrMP, [x0::AbsArrMP]])

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
struct MPSysLin
    A::AbsArrMP
    B::AbsArrMP
    C::AbsArrMP
    D::AbsArrMP
    x0::AbsArrMP

    # Constructor with Max-Plus matrices: A, B, C, D and x0
    function MPSysLin(A::AbsArrMP, B::AbsArrMP, C::AbsArrMP, D::AbsArrMP, x0::AbsArrMP)
        # show(stdout, (size2(A),size2(B),size2(C),size2(D),size2(x0)))

        (ma,na) = size2(A)
        (ma != na) &&
            error("Matrix A shall be squared")

        (mb,nb) = size2(B)
        (mb != na) && (mb != 0) && (nb != 0) &&
            error("The row dimension of B: $(mb) is not in accordance with dimensions of A: $(na)")

        (mc,nc) = size2(C)
        (nc != na) && (mc != 0) && (nc != 0) &&
            error("The column dimension of C: $(nc) is not in accordance with dimensions of A: $(na)")

        (mx0,nx0) = size2(x0)
        (mx0 != na) &&
            error("Dimensions of x0: $(mx0) do not in accordance with dimensions of A: $(na)")

        (md,nd) = size2(D)
        ((md != na) || (nd != na)) &&
            error("The column dimension of D: $(md)×$(nd) is not in accordance with dimensions of A: $(na)")

        new(A, B, C, D, x0)
    end

    # Constructor with Max-Plus matrices: A, B, C, D (and implicit x0 set as zeros)
    function MPSysLin(A::AbsArrMP, B::AbsArrMP, C::AbsArrMP, D::AbsArrMP)
        na = size(A,2)
        MPSysLin(A, B, C, D, spzeros(MP, na, 1))
    end

    # Constructor with Max-Plus matrices: A, B, C (and implicit D and x0 set as zeros)
    function MPSysLin(A::AbsArrMP, B::AbsArrMP, C::AbsArrMP)
        na = size(A,2)
        MPSysLin(A, B, C, spzeros(MP, na, na), spzeros(MP, na, 1))
    end
end # MPSysLin

# ==============================================================================

"""
    Base.:(==)(x::MPSysLin, y::MPSysLin)

Test whether two Max-Plus linear system are equal.
```
"""
function Base.:(==)(x::MPSysLin, y::MPSysLin)
    (x.A == y.A) && (x.B == y.B) && (x.C == y.C) && (x.D == y.D) && (x.x0 == y.x0)
end

# ==============================================================================

"""
    mpshow(io::IO, S::MPSysLin)

Base function for displaying a Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
function mpshow(io::IO, S::MPSysLin)
    (@printf io "Implicit dynamic linear Max-Plus system:\n")
    (@printf io "  x(n) = D*x(n) + A*x(n-1) + B*u(n)\n  y(n) = C*x(n)\n  x(0) = x0\n\nwith:")
    (@printf io "\nD = "); mpshow(io, S.D)
    (@printf io "\nA = "); mpshow(io, S.A)
    (@printf io "\nB = "); mpshow(io, S.B)
    (@printf io "\nC = "); mpshow(io, S.C)
    (@printf io "\nx0 = "); mpshow(io, S.x0)
end

"""
    show(io::IO, S::MPSysLin)

Display a Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
Base.show(io::IO, S::MPSysLin) = mpshow(io, S)

"""
    show(io::IO, ::MIME"text/plain", S::MPSysLin)

Display a Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
Base.show(io::IO, ::MIME"text/plain", S::MPSysLin) = mpshow(io, S)

# ==============================================================================

"""
    LaTeX(S::MPSysLin)

Return the LaTeX code as string from the given Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
function LaTeX(S::MPSysLin)
    "\\left\\{\\begin{array}{lcl}\nx_n & = & " *
    LaTeX(S.D) * " x_n \\oplus " *
    LaTeX(S.A) * " x_{n-1} \\oplus " *
    LaTeX(S.B) * " u_n\\\\ y_n & = & " *
    LaTeX(S.C) * " x_n\\\\ x_0 & = & " *
    LaTeX(S.x0) * "\\end{array}\\right."
end

"""
    LaTeX(io::IO, S::MPSysLin)

Display the LaTeX code from the given Max-Plus linear system.

# Examples
```julia-repl
julia>
```
"""
function LaTeX(io::IO, S::MPSysLin)
    (@printf io "\\left\\{\\begin{array}{lcl}\n")
    (@printf io "x_n & = & ")
    LaTeX(io, S.D)
    (@printf io " x_n \\oplus ")
    LaTeX(io, S.A)
    (@printf io " x_{n-1} \\oplus ")
    LaTeX(io, S.B)
    (@printf io " u_n\\\\")
    (@printf io "y_n & = & ")
    LaTeX(io, S.C)
    (@printf io " x_n\\\\")
    (@printf io "x_0 & = & ")
    LaTeX(io, S.x0)
    (@printf io "\\end{array}\\right.")
end

Base.show(io::IO, ::MIME"text/latex", S::MPSysLin) = LaTeX(io, S)

# ==============================================================================
# Parallel composition.
# From ScicosLab file: %mpls_a_mpls.sci

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
function Base.:(+)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; y.B],
             [x.C y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Series composition.
# From ScicosLab file: %mpls_m_mpls.sci

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
function Base.:(*)(y::MPSysLin, x::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc2 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; spzeros(MP, n2, nb1)],
             [spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Diagonal composition.
# From ScicosLab file: %mpls_g_mpls.sci

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
function Base.:(|)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nb2 = size(y.B, 2)
    nc1 = size(x.C, 1)
    nc2 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B spzeros(MP, n1, nb2); spzeros(MP, n2, nb1) y.B],
             [x.C spzeros(MP, nc1, n2); spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# computes [S1;S2] that is : inputs in common, concatenation of outputs.
# From ScicosLab file: %mpls_f_mpls.sci

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
function Base.:vcat(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nc1 = size(x.C, 1)
    nc2 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; y.B],
             [x.C spzeros(MP, nc1, n2); spzeros(MP, nc2, n1) y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# computes [S1,S2] that is concatenation of inputs, addition of outputs
# From ScicosLab file: %mpls_c_mpls.sci

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
function Base.:hcat(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nb2 = size(y.B, 2)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B spzeros(MP, n1, nb2); spzeros(MP, n2, nb1) y.B],
             [x.C y.C],
             [x.D spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# Feedback composition: computes star(S1*S2)*S1 in state-space form.
# From ScicosLab file: %mpls_v_mpls.sci

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
function Base.:(/)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc1 = size(y.C, 1)
    MPSysLin([x.A spzeros(MP, n1, n2); spzeros(MP, n2, n1) y.A],
             [x.B; spzeros(MP, n2, nb1)],
             [x.C spzeros(MP, nc1, n2)],
             [x.D x.B * y.C; y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# From ScicosLab file: mplssize.sci
function Base.:size(S::MPSysLin)
    [size(S.B, 1), size(S.B, 2), size(S.C, 1)]
end

# ==============================================================================
# From ScicosLab file: explicit.sci

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
function mpexplicit(S::MPSysLin)
    ds = mpstar(S.D)
    A = ds * S.A
    B = ds * S.B
    AC = [A; S.C]
    MPSysLin(AC, B, S.C,
             spzeros(MP, size(AC,1), size(AC,2)), S.x0)
end

# ==============================================================================

"""
    mpsimul(S::MPSysLin, u::AbsArrMP, history::Bool)

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
function mpsimul(S::MPSysLin, u::AbsArrMP, history::Bool)
    x = S.x0
    k = size(u, 1)
    if history
        Y = mpones(k, size(S.C, 1))
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
            Y[i,:] = S.C * x
        end
        Y
    else
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
        end
        Y = S.C * x
        Y[1,:]
    end
end

function mpsimul(S::MPSysLin, u::VecMP, history::Bool)
    mpsimul(S, map(x -> MP(x.λ), reshape(u, length(u), 1)), history)
end
