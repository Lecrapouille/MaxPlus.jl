# ##############################################################################
#
# State space representation of Max-Plus linear systems.
# Creation of Max-Plus dynamical linear systems in implicit state form:
#    X(n) = D.X(n) + A.X(n-1) + B.U(n),
#    Y(n) = C.X(n)
#
# ##############################################################################

# Utility function: convert a vector to a dense matrix
v2m(v::Vector) = reshape(v, length(v), 1)
v2m(v::Array) = v

# Utility function: return an unified size for a matrix or a vector
size2(A::SpaMP) = (size(A, 1), size(A, 2))
size2(A::Array) = (size(A, 1), size(A, 2))

# Max-Plus dense matrix of zeros
mpfullzeros(m::Int64, n::Int64) = full(spzeros(MP, m, n))

# ==============================================================================
# Implicit dynamic linear Max-Plus system.

"""
    MPSysLin(A::ArrMP, B::ArrMP, C::ArrMP, D::ArrMP, x0::ArrMP)

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
                    full(mpzeros(3, 1)))           # Max-Plus column vector x0

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 3×3 Matrix{MP{Int64}}:
  0   .   .
  .   0   .
  .   .   0

A = 3×3 Matrix{MP{Int64}}:
  1   2   3
  4   5   6
  7   8   9

B = 3-element Vector{MP{Int64}}:
  0
  0
  0

C = 1×3 Matrix{MP{Int64}}:
  0   0   0

x0 = 3-element Vector{MP{Int64}}:
  .
  .
  .

julia> typeof(S)
MPSysLin{Int64}

julia> S.A
3×3 Matrix{MP{Float64}}:
  1   2   3
  4   5   6
  7   8   9
```
"""
struct MPSysLin
    A::ArrMP
    B::ArrMP
    C::ArrMP
    D::ArrMP
    x0::ArrMP

    ### Dense Max-Plus matrices

    # Constructor with dense Max-Plus matrices: A, B, C, D and x0
    function MPSysLin(A::ArrMP, B::ArrMP, C::ArrMP, D::ArrMP, x0::ArrMP)
        (ma,na) = size2(A)
        (ma != na) && error("Matrix A shall be squared")
        (mb,nb) = size2(B)
        ((mb != na) && (mb != 0)) && error("The row dimension of B: $(mb) is not in accordance with dimensions of A: $(na)")
        (mc,nc) = size2(C)
        ((nc != na) && (mc != 0)) && error("The column dimension of C: $(nc) is not in accordance with dimensions of A: $(na)")
        (mx0,nx0) = size2(x0)
        ((mx0 != na) || (nx0 != min(na, 1))) && error("Dimensions of x0: $(mx0) do not in accordance with dimensions of A: $(na)")
        (md,nd) = size2(D)
        ((md != na) || (nd != na)) && error("The column dimension of D: $(md)×$(nd) is not in accordance with dimensions of A: $(na)")
        new(A, v2m(B), v2m(C), D, v2m(x0))
    end

    # Constructor with dense Max-Plus matrices: A, B, C, D and implicit x0 (set as zeros)
    function MPSysLin(A::ArrMP, B::ArrMP, C::ArrMP, D::ArrMP)
        new(A, B, C, D, mpfullzeros(T, size(A,2), 1))
    end

    # Constructor with dense Max-Plus matrices: A, B, C and implicit D and x0 (set as zeros)
    function MPSysLin(A::ArrMP, B::ArrMP, C::ArrMP)
        na = size(A,2)
        new(A, B, C, mpeye(T, na, na), mpfullzeros(T, na, 1))
    end

    ### Sparse Max-Plus matrices

    # Constructor with sparse Max-Plus matrices: A, B, C, D and x0
    function MPSysLin(A::SpaMP, B::SpaMP, C::SpaMP, D::SpaMP, x0::SpaMP)
       new(full(A), full(B), full(C), full(D), full(x0))
    end

    # Constructor with sparse Max-Plus matrices: A, B, C, D and implicit x0 (set as zeros)
    function MPSysLin(A::SpaMP, B::SpaMP, C::SpaMP, D::SpaMP)
        new(full(A), full(B), full(C), full(D), mpfullzeros(T, size(A,2), 1))
    end

    # Constructor with sparse Max-Plus matrices: A, B, C and implicit D and x0 (set as zeros)
    function MPSysLin(A::SpaMP, B::SpaMP, C::SpaMP)
        na = size(A,2)
        new(full(A), full(B), full(C), mpeye(T, na, na), mpfullzeros(T, na, 1))
    end
end # MPSysLin

# ==============================================================================

function Base.:(==)(x::MPSysLin, y::MPSysLin)
    (x.A == y.A) && (x.B == y.B) && (x.C == y.C) && (x.D == y.D) && (x.x0 == y.x0)
end

# ==============================================================================

function Base.show(io::IO, S::MPSysLin)
    (@printf io "Implicit dynamic linear Max-Plus system:\n")
    (@printf io "  x(n) = D*x(n) + A*x(n-1) + B*u(n)\n  y(n) = C*x(n)\n  x(0) = x0\n\nwith:")
    (@printf io "\nD = "); mpshow(io, S.D)
    (@printf io "\nA = "); mpshow(io, S.A)
    (@printf io "\nB = "); mpshow(io, S.B)
    (@printf io "\nC = "); mpshow(io, S.C)
    (@printf io "\nx0 = "); mpshow(io, S.x0)
end

# ==============================================================================

function LaTeX(io::IO, S::MPSysLin)
    (@printf io "\\begin{array}{lcl}\n")
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
    (@printf io "\\end{array}")
end

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
                     full(mpzeros(3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(mpzeros(3, 1)));

julia> S1 + S2

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
function Base.:(+)(x::MPSysLin, y::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    MPSysLin([x.A mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.A],
             [x.B; y.B],
             [x.C y.C],
             [x.D mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.D],
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
                     full(mpzeros(3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(mpzeros(3, 1)));

julia> S1 * S2

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 6×6 Matrix{MP{Float64}}:
  0   .   .   .   .   .
  .   0   .   .   .   .
  .   .   0   .   .   .
  1   1   1   0   .   .
  2   2   2   .   0   .
  3   3   3   .   .   0

A = 6×6 Matrix{MP{Float64}}:
  10   11   12   .   .   .
  13   14   15   .   .   .
  16   17   18   .   .   .
   .    .    .   1   2   3
   .    .    .   4   5   6
   .    .    .   7   8   9

B = 6-element Vector{MP{Float64}}:
  0
  0
  0
  .
  .
  .

C = 1×6 Matrix{MP{Float64}}:
  .   .   .   4   5   6

x0 = 6-element Vector{MP{Float64}}:
  .
  .
  .
  .
  .
  .
```
"""
function Base.:(*)(y::MPSysLin, x::MPSysLin)
    n1 = size(x.A, 1)
    n2 = size(y.A, 1)
    nb1 = size(x.B, 2)
    nc2 = size(y.C, 1)
    MPSysLin([x.A mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.A],
             [x.B; mpfullzeros(T, n2, nb1)],
             [mpfullzeros(T, nc2, n1) y.C],
             [x.D mpfullzeros(T, n1, n2); y.B * x.C y.D],
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
                     full(mpzeros(3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(mpzeros(3, 1)));

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
    MPSysLin([x.A mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.A],
             [x.B mpfullzeros(T, n1, nb2); mpfullzeros(T, n2, nb1) y.B],
             [x.C mpfullzeros(T, nc1, n2); mpfullzeros(T, nc2, n1) y.C],
             [x.D mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.D],
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
                     full(mpzeros(3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(mpzeros(3, 1)));

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
    MPSysLin([x.A mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.A],
             [x.B; y.B],
             [x.C mpfullzeros(T, nc1, n2); mpfullzeros(T, nc2, n1) y.C],
             [x.D mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.D],
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
                     full(mpzeros(3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(mpzeros(3, 1)));

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
    MPSysLin([x.A mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.A],
             [x.B mpfullzeros(T, n1, nb2); mpfullzeros(T, n2, nb1) y.B],
             [x.C y.C],
             [x.D mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.D],
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
                     full(mpzeros(3, 1)));

julia> S2 = MPSysLin(MP([10.0 11 12;  13 14 15;  16 17 18]),
                     MP([    0.0;      0;      0]),
                     MP([    0.0       0       0]),
                     mpeye(3, 3),
                     full(mpzeros(3, 1)));

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
    MPSysLin([x.A mpfullzeros(T, n1, n2); mpfullzeros(T, n2, n1) y.A],
             [x.B; mpfullzeros(T, n2, nb1)],
             [x.C mpfullzeros(T, nc1, n2)],
             [x.D x.B * y.C; y.B * x.C y.D],
             [x.x0; y.x0])
end

# ==============================================================================
# From ScicosLab file: %mpls_m_talg.sci
# S1 * k

"""
    Base.:(*)(S::MPSysLin, k::U) where U

Diagonal composition of two Max-Plus linear systems.

# Examples
```julia-repl
julia> S = MPSysLin(MP([1.0 2 3;  4 5 6;  7 8 9]),
                    MP([    0.0;      0;      0]),
                    MP([    0.0       0       0]),
                    mpeye(3, 3),
                    full(mpzeros(3, 1)));

julia> 3 * S    # or S * 3 or MP(3.0) * S or ...

Implicit dynamic linear Max-Plus system:
  x(n) = D*x(n) + A*x(n-1) + B*u(n)
  y(n) = C*x(n)
  x(0) = x0

with:
D = 3×3 Matrix{MP{Float64}}:
  0   .   .
  .   0   .
  .   .   0

A = 3×3 Matrix{MP{Float64}}:
  1   2   3
  4   5   6
  7   8   9

B = 3-element Vector{MP{Float64}}:
  0
  0
  0

C = 1×3 Matrix{MP{Float64}}:
  3   3   3

x0 = 3-element Vector{MP{Float64}}:
  .
  .
  .
```
"""
function Base.:(*)(S::MPSysLin, k::U) where U
    MPSysLin(S.A, S.B * k, S.C, S.D, S.x0)
end

function Base.:(*)(k::U, S::MPSysLin) where U
    MPSysLin(S.A, S.B, S.C * k, S.D, S.x0)
end

# %mpls_m_talg.sci
function Base.:(*)(S::MPSysLin, M::SpaMP)
    MPSysLin(S.A, S.B * M, S.C, S.D, S.x0)
end

function Base.:(*)(M::SpaMP, S::MPSysLin)
    MPSysLin(S.A, M * S.B, S.C, S.D, S.x0)
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
""" # TODO
function mpexplicit(S::MPSysLin)
    ds = mpstar(S.D)
    bs = ds * S.B
    ac = [ds * S.A; S.C]
    zerocol = map(x -> Bool(x.λ), mpones(T, 1, size(ac, 1)) * (ac .!= mpzero(T)))
    keep = findall(zerocol[1,:])
    c = [S.C[i] for i in keep]
    MPSysLin([ac[i, j] for i in keep, j in keep],
             [bs[i] for i in keep],
             reshape(c, 1, length(c)))
end

# ==============================================================================

"""
    mpsimul(S::MPSysLin, u::ArrMP, history::Bool)

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
function mpsimul(S::MPSysLin, u::ArrMP, history::Bool)
    x = S.x0
    k = size(u, 1)
    if history
        Y = mpones(T, size(S.C, 1), k)
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
            Y[:,i] = S.C * x
        end
        Y
    else
        for i = 1:k
            x = S.A * x + S.B * u[i,:]
        end
        Y = S.C * x
        Y[:,1]
    end
end

function mpsimul(S::MPSysLin, u::VecMP, history::Bool)
    mpsimul(S, map(x -> MP(T(x.λ)), reshape(u, length(u), 1)), history)
end
