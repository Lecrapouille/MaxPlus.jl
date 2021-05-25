# ==============================================================================
# Check type of Max-Plus and type promotions on scalars
# ==============================================================================

a = MP(1.0)
@test typeof(a) == MP
@test a isa MP
@test a.λ isa Float64
@test a.λ == 1.0

b = MP(1)
@test typeof(b) == MP
@test b isa MP
@test b.λ isa Float64
@test b.λ == 1

c = MP(-1)
@test c isa MP
@test c.λ isa Float64
@test c.λ == -1.0

d = -MP(1)
@test d isa MP
@test d.λ isa Float64
@test d.λ == -1.0

# ==============================================================================
# Max-Plus constructor
# ==============================================================================

c = MP(1.0)
d = MP(c)
@test d isa MP

# ==============================================================================
# Max-Plus constructor pathological case
# ==============================================================================

c = MP(NaN)
@test c isa MP
@test c.λ isa Float64
@test c.λ == -Inf

# ==============================================================================
# Max-Plus comparaisons
# ==============================================================================

# Max-Plus scalars

b = MP(3.0)
@test (b == b) == true
@test (b != b) == (b ≠ b) == false
@test (b >= b) == (b ≥ b) == true
@test (b <= b) == (b ≤ b) == true
@test (b > b) == false
@test (b < b) == false

# Max-Plus vs. non Max-Plus comparaison

@test (b == 3.0) == true
@test (b == 4.0) == false
@test (b != 4.0) == true
@test (b < 4.0) == true
@test (b <= 4.0) == true
@test (b > 4.0) == false
@test (b >= 4.0) == false
@test (4.0 > b) == true
@test (4.0 >= b) == true

# Dense Max-Plus matrix comparaison

B = [MP(1) MP(2); MP(3) MP(4)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# Sparse Max-Plus matrix comparaison

S = sparse([MP(1) MP(2); MP(3) MP(4)])
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# Sparse/Dense Max-Plus matrix comparaison

@test (S == B) == (B == S) == true

# Sort is using isless()
v = MP([3, 1, 2]);
@test sort(v) == MP([1, 2, 3])

# ==============================================================================
# Max-Plus constants: absorbing and neutral
# ==============================================================================

@test mp0 isa MP
@test mp0.λ isa Float64
@test mp0.λ == -Inf
@test iszero(mp0) == true
@test isone(mp0) == false

@test ε isa MP
@test ε.λ isa Float64
@test ε.λ == -Inf
@test iszero(ε) == true
@test isone(ε) == false

@test mp1 isa MP
@test mp1.λ isa Float64
@test mp1.λ == 0.0
@test iszero(mp1) == false
@test isone(mp1) == true

@test mpe isa MP
@test mpe.λ isa Float64
@test mpe.λ == 0.0
@test iszero(mpe) == false
@test isone(mpe) == true

@test mptop isa MP
@test mptop.λ isa Float64
@test mptop.λ == Inf

# ==============================================================================
# Max-Plus type promotion/contamination on dense/sparse matrix
# ==============================================================================

# Dense matrix

A = [1 2; 3 4]
@test A isa Matrix{Int64}
@test MP(A) isa Matrix{MP}

B = [MP(1) MP(2); MP(3) MP(4)]
@test B isa Matrix{MP}

C = [MP(1) 2; 3 4]
@test C isa Matrix{MP}

D = [MP(1.0) 2; 3 4]
@test D isa Matrix{MP}

E = [MP(1) 2.0; 3 4]
@test E isa Matrix{MP}

F = [2 MP(1.0); 3 4]
@test F isa Matrix{MP}

# Sparse matrix

sA = sparse([1 2; 3 4])
@test sA isa SparseMatrixCSC{Int64, Int64}
@test MP(sA) isa SparseMatrixCSC{MP, Int64}

sB = sparse([MP(1) MP(2); MP(3) MP(4)])
@test sB isa SparseMatrixCSC{MP, Int64}

sC = sparse([MP(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MP, Int64}

sD = sparse([MP(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MP, Int64}

sE = sparse([MP(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MP, Int64}

sF = sparse([2 MP(1.0); 3 4])
@test sF isa SparseMatrixCSC{MP, Int64}

# Dense/Sparse matrix comparaison

@test (A == sA) == (sA == A)
@test (B == sB) == (sB == B)
@test (C == sC) == (sC == C)
@test (D == sD) == (sD == D)
@test (E == sE) == (sE == E)
@test (F == sF) == (sF == F)

# ==============================================================================
# Bug with Julia with SparseMatrixCSC and operator== which confused zero() and 0.0
# ==============================================================================

A = MP(sparse([1, 2], [1, 2], [0.0, 0.0]))
B = mpzeros(2,2)
@test A.nzval == MP([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mp0, mp0])
BB = mpzeros(2,2)
@test AA.nzval == MP([-Inf, -Inf])
@test (AA == BB) == true

# ==============================================================================
# Max-Plus sparse construction
# ==============================================================================

# Using SparseArray.sparse

spA = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test size(spA.nzval,1) == 3
@test spA.nzval == [mp0; 2.0; 0.0]

spB = sparse(MP([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test size(spB.nzval,1) == 3
@test spB.nzval == [mp0; 2.0; 0.0]
@test spB == spA

spC = MP(spA)
@test size(spC.nzval,1) == 3
@test spC.nzval == [mp0; 2.0; 0.0]
@test spC == spA

# Using Max-Plus

spA = sparse(MP([-Inf 0; 0 -Inf]))
@test findnz(spA) == ([2, 1], [1, 2], MP[0.0, 0.0])
spB = MP(sparse([-Inf 0; 0 -Inf]))
@test findnz(spB) == ([1, 2], [1, 2], MP[mp0, mp0])
@test spA != spB
spC = sparse(MP([4 0; 7 -Inf]))
@test findnz(spC) == ([1, 2, 1], [1, 1, 2], MP[4, 7, 0])

# ==============================================================================
# Max-Plus range construction
# ==============================================================================

V = MP(1:3)
@test V == MP([1; 2; 3])
V = MP(1.0:0.5:3.0)
@test V == MP([1.0; 1.5; 2.0; 2.5; 3.0])

# ==============================================================================
# Max-Plus scalar to class algebra scalar conversion
# ==============================================================================

# Scalar

@test plustimes(MP(2.0)) isa Float64
@test plustimes(MP(2.0)) == 2.0

# Dense Matrix

@test plustimes(C) isa Matrix{Float64}
@test plustimes(D) isa Matrix{Float64}
@test plustimes(E) isa Matrix{Float64}
@test plustimes(F) isa Matrix{Float64}

# Sparse Matrix

@test plustimes(sC) isa SparseMatrixCSC{Float64, Int64}
@test plustimes(sD) isa SparseMatrixCSC{Float64, Int64}
@test plustimes(sE) isa SparseMatrixCSC{Float64, Int64}
@test plustimes(sF) isa SparseMatrixCSC{Float64, Int64}

# Max-Plus Sparse Matrix to Max-Plus Dense Matrix

@test full(sC) isa Matrix{MP}
@test full(sD) isa Matrix{MP}
@test full(sE) isa Matrix{MP}
@test full(sF) isa Matrix{MP}

@test full(sC) == C
@test full(sD) == D
@test full(sE) == E
@test full(sF) == F

# dense() is a anlias for full()

@test dense(sC) isa Matrix{MP}
@test dense(sD) isa Matrix{MP}
@test dense(sE) isa Matrix{MP}
@test dense(sF) isa Matrix{MP}

@test dense(sC) == C
@test dense(sD) == D
@test dense(sE) == E
@test dense(sF) == F

# Max-Plus Sparse Matrix to classic algebra Dense Matrix

@test array(sC) isa Matrix{Float64}
@test array(sD) isa Matrix{Float64}
@test array(sE) isa Matrix{Float64}
@test array(sF) isa Matrix{Float64}

# Max-Plus sparse array to Max-Plus dense array

Z = dense(mpzeros(2,2))
@test typeof(Z) == Matrix{MP}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to Max-Plus dense array

Z = full(mpzeros(2,2))
@test typeof(Z) == Matrix{MP}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to non Max-Plus dense array

Z = array(mpzeros(2,2))
@test typeof(Z) == Matrix{Float64}
@test Z == [-Inf -Inf; -Inf -Inf]

# ==============================================================================
# Max-Plus zero
# ==============================================================================

@test zero(MP) == MP(-Inf) == -Inf
@test zero(MP) == mp0 == ε == -Inf
@test zero(MP(42.0)) == mp0 == ε == -Inf
@test zero(MP(42)) == mp0 == ε == -Inf

# ==============================================================================
# Max-Plus one
# ==============================================================================

@test one(MP) == MP(0.0) == 0.0
@test one(MP) == mp1 == mpe == 0.0
@test one(MP(42.0)) == mp1 == mpe == 0.0
@test one(MP) == mp1 == mpe == 0.0
@test one(MP(42)) == mp1 == mpe == 0

# ==============================================================================
# Max-Plus distributivity of operations
# ==============================================================================

a = MP(5.0)
b = MP(3.0)
c = MP(1.0)
@test (a + b) + c == a + (b + c) == a + b + c == MP(max(1.0, 3.0, 5.0)) == MP(5.0)
@test (a * b) * c == a * (b * c) == a * b * c == MP(1.0 + 3.0 + 5.0) == MP(9.0)
@test (a + b) * c == MP(max(5.0, 3.0) + 1.0) == MP(6.0)
@test (a * c) + (b * c) == MP(max(5.0 + 1.0, 3.0 + 1.0)) == MP(6.0)

# ==============================================================================
# mp0, mp1 and mptop operations
# ==============================================================================

@test -mp0 == mp0 == -ε == ε == -Inf

@test mp0 + mp0 == mp0 == ε + ε == ε == -Inf
@test mp0 + mp1 == ε + mpe == mpe + ε == mpe == 0
@test mp0 + mptop == mptop
@test mp1 + mp0 == mpe + ε == mpe == mp1
@test mp1 + mp1 == mpe + mpe == mpe == mp1
@test mp1 + mptop == mptop
@test mptop + mp0 == mptop
@test mptop + mp1 == mptop
@test mptop + mptop == mptop

@test mp0 * mp0 == ε * ε == mp0
@test mp1 * mp0 == mpe * ε == mp0
@test mptop * mp0 == mptop * ε == mp0
@test mp0 * mp1 == ε * mpe == mp0
@test mp1 * mp1 == mpe * mpe == mpe
@test mptop * mp1 == mptop
@test mp0 * mptop == mp0
@test mp1 * mptop == mptop
@test mptop * mptop == mptop

@test_broken mp0 / mp0 == mptop # FIXME
@test mp1 / mp0 == mptop
@test mptop / mp0 == mptop
@test mp0 / mp1 == mp0
@test mp1 / mp1 == mp1
@test mptop / mp1 == mptop
@test mp0 / mptop == mp0
@test mp1 / mptop == mp0
@test_broken mptop / mptop == mptop # FIXME

@test mp0 - mp1 == mp0
@test mp1 - mp1 == mp1
@test mptop - mp1 == mptop
@test mp0 - mp0 == mp0
@test_broken mp1 - mp0 == mp1 # FIXME
@test mptop - mp0 == mptop
@test mp0 - mptop == mp0
@test mp1 - mptop == mp0
@test mptop - mptop == mp0 # Scilab: NaN

# ==============================================================================
# Max-Plus operations and neutral elements and commutativity
# ==============================================================================

b = MP(3.0)
@test b + mp0 == mp0 + b == MP(max(3.0, -Inf)) == MP(max(-Inf, 3.0)) == b == MP(3.0)
@test b * mp1 == mp1 * b == MP(3.0 + 0.0) == MP(0.0 + 3.0) == b == MP(3.0)

@test b * b == MP(3.0 + 3.0) == MP(6.0)
@test b + b == MP(max(3.0, 3.0)) == MP(3.0)

a = MP(5.0)
b = MP(3.0)
@test a + b == MP(max(5.0, 3.0)) == b + a == MP(max(3.0, 5.0)) == MP(5.0)
@test a * b == MP(5.0 + 3.0) == b * a == MP(3.0 + 5.0) == MP(8.0)

# ==============================================================================
# Max-Plus between objects of different types
# ==============================================================================

b = MP(3.0)
@test b + 3.0 == 3.0 + b == b == MP(max(3.0, 3.0)) == MP(3.0)
@test b + 3   == 3   + b == b == MP(max(3.0, 3.0)) == MP(3.0)
@test b * 3.0 == 3.0 * b == MP(3.0 + 3.0) == MP(6.0)
@test b * 3   == 3   * b == MP(3.0 + 3.0) == MP(6.0)

# ==============================================================================
# Matrix ones, eye, zeros constructions
# ==============================================================================

# Matrix of ones

@test mpones(2) isa Vector{MP}
@test mpones(2) == [mp1; mp1]
@test mpones(2,5) isa Matrix{MP}
@test mpones(2,5) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]
@test mpones([1 2 3 4 5; 6 7 8 9 10]) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]
@test mpones(MP([1 2 3 4 5; 6 7 8 9 10])) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]

# Identity matrix

@test mpeye(2) isa Matrix{MP}
@test mpeye(2) == [mp1 mp0; mp0 mp1]
@test mpeye(2,5) isa Matrix{MP}
@test mpeye(2,5) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]
@test mpeye([1 2 3 4 5; 6 7 8 9 10]) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]
@test mpeye(MP([1 2 3 4 5; 6 7 8 9 10])) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]

# Matrix of zeros

@test mpzeros(2) isa SparseVector{MP, Int64}
@test mpzeros(2).nzval == MP([])
@test mpzeros(2,3) isa SparseMatrixCSC{MP, Int64}
@test mpzeros(2,3).nzval == MP([])
@test mpzeros([1 2 3; 6 7 8]).nzval == MP([])
@test mpzeros(MP([1 2 3; 6 7 8])).nzval == MP([])

# ==============================================================================
# Matrix ones, eye, zeros operations
# ==============================================================================

I = mpeye(4,4)
Z = mpzeros(4,4)
O = mpones(4,4)
A = MP(rand(4,4))

@test (A * I) == (I * A) == A
@test (A + Z) == (Z + A) == A
@test (A * Z) == (Z * A) == Z
@test (A + O) == (O + A) == A
@test (A * O) != (O * A) != A

# ==============================================================================
# Sparse Max-Plus element insertion
# ==============================================================================

# Note: ScicosLab will return isempty == false and length = 2
# Note: NSP, like Julia will return isempty == false and length = 4

O = mpzeros(2,2)
@test length(O) == 4
@test isempty(O) == false

# Zero elements are not inserted
O[1,1] = mp0
@test length(O) == 4
@test isempty(O) == false
@test O.nzval == []
@test isempty(nonzeros(O))

# Insert fake zero element
O[1,1] = MP(0.0)
@test O[1,1] == MP(0.0)
@test O.nzval == MP([0.0])
@test !isempty(nonzeros(O))

# Allow replacing a non-zero() element by a zero() element
O[1,1] = mp0
@test O[1,1] == MP(-Inf)
@test O.nzval == MP([-Inf])
@test !isempty(nonzeros(O))

# ==============================================================================
# Max-Plus element by element operations on matrices
# ==============================================================================

# Dense matrix

A = MP([1.0 2.0; 3.0 4.0])
@test 2.0 .+ A == A .+ 2.0 == [MP(2.0) + MP(1.0) MP(2.0) + MP(2.0); MP(2.0) + MP(3.0) MP(2.0) + MP(4.0)] == MP([2.0 2.0; 3.0 4.0])
@test 2.0 .* A == A .* 2.0 == [MP(2.0) * MP(1.0) MP(2.0) * MP(2.0); MP(2.0) * MP(3.0) MP(2.0) * MP(4.0)] == MP([3.0 4.0; 5.0 6.0])

# Sparse matrix

sA = sparse(A)
@test 2.0 .+ sA == sA .+ 2.0 == [MP(2.0) + MP(1.0) MP(2.0) + MP(2.0); MP(2.0) + MP(3.0) MP(2.0) + MP(4.0)] == MP([2.0 2.0; 3.0 4.0])
@test 2.0 .* sA == sA .* 2.0 == [MP(2.0) * MP(1.0) MP(2.0) * MP(2.0); MP(2.0) * MP(3.0) MP(2.0) * MP(4.0)] == MP([3.0 4.0; 5.0 6.0])

# ==============================================================================
# Max-Plus matrix dense/sparse multiplication
# ==============================================================================

A = [ε ε 3; -2 ε ε; ε ε 0]
B = [ε 2 ε; ε ε 0; -3 ε ε]
D = [3 ε ε; ε -2 ε; ε ε 0]
P = [ε ε 0; 0 ε ε; ε ε 0]

@test A == (D * P)
@test A == (sparse(D) * P)
@test A == (D * sparse(P))
@test A == (sparse(D) * sparse(P))
@test sparse(A) == (sparse(D) * sparse(P))

# ==============================================================================
# Max-Plus matrix addition
# ==============================================================================

A = MP([3.0 4; 5 6])
B = MP([1.0 2; 3 4])
C = MP([3.0 4; 5 6])
@test (A + B) == (B + A) == C

# ==============================================================================
# Max-Plus matrix division
# ==============================================================================

A = MP([3.0 4; 5 6])
B = MP([0.0 -2; 2 0])
@test_broken (A / A) == C # FIXME

# ==============================================================================
# Max-Plus trace
# ==============================================================================

A = [5 ε 5; ε 6 3; 11 12 11]

@test mptrace(A) == MP(11.0)
@test mptrace([]) == mp0
@test mptrace(mpeye(2,2)) == tr(mpeye(2,2)) == 0.0
@test mptrace(mpeye(2,5)) == 0.0
@test mptrace(mpzeros(2,2)) == mptrace(full(mpzeros(2,2))) == tr(mpzeros(2,2)) == mp0
@test mptrace(mpzeros(2,5)) == mptrace(full(mpzeros(2,5))) == mp0
@test mptrace([1.0 2.0; 3.0 4.0]) == tr(MP([1.0 2.0; 3.0 4.0])) == MP(1.0) + MP(4.0) == MP(4.0)
@test mptrace(sparse(MP([1.0 2.0; 3.0 4.0]))) == mptrace(sparse([1.0 2.0; 3.0 4.0])) == MP(4.0)
@test mptrace(mpzeros(2,2)) == mp0
@test sign(MP(0)) == MP(0)
@test sign(MP(5)) == MP(1)
@test sign(MP(-5)) == MP(-1)
@test sign(mp0) == MP(-1)

# ==============================================================================
# Max-Plus norm
# ==============================================================================

@test mpnorm(MP([1.0 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399.0)
@test mpnorm(MP([1 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399)
@test mpnorm(sparse(MP([1 20 2; 30 400 4; 4 50 10]))) == MP(400 - 1) == MP(399)
@test mpnorm([mp0 1; 10 mp0]) == MP(10.0) - mp0 == mptop

# ==============================================================================
# Max-Plus star
# From https://jpquadrat.github.io/TPALGLIN.pdf
# ==============================================================================

# Scalars

@test mpstar(MP(2)) == mptop
@test mpstar(MP(1.0)) == mptop
@test mpstar(MP(-1.0)) == mp1
@test mpstar(mp0) == mp1
@test mpstar(mp1) == mp1
@test mpstar(mptop) == mptop

# Matrices

@test_throws ErrorException("Matrix shall be squared") mpstar(mpeye(3,2))
@test mpstar(mpeye(2,2)) == mpeye(2,2)
@test mpstar(full(mpzeros(2,2))) == mpeye(2,2)

#

@test mpstar(MP([1 2; 3 4])) == [mptop mptop; mptop mptop]
A = MP([-3 -2; -1 0]); B = mpstar(A)
@test B == mpeye(2,2) + A
@test B == B + A^2
@test B == B + A^3

#

@test mpstar(MP([-1 2; mp0 -3])) == MP([0 2; -Inf 0])
A = [mp0 2 3; -2 -10 -1; -5 -2 mp1]
@test mpstar(A) == MP([0 2 3; -2 0 1; -4 -2 0])
@test mpstar(mpstar(A)) == mpstar(A)
@test mpstar(A) == (A^0 + A)^2
@test (A^0 + A)^2 == (A^0 + A)^3

# FIXME KO: Mixing +inf and -inf

# @test mpstar(MP([2 3; mp0 -1])) == [mptop mptop; mp0 mp1]

# Random large matrix

A = MP(rand(64,64))
@test mpstar(A) == fill(mptop, 64,64)

# FIXME KO
B = (((mpones(1, size(A,1)) * A * mpones(size(A,2), 1))[1,1])^-1) * A
@test_broken maximum(plustimes(B)) == 0.0

# ==============================================================================
# Max-Plus plus
# ==============================================================================

# Scalars

@test mpplus(MP(2)) == mptop
@test mpplus(MP(1.0)) == mptop
@test mpplus(MP(-1.0)) == MP(-1)
@test mpplus(mp0) == mp0
@test mpplus(mp1) == mp1
@test mpplus(mptop) == mptop

# Matrices

@test_throws ErrorException("Matrix shall be squared") mpplus(mpeye(3,2))
A = [mp0 2 3; -2 -10 -1; -5 -2 mp1]
B = mpplus(A)
@test B == MP([0 2 3; -2 0 1; -4 -2 0])
@test B == A * mpstar(A)

# FIXME donne le bon resultat dans REPL
A[2,1] = MP(-10)
@test_broken mpplus(A) == MP([-2 2 3; -6 -3 -1; -5 -2 0])

# What happens if a circuit has strictly positive weight ?
A[3,1] = 6
@test_broken mpplus(A) == fill(mptop, 2,3) # FIXME

# ==============================================================================
# Max-Plus a star b
# ==============================================================================

# TODO astarb

# ==============================================================================
# Max-Plus Howard
# ==============================================================================

S = sparse(MP([1 2; 3 4]))
λ,v = howard(S)
@test (λ,v) == (MP[4, 4], MP[2, 4])
@test (S * v) == (λ[1] * v)

S = sparse([mp0 2 mp0; mp1 mp0 mp0; mp0 mp0 2])
λ,v = howard(S)
@test (λ,v) == (MP[1, 1, 2], MP[1, 0, 2])
@test (S / Matrix(Diagonal(λ))) * v == v

S = sparse([1 1; mp0 2])
λ,v = howard(S)
@test (λ,v) == (MP[2, 2], MP[1, 2])
@test (S * v) == (λ[1] * v)
@test (S * [0; mp0]) == (MP(1) * [0; mp0])

S = sparse([2 1; mp0 mp1])
λ,v = howard(S)
@test (λ,v) == (MP[2, 0], MP[2, 0])
@test (S / Matrix(Diagonal(λ))) * v == v

# ==============================================================================
# Max-Plus power operator
# ==============================================================================

# Scalar

@test MP(2)^4   == MP(2 * 4) == MP(8)
@test MP(0)^0   == MP(0 * 0) == MP(0)
@test MP(2)^0   == MP(2 * 0) == MP(0)
@test MP(2)^-3  == MP(2)^(-3) == MP(2 * -3) == MP(-6)
@test MP(2)^0.5 == MP(2 * 0.5) == MP(1.0)
@test MP(2)^(-0.5) == MP(2 * -0.5) == MP(-1.0)
@test ε^0 == MP(0.0)
@test ε^2 == ε
@test_broken ε^(-2) == ε # FIXME
@test inv(MP(5)) == MP(5)^-1 == MP(-5)

# Matrix

A = MP([4 3; 7 -Inf])
@test A^0 == mpeye(2,2)
@test A * A == A^2 == MP([10.0 7.0; 11.0 10.0])
@test A * A * A == A^3 == MP([14.0 13.0; 17.0 14.0])
@test_broken A^-1

# ==============================================================================
# Max-Plus other operations
# ==============================================================================

@test abs2(MP(3.0)) == MP(6.0)
@test abs(MP(-3.0)) == MP(3.0)
@test abs(mp0) == mptop
@test float(MP(2)) == MP(2)
@test log(MP(2)) == MP(log(2.0)) == log(2.0)
@test round(MP(1.7)) == MP(round(1.7)) == round(1.7)
@test floor(MP(1.7)) == MP(floor(1.7)) == floor(1.7)

@test b / b == MP(3.0 - 3.0) == MP(0.0)
@test b - b == MP(3.0 - 3.0) == MP(0.0)
@test MP(3.0) / MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)
@test MP(3.0) - MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)

@test MP(3) \ MP(6) == MP(3)
@test MP(3) \ mp0 == mp0
@test MP(3) \ mp1 == MP(-3)
@test MP(3) \ mptop == mptop
@test mp0 \ mp1 == mptop

# ==============================================================================
# Max-Plus inverse
# ==============================================================================

A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
@test inv(A) == A^-1 == [mp0 -2 mp0; -1 mp0 mp0; mp0 mp0 -3]
@test A * inv(A) == inv(A) * A == mpeye(3,3)

A = [mp0 1 mp0; 2 mp0 mp0]
@test inv(A) == A^-1 == [mp0 -2; -1 mp0; mp0 mp0]
@test A * inv(A) == mpeye(2,2)
@test inv(A) * A == [0 mp0 mp0; mp0 0 mp0; mp0 mp0 mp0]
@test A * inv(A) != inv(A) * A

@test_throws ErrorException("The matrix cannot be inversed") inv(MP([1 2; 3 4]))

# ==============================================================================
# Max-Plus residu
# ==============================================================================

A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
B = [3 mp0 mp0; mp0 mp0 4; mp0 5 mp0]
x = A \ B
@test A * x == B
@test A \ A == mpeye(3,3)

# ==============================================================================
# Max-Plus min operator
# ==============================================================================

@test min(MP(3.0), mp0) == min(mp0, MP(3.0)) == mp0
@test min(MP(3.0), mp1) == min(mp1, MP(3.0)) == mp1
@test min(MP(1), MP(2)) == min(1, MP(2)) == min(MP(1), 2) == MP(1)
@test min(MP([10 1; 10 1]), MP([4 5; 6 5])) == MP([4 1; 6 1])
@test min(sparse(MP([10.0 1; 10 1])), sparse(MP([4.0 5; 6 5]))) == sparse(MP([4.0 1; 6 1]))
@test min(sparse(MP([10 1; 10 1])), sparse(MP([4 5; 6 5]))) == sparse(MP([4 1; 6 1]))
@test min(sparse(mpeye(2,2)), mpzeros(2,2)) == mpzeros(2,2)
@test min(mpeye(2,2), mpones(2,2)) == mpeye(2,2)

# ==============================================================================
# Max-Plus map on sparse matrix
# ==============================================================================

f = v -> v * 3
@test mpsparse_map(f, sparse(MP([5 2; 3 4]))) == sparse(MP([8 5; 6 7]))

# ==============================================================================
# Max-Plus display
# ==============================================================================

using Suppressor

mp_change_display(0)
result = @capture_out show(stdout, mp0)
@test result == "-Inf"
result = @capture_out show(stdout, mp1)
@test result == "0.0"
result = @capture_out show(stdout, MP(4.0))
@test result == "4.0"
result = @capture_out show(stdout, mpzero())
@test result == "-Inf"
result = @capture_out show(stdout, mpone())
@test result == "0.0"

mp_change_display(1)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, mpzero())
@test result == "."
result = @capture_out show(stdout, mpone())
@test result == "0"

mp_change_display(2)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, mpzero())
@test result == "."
result = @capture_out show(stdout, mpone())
@test result == "e"

mp_change_display(3)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MP(4))
@test result == "4"
result = @capture_out show(stdout, mpzero())
@test result == "ε"
result = @capture_out show(stdout, mpone())
@test result == "0"

mp_change_display(4)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, mpzero())
@test result == "ε"
result = @capture_out show(stdout, mpone())
@test result == "e"

A = MP([4.5 0.0; 7.0 -Inf])
result = @capture_out show(stdout, A)
@test result == "MP[4.5 e; 7 ε]"
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & e \\\\\n7 & \\varepsilon \\\\\n\\end{array}\n\\right]\n"
mp_change_display(0)
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & 0 \\\\\n7 & -\\infty \\\\\n\\end{array}\n\\right]\n"
