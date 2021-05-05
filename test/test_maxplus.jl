# ==============================================================================
# Check type of Max-Plus and type promotions on scalars
# ==============================================================================

a = MP(1.0)
@test typeof(a) == MP{Float64}
@test a isa MP{Float64}
@test a.λ isa Float64
@test a.λ == 1.0

b = MP(1)
@test typeof(b) == MP{Int64}
@test b isa MP{Int64}
@test b.λ isa Int64
@test b.λ == 1

c = MP{Float64}(1)
@test c isa MP{Float64}
@test c.λ isa Float64
@test c.λ == 1.0

# ==============================================================================
# Max-Plus constructor
# ==============================================================================

c = MP(1.0)
d = MP(c)
@test d isa MP{Float64}

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

# ==============================================================================
# Max-Plus constants: absorbing and neutral
# ==============================================================================

@test mp0 isa MP{Float64}
@test mp0.λ isa Float64
@test mp0.λ == -Inf
@test iszero(mp0) == true
@test isone(mp0) == false

@test ϵ isa MP{Float64}
@test ϵ.λ isa Float64
@test ϵ.λ == -Inf
@test iszero(ϵ) == true
@test isone(ϵ) == false

@test mp1 isa MP{Float64}
@test mp1.λ isa Float64
@test mp1.λ == 0.0
@test iszero(mp1) == false
@test isone(mp1) == true

@test mpe isa MP{Float64}
@test mpe.λ isa Float64
@test mpe.λ == 0.0
@test iszero(mpe) == false
@test isone(mpe) == true

@test mptop isa MP{Float64}
@test mptop.λ isa Float64
@test mptop.λ == Inf

# ==============================================================================
# Max-Plus type promotion/contamination on dense/sparse matrix
# ==============================================================================

# Dense matrix

A = [1 2; 3 4]
@test A isa Matrix{Int64}
@test MP(A) isa Matrix{MP{Int64}}

B = [MP(1) MP(2); MP(3) MP(4)]
@test B isa Matrix{MP{Int64}}

C = [MP(1) 2; 3 4]
@test C isa Matrix{MP{Int64}}

D = [MP(1.0) 2; 3 4]
@test D isa Matrix{MP{Float64}}

E = [MP(1) 2.0; 3 4]
@test E isa Matrix{MP{Int64}}

F = [2 MP(1.0); 3 4]
@test F isa Matrix{MP{Float64}}

# Sparse matrix

sA = sparse([1 2; 3 4])
@test sA isa SparseMatrixCSC{Int64, Int64}
@test MP(sA) isa SparseMatrixCSC{MP{Int64}, Int64}

sB = sparse([MP(1) MP(2); MP(3) MP(4)])
@test sB isa SparseMatrixCSC{MP{Int64}, Int64}

sC = sparse([MP(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MP{Int64}, Int64}

sD = sparse([MP(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MP{Float64}, Int64}

sE = sparse([MP(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MP{Int64}, Int64}

sF = sparse([2 MP(1.0); 3 4])
@test sF isa SparseMatrixCSC{MP{Float64}, Int64}

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
B = mpzeros(Float64, 2,2)
@test A.nzval == MP([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mp0, mp0])
BB = mpzeros(Float64, 2,2)
@test AA.nzval == MP([-Inf, -Inf])
@test (AA == BB) == true

# ==============================================================================
# Max-Plus sparse construction
# ==============================================================================

# Using SparseArray.sparse

spA = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test size(spA.nzval,1) == 3
@test spA.nzval == [mp0; 2.0; 0.0]

spB = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]), keepzeros=false)
@test size(spB.nzval,1) == 2
@test spB.nzval == [2.0; 0.0]
@test spB == spA

spC = MP(spA)
@test size(spC.nzval,1) == 3
@test spC.nzval == [mp0; 2.0; 0.0]
@test spC == spA

# Using Max-Plus

spA = mpsparse([-Inf 0; 0 -Inf])
@test findnz(spA) == ([2, 1], [1, 2], MP{Float64}[0.0, 0.0])
spB = mpsparse([-Inf 0; 0 -Inf], keepzeros=false)
@test findnz(spB) == (Int64[], Int64[], MP{Float64}[])
@test spA != spB

spC = mpsparse(MP([4 0; 7 -Inf]))
@test findnz(spC) == ([1, 2, 1], [1, 1, 2], MP{Float64}[4, 7, 0])

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

@test plustimes(C) isa Matrix{Int64}
@test plustimes(D) isa Matrix{Float64}
@test plustimes(E) isa Matrix{Int64}
@test plustimes(F) isa Matrix{Float64}

# Sparse Matrix

@test plustimes(sC) isa SparseMatrixCSC{Int64, Int64}
@test plustimes(sD) isa SparseMatrixCSC{Float64, Int64}
@test plustimes(sE) isa SparseMatrixCSC{Int64, Int64}
@test plustimes(sF) isa SparseMatrixCSC{Float64, Int64}

# Max-Plus Sparse Matrix to Max-Plus Dense Matrix

@test full(sC) isa Matrix{MP{Int64}}
@test full(sD) isa Matrix{MP{Float64}}
@test full(sE) isa Matrix{MP{Int64}}
@test full(sF) isa Matrix{MP{Float64}}

@test full(sC) == C
@test full(sD) == D
@test full(sE) == E
@test full(sF) == F

# dense() is a anlias for full()

@test dense(sC) isa Matrix{MP{Int64}}
@test dense(sD) isa Matrix{MP{Float64}}
@test dense(sE) isa Matrix{MP{Int64}}
@test dense(sF) isa Matrix{MP{Float64}}

@test dense(sC) == C
@test dense(sD) == D
@test dense(sE) == E
@test dense(sF) == F

# Max-Plus Sparse Matrix to classic algebra Dense Matrix

@test array(sC) isa SparseMatrixCSC{Int64, Int64}
@test array(sD) isa SparseMatrixCSC{Float64, Int64}
@test array(sE) isa SparseMatrixCSC{Int64, Int64}
@test array(sF) isa SparseMatrixCSC{Float64, Int64}

# Convert to Min-Plus algebra

@test minplus(MP(3.0)) == MP(3.0)
@test minplus(mp1) == mp1
@test minplus(mp0) == mptop
@test minplus(mptop) == mp0
@test minplus(MP([Inf 0.0; 7 -Inf])) == MP([-Inf 0.0; 7 Inf])
A = MP([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
@test minplus(A) == MP([0 3 -Inf 1; 1 2 2 Inf; Inf -Inf 1 0])
@test minplus(minplus(A)) == A

# Max-Plus sparse array to Max-Plus dense array

Z = dense(mpzeros(Float64, 2,2))
@test typeof(Z) == Matrix{MP{Float64}}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to Max-Plus dense array

Z = full(mpzeros(Float64, 2,2))
@test typeof(Z) == Matrix{MP{Float64}}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to non Max-Plus dense array

Z = array(mpzeros(Float64, 2,2))
@test typeof(Z) == SparseMatrixCSC{Float64,Int64}
@test Z == [-Inf -Inf; -Inf -Inf]

# ==============================================================================
# Max-Plus zero
# ==============================================================================

@test zero(MP{Float64}) == MP(-Inf) == -Inf
@test zero(MP{Float64}) == mp0 == ϵ == -Inf
@test zero(MP(42.0)) == mp0 == ϵ == -Inf
@test zero(MP{Int64}) == -9223372036854775808
@test zero(MP(42)) == -9223372036854775808

# ==============================================================================
# Max-Plus one
# ==============================================================================

@test one(MP{Float64}) == MP(0.0) == 0.0
@test one(MP{Float64}) == mp1 == mpe == 0.0
@test one(MP(42.0)) == mp1 == mpe == 0.0
@test one(MP{Int64}) == mp1 == mpe == 0.0
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

@test mp0 + mp0 == mp0 == ϵ + ϵ == ϵ == -Inf
@test mp0 + mp1 == ϵ + mpe == mpe + ϵ == mpe == 0
@test mp0 + mptop == mptop
@test mp1 + mp0 == mpe + ϵ == mpe == mp1
@test mp1 + mp1 == mpe + mpe == mpe == mp1
@test mp1 + mptop == mptop

@test mp0 * mp0 == ϵ * ϵ == mp0
@test mp0 * mp1 == ϵ * mpe == mp0
@test mp1 * mp0 == mpe * ϵ == mp0
@test mp1 * mp1 == mpe * mpe == mpe
@test mp1 * mptop == mptop

@test mp0 - mp1 == ϵ - mpe == mp0
@test mp1 - mp1 == mpe - mpe == mpe

@test mptop - mp0 == mptop
@test mptop - mp1 == mptop
@test mp0 - mptop == mp0
@test mp1 - mptop == mp0

# Border cases
# FIXME @test mp0 * mptop == mp0
# FIXME @test mp0 - mp0 == mp0
# FIXME @test mp1 - mp0 == 0

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

@test mpones(Int64, 2) isa Vector{MP{Int64}}
@test mpones(Int64, 2) == [mp1; mp1]
@test mpones(Float64, 2) isa Vector{MP{Float64}}
@test mpones(Float64, 2) == [mp1; mp1]
@test mpones(2) isa Vector{MP{Float64}}
@test mpones(2) == [mp1; mp1]
@test mpones(Int64, 2,5) isa Matrix{MP{Int64}}
@test mpones(Int64, 2,5) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]
@test mpones(Float64, 2,5) isa Matrix{MP{Float64}}
@test mpones(Float64, 2,5) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]
@test mpones(2,5) isa Matrix{MP{Float64}}
@test mpones(2,5) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]

# Identity matrix

@test mpeye(Int64, 2) isa Matrix{MP{Int64}}
@test mpeye(Int64, 2) == [mpone(Int64) mpzero(Int64); mpzero(Int64) mpone(Int64)]
@test mpeye(Float64, 2) isa Matrix{MP{Float64}}
@test mpeye(Float64, 2) == [mp1 mp0; mp0 mp1]
@test mpeye(2) isa Matrix{MP{Float64}}
@test mpeye(2) == [mp1 mp0; mp0 mp1]
@test mpeye(Float64, 2,5) isa Matrix{MP{Float64}}
@test mpeye(Float64, 2,5) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]
@test mpeye(2,5) isa Matrix{MP{Float64}}
@test mpeye(2,5) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]
@test mpeye(2) isa Matrix{MP{Float64}}
@test mpeye(2) == [mp1 mp0; mp0 mp1]

# Matrix of zeros

@test mpzeros(Float64, 2) isa SparseVector{MP{Float64}, Int64}
@test mpzeros(Float64, 2).nzval == MP([])
@test mpzeros(2) isa SparseVector{MP{Float64}, Int64}
@test mpzeros(2).nzval == MP([])
@test mpzeros(Float64, 2,3) isa SparseMatrixCSC{MP{Float64}, Int64}
@test mpzeros(Float64, 2,3).nzval == MP([])
@test mpzeros(2,3) isa SparseMatrixCSC{MP{Float64}, Int64}
@test mpzeros(2,3).nzval == MP([])

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
# Max-Plus matrix multiplication
# ==============================================================================

A = [ϵ ϵ 3; -2 ϵ ϵ; ϵ ϵ 0]
B = [ϵ 2 ϵ; ϵ ϵ 0; -3 ϵ ϵ]
D = [3 ϵ ϵ; ϵ -2 ϵ; ϵ ϵ 0]
P = [ϵ ϵ 0; 0 ϵ ϵ; ϵ ϵ 0]

@test A == (D * P)
@test A == (mpsparse(D) * P)
@test A == (D * mpsparse(P))
@test A == (mpsparse(D) * mpsparse(P))
@test mpsparse(A) == (mpsparse(D) * mpsparse(P))

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
# FIXME @test (A / A) == C

# ==============================================================================
# Max-Plus trace
# ==============================================================================

A = [5 ϵ 5; ϵ 6 3; 11 12 11]

@test mptrace(A) == MP(11.0)
@test mptrace([]) == mp0
@test mptrace(mpeye(Float64, 2,2)) == tr(mpeye(Float64, 2,2)) == 0.0
@test mptrace(mpeye(Float64, 2,5)) == 0.0
@test mptrace(mpzeros(Float64, 2,2)) == mptrace(full(mpzeros(Float64, 2,2))) == tr(mpzeros(Float64, 2,2)) == mp0
@test mptrace(mpzeros(Float64, 2,5)) == mptrace(full(mpzeros(Float64, 2,5))) == mp0
@test mptrace([1.0 2.0; 3.0 4.0]) == tr(MP([1.0 2.0; 3.0 4.0])) == MP(1.0) + MP(4.0) == MP(4.0)

# ==============================================================================
# Max-Plus norm
# ==============================================================================

@test mpnorm(MP([1.0 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399.0)
@test mpnorm(MP([1 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399)
@test mpnorm(mpsparse([1 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399)
@test mpnorm([mp0 1; 10 mp0]) == MP(10.0) - mp0 == mptop

# ==============================================================================
# Max-Plus star
# ==============================================================================
@test mpstar(MP(1.0)) == mptop
@test mpstar(MP(-1.0)) == mp1
@test mpstar(MP([1.0 2; 3 4])) == [mptop mptop; mptop mptop]
A=MP([-3.0 -2; -1 0]); B = mpstar(A)
@test B == mpeye(Float64, 2,2) + A
@test B == B + A^2
@test B == B + A^3

@test mpstar(mpeye(Float64, 2,2)) == mpeye(Float64, 2,2)

# TODO astarb

# ==============================================================================
# Max-Plus power operator
# ==============================================================================

# Scalar

@test MP(2)^4   == MP(2 * 4) == MP(8)
@test MP(0)^0   == MP(0 * 0) == MP(0)
@test MP(2)^0   == MP(2 * 0) == MP(0)
@test MP(2)^-3  == MP(2)^(-3) == MP(2 * -3) == MP(-6)
@test MP(2)^0.5 == MP(2 * 0.5) == MP(1.0)
@test ϵ^0 == MP(0.0)
@test ϵ^2 == ϵ
# FIXME @test ϵ^(-2) == ϵ

# Matrix

A = MP([4 3; 7 -Inf])
@test A^0 == mpeye(Float64, 2,2)
@test A * A == A^2 == MP([10.0 7.0; 11.0 10.0])
@test A * A * A == A^3 == MP([14.0 13.0; 17.0 14.0])
# FIXME A^-1

# ==============================================================================
# Max-Plus other operations
# ==============================================================================

@test abs2(MP(3.0)) == MP(6.0)

@test b / b == MP(3.0 - 3.0) == MP(0.0)
@test b - b == MP(3.0 - 3.0) == MP(0.0)
@test MP(3.0) / MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)
@test MP(3.0) - MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)

# ==============================================================================
# Max-Plus min operator
# ==============================================================================

@test min(MP(3.0), mp0) == mp0
@test min(MP(3.0), mp1) == mp1
@test min(MP(1), MP(2)) == MP(1)
@test min(MP([10 1; 10 1]), MP([4 5; 6 5])) == MP([4 1; 6 1])
@test min(mpsparse([10.0 1; 10 1]), mpsparse([4.0 5; 6 5])) == mpsparse([4.0 1; 6 1])
@test min(mpsparse([10 1; 10 1]), mpsparse([4 5; 6 5])) == mpsparse([4 1; 6 1])
@test min(sparse(mpeye(2,2)), mpzeros(2,2)) == mpzeros(2,2)
@test min(mpeye(2,2), mpones(2,2)) == mpeye(2,2)

@test max(MP(3.0), mp0) == MP(3.0)
@test max(MP(3.0), mp1) == MP(3.0)
@test max(MP(1), MP(2)) == MP(2)
@test max(MP([10 1; 10 1]), MP([4 5; 6 5])) == MP([10 5; 10 5])
@test max(mpsparse([10.0 1; 10 1]), mpsparse([4.0 5; 6 5])) == mpsparse([10.0 5; 10 5])
@test max(mpsparse([10 1; 10 1]), mpsparse([4 5; 6 5])) == MP([10 5; 10 5])
@test max(sparse(mpeye(2,2)), mpzeros(2,2)) == mpeye(2,2)
@test max(mpeye(2,2), mpones(2,2)) == mpones(2,2)

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
result = @capture_out show(stdout, mpzero(Int64))
@test result == "-9223372036854775808"

mp_change_display(1)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, mpzero(Int64))
@test result == "."

mp_change_display(2)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, mpzero(Int64))
@test result == "."

mp_change_display(3)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MP(4))
@test result == "4"
result = @capture_out show(stdout, mpzero(Int64))
@test result == "ε"

mp_change_display(4)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, mpzero(Int64))
@test result == "ε"

A = MP([4.5 0.0; 7.0 -Inf])
result = @capture_out show(stdout, A)
@test result == "MP{Float64}[4.5 e; 7 ε]"
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & e \\\\\n7 & \\varepsilon \\\\\n\\end{array}\n\\right]\n"
mp_change_display(0)
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & 0 \\\\\n7 & -\\infty \\\\\n\\end{array}\n\\right]\n"
