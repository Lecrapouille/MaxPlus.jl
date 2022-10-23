# ==============================================================================
# Unit tests on Min-Plus type: type promotions, scalars, matrices and operations
# ==============================================================================

# ==============================================================================
# Scalar constructor
# ==============================================================================

a = MI(1.0)
@test typeof(a) == MI
@test a isa MI
@test a.λ isa Float64
@test a.λ == 1.0
@test isnan(a) == false
@test isinf(a) == false

b = MI(1)
@test typeof(b) == MI
@test b isa MI
@test b.λ isa Float64
@test b.λ == 1.0
@test isnan(b) == false
@test isinf(b) == false

c = MI(-1)
@test c isa MI
@test c.λ isa Float64
@test c.λ == -1.0
@test isnan(c) == false
@test isinf(c) == false

d = -MI(1)
@test d isa MI
@test d.λ isa Float64
@test d.λ == -1.0
@test isnan(d) == false
@test isinf(d) == false

e = MI(NaN)
@test e isa MI
@test e.λ isa Float64
@test e.λ == Inf
@test isnan(e) == false
@test isinf(e) == true

f = Tropical(MI, NaN)
@test f isa MI
@test f.λ isa Float64
@test f.λ == Inf
@test isnan(f) == false
@test isinf(f) == true

g = Tropical(MI, 42.0)
@test g isa MI
@test g.λ isa Float64
@test g.λ == 42.0
@test isnan(g) == false
@test isinf(g) == false

# ==============================================================================
# Scalar copy constructor
# ==============================================================================

a = MI(MI(1))
@test a isa MI
@test a.λ == 1.0

b = MI(MP(3))
@test b isa MI
@test b.λ == 3.0

# ==============================================================================
# Neutral / absorbing constants
# ==============================================================================

@test mi0 isa MI
@test mi0.λ isa Float64
@test mi0.λ == Inf
@test iszero(mi0) == true
@test isone(mi0) == false

@test mi1 isa MI
@test mi1.λ isa Float64
@test mi1.λ == 0.0
@test iszero(mi1) == false
@test isone(mi1) == true

@test mitop isa MI
@test mitop.λ isa Float64
@test mitop.λ == -Inf
@test iszero(mitop) == false
@test isone(mitop) == false

# ==============================================================================
# Boolean constructor (used for identity matrices)
# ==============================================================================

@test MI(true) == mi1
@test MI(false) == mi0

# ==============================================================================
# Signs
# ==============================================================================

@test sign(MI(0)) == 0
@test sign(MI(5)) == 1
@test sign(MI(-5)) == -1
@test sign(mi1) == 0
@test sign(mi0) == 1

# ==============================================================================
# Scalar comparaisons
# ==============================================================================

# Min-Plus scalars

a = MI(3.0)
@test (a == a) == true
@test (a != a) == (a ≠ a) == false
@test (a >= a) == (a ≥ a) == true
@test (a <= a) == (a ≤ a) == true
@test (a > a) == false
@test (a < a) == false

# Min-Plus vs. non Max-Plus comparaison

@test (a == 3.0) == true
@test (a == 4.0) == false
@test (a != 4.0) == true
@test (a < 4.0) == true
@test (a <= 4.0) == true
@test (a > 4.0) == false
@test (a >= 4.0) == false
@test (4.0 > a) == true
@test (4.0 >= a) == true

# ==============================================================================
# Min-Plus to classic algebra conversion
# ==============================================================================

# Scalar

@test plustimes(MI(2.0)) isa Float64
@test plustimes(MI(2.0)) == 2.0
@test plustimes(MI(-2.0)) isa Float64
@test plustimes(MI(-2.0)) == -2.0

@test float(MI(2.0)) isa Float64
@test float(MI(2.0)) == 2.0
@test float(MI(-2.0)) isa Float64
@test float(MI(-2.0)) == -2.0

@test Float64(MI(2.0)) isa Float64
@test Float64(MI(2.0)) == 2.0
@test Float64(MI(-2.0)) isa Float64
@test Float64(MI(-2.0)) == -2.0

# Array

@test plustimes(MI([2 3; 4 5])) isa Matrix{Float64}
@test plustimes(MI([2 3; 4 5])) == [2.0 3.0; 4.0 5.0]
@test float(MI([2 3; 4 5])) isa Matrix{Float64}
@test float(MI([2 3; 4 5])) == [2.0 3.0; 4.0 5.0]

# ==============================================================================
# zero / one elements
# ==============================================================================

@test zero(MI) == MI(Inf) == Inf
@test zero(MI) == mi0 == Inf
@test zero(MI(42.0)) == mi0 == Inf
@test zero(MI(42)) == mi0 == Inf
@test zero(MI) == MI(Inf) == Inf
@test zero(eltype(MI([1 2]))) == mi0

@test one(MI) == MI(0.0) == 0.0
@test one(MI) == mi1 == 0.0
@test one(MI(42.0)) == mi1 == 0.0
@test one(MI) == mi1 == 0.0
@test one(MI(42)) == mi1 == 0
@test one(eltype(MI([1 2]))) == mi1

# ==============================================================================
# Operations on scalars and neutral elements and commutativity
# ==============================================================================

b = MI(3.0)
@test b + mi0 == mi0 + b == MI(min(3.0, Inf)) == MI(min(Inf, 3.0)) == b == MI(3.0)
@test b * mi1 == mi1 * b == MI(3.0 + 0.0) == MI(0.0 + 3.0) == b == MI(3.0)

@test b * b == MI(3.0 + 3.0) == MI(6.0)
@test b + b == MI(min(3.0, 3.0)) == MI(3.0)

a = MI(5.0)
b = MI(3.0)
@test a + b == MI(min(5.0, 3.0)) == b + a == MI(min(3.0, 5.0)) == MI(3.0)
@test a * b == MI(5.0 + 3.0) == b * a == MI(3.0 + 5.0) == MI(8.0)

# ==============================================================================
# Operations on scalars of different types
# ==============================================================================

b = MI(3.0)
@test b + 3.0 == 3.0 + b == b == MI(min(3.0, 3.0)) == MI(3.0)
@test b + 3  == 3 + b == b == MI(min(3.0, 3.0)) == MI(3.0)
@test b * 3.0 == 3.0 * b == MI(3.0 + 3.0) == MI(6.0)
@test b * 3  == 3 * b == MI(3.0 + 3.0) == MI(6.0)

# ==============================================================================
#  Distributivity of operations
# ==============================================================================

a = MI(5.0)
b = MI(3.0)
c = MI(1.0)
@test (a + b) + c == a + (b + c) == a + b + c == MI(min(1.0, 3.0, 5.0)) == MI(1.0)
@test (a * b) * c == a * (b * c) == a * b * c == MI(1.0 + 3.0 + 5.0) == MI(9.0)
@test (a + b) * c == MI(min(5.0, 3.0) + 1.0) == MI(4.0)
@test (a * c) + (b * c) == MI(min(5.0 + 1.0, 3.0 + 1.0)) == MI(4.0)

# ==============================================================================
# Residuation: not implemented
# ==============================================================================

# ==============================================================================
# Forbiden minus operator
# ==============================================================================

@test_throws ErrorException("Minus operator does not exist in (min,+) algebra") MI(3) - MI(3)
@test_throws ErrorException("Minus operator does not exist in (min,+) algebra") MI(3) - 3
@test_throws ErrorException("Minus operator does not exist in (min,+) algebra") 3 - MI(3)

# ==============================================================================
# Mixture of Max-Plus and Min-Plus is forbidden
# ==============================================================================

@test_throws ErrorException("Cannot promote (min,+) to (max,+)") MI(3) + MP(3)
@test_throws ErrorException("Cannot promote (min,+) to (max,+)") MI(3) * MP(3)

# ==============================================================================
# Scalar power operator: not implemented
# ==============================================================================

# ==============================================================================
# Max-Plus min operator on scalars and on matrices
# ==============================================================================

@test min(MP(3.0), mp0) == min(mp0, MP(3.0)) == mp0
@test min(MP(3.0), mp1) == min(mp1, MP(3.0)) == mp1
@test min(MP(1), MP(2)) == min(1, MP(2)) == min(MP(1), 2) == MP(1)
@test min(MP([10 1; 10 1]), MP([4 5; 6 5])) == MP([4 4; 4 4])
@test min(sparse(MP([10.0 1; 10 1])), sparse(MP([4.0 5; 6 5]))) == sparse(MP([4 4; 4 4]))
@test min(MP([10 1; 10 1]), sparse(MP([4 5; 6 5]))) == sparse(MP([4 4; 4 4]))
@test min(sparse(eye(MP,2,2)), spzeros(MP,2,2)) == spzeros(MP,2,2)
@test min(eye(MP,2,2), ones(MP,2,2)) == ones(MP,2,2)

# ==============================================================================
# Border cases operations
# ==============================================================================

@test mi0 + mi0 == mi0
@test mi0 + mi1 == 0
@test mi0 + mitop == mitop
@test mi1 + mi0 == mi1
@test mi1 + mi1 == mi1
@test mi1 + mitop == mitop
@test mitop + mi0 == mitop
@test mitop + mi1 == mitop
@test mitop + mitop == mitop

@test mi0 * mi0 == mi0
@test mi0 * mi1 == mi0
@test_broken mi0 * mitop == mitop
@test mi1 * mi0 == mi0
@test mi1 * mi1 == mi1
@test mi1 * mitop == mitop
@test_broken mitop * mi0 == mitop
@test mitop * mi1 == mitop
@test mitop * mitop == mitop

@test_broken mi0 / mi0 == mi0
@test mi0 / mi1 == mi0
@test mi0 / mitop == mi0
@test_broken mi1 / mi0 == mi1
@test mi1 / mi1 == mi1
@test_broken mi1 / mitop == mi1
@test_broken mitop / mi0 == mi1
@test_broken mitop / mi1 == mi1
@test_broken mitop / mitop == mi1

@test_broken mi0 \ mi0 == mi0
@test_broken mi0 \ mi1 == mi0
@test_broken mi0 \ mitop == mi0
@test_broken mi1 \ mi0 == mi1
@test mi1 \ mi1 == mi1
@test_broken mi1 \ mitop == mi1
@test_broken mitop \ mi0 == mi1
@test_broken mitop \ mi1 == mi1
@test_broken mitop \ mitop == mi1

# ==============================================================================
# Other operations
# ==============================================================================

@test abs2(MI(3.0)) == MI(6.0)
@test abs(MI(-3.0)) == MI(3.0)
@test abs(mi0) == mi0
@test float(MI(2)) == MI(2)
@test round(MI(1.7)) == MI(round(1.7)) == round(1.7)
@test floor(MI(1.7)) == MI(floor(1.7)) == floor(1.7)
@test ceil(MI(1.7)) == MI(ceil(1.7)) == ceil(1.7)
@test trunc(MI(1.7)) == MI(trunc(1.7)) == trunc(1.7)

# ==============================================================================
# Matrix comparaisons
# ==============================================================================

# Dense Min-Plus matrix comparaison

B = [MI(1) MI(2); MI(3) MI(4)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# Sparse Min-Plus matrix comparaison

S = sparse([MI(1) MI(2); MI(3) MI(4)])
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# Sparse/Dense Min-Plus matrix comparaison

@test (S == B) == (B == S) == true

# Sort is using isless()
v = MI([3, 1, 2]);
@test sort(v) == MI([1, 2, 3])

# ==============================================================================
# Vector comparaisons
# ==============================================================================

# Dense vector comparaison

B = [MI(1); MI(3)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true; true]
@test (B .<= B) == (B .≤ B) == [true; true]
@test (B .> B) == [false; false]
@test (B .< B) == [false; false]

# Sparse vector comparaison

S = sparse([MI(1); MI(3)])
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true; true]
@test (B .<= B) == (B .≤ B) == [true; true]
@test (B .> B) == [false; false]
@test (B .< B) == [false; false]

# Sparse/Dense matrix comparaison

@test (S == B) == (B == S) == true

# Sort is using isless()
v = MI([3, 1, 2]);
@test sort(v) == MI([1, 2, 3])

# ==============================================================================
# Range construction
# ==============================================================================

A = MI(1:3)
@test A == [MI(1); MI(2); MI(3)]

B = MI(1.0:0.5:3.0)
@test B == [MI(1.0); MI(1.5); MI(2.0); MI(2.5); MI(3.0)]

# ==============================================================================
# Dense matrix constructor and type promotion/contamination
# ==============================================================================

A = [MI(1) MI(2); MI(3) MI(4)]
@test A isa Matrix{MI}

B = MI([1 2; 3 4])
@test B isa Matrix{MI}

C = [MI(1) 2; 3 4]
@test C isa Matrix{MI}

D = [MI(1.0) 2; 3 4]
@test D isa Matrix{MI}

E = [MI(1) 2.0; 3 4]
@test E isa Matrix{MI}

F = [1 MI(2.0); 3 4]
@test F isa Matrix{MI}

@test A == B == C == D == E == F

# ==============================================================================
# Sparse constructor
# ==============================================================================

# Min-Plus: Using SparseArray sparse(MI)

sA = sparse([MI(1) MI(2); MI(3) MI(4)])
@test sA isa SparseMatrixCSC{MI, Int64}

sB = MI(sparse([1 2; 3 4]))
@test sB isa SparseMatrixCSC{MI, Int64}

sC = sparse([MI(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MI, Int64}

sD = sparse([MI(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MI, Int64}

sE = sparse([MI(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MI, Int64}

sF = sparse([1 MI(2.0); 3 4])
@test sF isa SparseMatrixCSC{MI, Int64}

sG = MI([1; 2; 1; 2], [1; 1; 2; 2], [1; 3; 2; 4])
@test sG isa SparseMatrixCSC{MI, Int64}

sH = MI([1, 2, 1, 2], [1, 1, 2, 2], [1, 3, 2, 4])
@test sH isa SparseMatrixCSC{MI, Int64}

@test sA == sB == sC == sD == sE == sF == sG == sH

# Using SparseArray.sparse (I, J, D) vectors

spA = MI(sparse([1, 2, 3], [1, 2, 3], [Inf, 2, 0]))
@test size(spA.nzval,1) == 3
@test spA.nzval == [mi0; 2.0; 0.0]

spB = MI(dropzeros(sparse([1, 2, 3], [1, 2, 3], [Inf, 2, 0])))
@test size(spB.nzval,1) == 2
@test spB.nzval == [mi0; 2.0]

spC = MI([1, 2, 3], [1, 2, 3], [Inf, 2, 0])
@test size(spC.nzval,1) == 2
@test spC.nzval == [2.0; 0.0]
@test spC != spB

spC = MI(spA)
@test size(spC.nzval,1) == 3
@test spC.nzval == [mi0; 2.0; 0.0]
@test spC == spA

spD = sparse(spB)
@test size(spD.nzval,1) == 2
@test spD.nzval == [mi0; 2.0]

spE = MI([1, 2, 3], [1, 2, 3], MI([Inf, 2, 0]))
@test size(spE.nzval,1) == 2
@test spE.nzval == [2.0; 0.0]

# Using MI(sparse)

spA = sparse(MI([Inf 0; 0 Inf]))
@test findnz(spA) == ([2, 1], [1, 2], MI([0.0, 0.0]))
spB = MI(sparse([Inf 0; 0 Inf]))
@test findnz(spB) == ([1, 2], [1, 2], MI([mi0, mi0]))
@test spA != spB
spC = sparse(MI([4 0; 7 Inf]))
@test findnz(spC) == ([1, 2, 1], [1, 1, 2], MI([4, 7, 0]))

# ==============================================================================
# Matrix ones, eye, zeros constructions
# ==============================================================================

# Min-Plus Matrix of ones

@test ones(MI, 2) isa Vector{MI}
@test ones(MI, 2) == [mi1; mi1]
@test ones(MI, 2,5) isa Matrix{MI}
@test ones(MI, 2,5) == [mi1 mi1 mi1 mi1 mi1; mi1 mi1 mi1 mi1 mi1]
@test ones(MI([1 2; 3 4])) isa Matrix{MI}
@test ones(MI([1 2; 3 4])) == [mi1 mi1; mi1 mi1]

# Min-Plus Identity dense matrix

@test eye(MI, 2) isa Matrix{MI}
@test eye(MI, 2) == [mi1 mi0; mi0 mi1]
@test eye(MI, 2,5) isa Matrix{MI}
@test eye(MI, 2,5) == [mi1 mi0 mi0 mi0 mi0; mi0 mi1 mi0 mi0 mi0]
@test eye(MI([1 2; 3 4])) isa Matrix{MI}
@test eye(MI([1 2; 3 4])) == [mi1 mi0; mi0 mi1]

# Min-Plus Identity sparse matrix

@test speye(MI, 2) isa SparseMatrixCSC{MI}
@test speye(MI, 2) == sparse([mi1 mi0; mi0 mi1])
@test speye(MI, 2,5) isa SparseMatrixCSC{MI}
@test speye(MI, 2,5) == sparse([mi1 mi0 mi0 mi0 mi0; mi0 mi1 mi0 mi0 mi0])
@test speye(MI([1 2; 3 4])) isa SparseMatrixCSC{MI}
@test speye(MI([1 2; 3 4])) == sparse([mi1 mi0; mi0 mi1])

# Min-Plus Matrix of zeros

@test spzeros(MI, 2) isa SparseVector{MI, Int64}
@test spzeros(MI, 2).nzval == MI([])
@test spzeros(MI, 2,3) isa SparseMatrixCSC{MI, Int64}
@test spzeros(MI, 2,3).nzval == MI([])
#FIXME broken with Julia 1.6
#@test spzeros(MI([1 2; 3 4])) isa SparseMatrixCSC{MIMP, Int64}
#@test spzeros(MI([1 2; 3 4])).nzval == MI([])

@test zeros(MI, 2) isa Vector{MI}
@test zeros(MI, 2) == MI([mi0; mi0])
@test zeros(MI, 2,3) isa Matrix{MI}
@test zeros(MI, 2,3) == MI([mi0 mi0 mi0; mi0 mi0 mi0])
@test zeros(MI([1 2; 3 4])) isa Matrix{MI}
@test zeros(MI([1 2; 3 4])) == [mi0 mi0; mi0 mi0]

# ==============================================================================
# Matrix dense/sparse multiplication
# ==============================================================================

A = [mi0 mi0 3; -2 mi0 mi0; mi0 mi0 0]
B = [mi0 2 mi0; mi0 mi0 0; -3 mi0 mi0]
D = [3 mi0 mi0; mi0 -2 mi0; mi0 mi0 0]
P = [mi0 mi0 0; 0 mi0 mi0; mi0 mi0 0]

@test A == (D * P)
@test A == (sparse(D) * P)
@test A == (D * sparse(P))
@test A == (sparse(D) * sparse(P))
@test sparse(A) == (sparse(D) * sparse(P))

# ==============================================================================
# Matrix addition
# ==============================================================================

A = MP([3.0 4; 5 6])
B = MP([1.0 2; 3 4])
C = MP([3.0 4; 5 6])
@test (A + B) == (B + A) == C

# ==============================================================================
# Matrix ones, eye, zeros operations
# ==============================================================================

I = eye(MI,4,4)
Z = spzeros(MI,4,4)
O = ones(MI,4,4)
A = MI(rand(4,4))

@test (A * I) == (I * A) == A
@test (A + Z) == (Z + A) == A
@test (A * Z) == (Z * A) == Z
@test (A + O) == (O + A) == O
@test (A * O) != (O * A) != A

# ==============================================================================
# Max-Plus element by element operations on matrices
# ==============================================================================

# Dense matrix

A = MI([1.0 2.0; 3.0 4.0])
@test 2.0 .+ A == A .+ 2.0 == [MI(2.0) + MI(1.0) MI(2.0) + MI(2.0); MI(2.0) + MI(3.0) MI(2.0) + MI(4.0)] == MI([1.0 2.0; 2.0 2.0])
@test 2.0 .* A == A .* 2.0 == [MI(2.0) * MI(1.0) MI(2.0) * MI(2.0); MI(2.0) * MI(3.0) MI(2.0) * MI(4.0)] == MI([3.0 4.0; 5.0 6.0])

# Sparse matrix

sA = sparse(A)
@test 2.0 .+ sA == sA .+ 2.0 == [MI(2.0) + MI(1.0) MI(2.0) + MI(2.0); MI(2.0) + MI(3.0) MI(2.0) + MI(4.0)] == MI([1.0 2.0; 2.0 2.0])
@test 2.0 .* sA == sA .* 2.0 == [MI(2.0) * MI(1.0) MI(2.0) * MI(2.0); MI(2.0) * MI(3.0) MI(2.0) * MI(4.0)] == MI([3.0 4.0; 5.0 6.0])

# ==============================================================================
# Min-Plus type promotion/contamination on dense/sparse matrix
# ==============================================================================

# Dense matrix

A = [1 2; 3 4]
@test A isa Matrix{Int64}
@test MI(A) isa Matrix{MI}

B = [MI(1) MI(2); MI(3) MI(4)]
@test B isa Matrix{MI}

C = [MI(1) 2; 3 4]
@test C isa Matrix{MI}

D = [MI(1.0) 2; 3 4]
@test D isa Matrix{MI}

E = [MI(1) 2.0; 3 4]
@test E isa Matrix{MI}

F = [1 MI(2.0); 3 4]
@test F isa Matrix{MI}

# Sparse matrix

sA = sparse([1 2; 3 4])
@test sA isa SparseMatrixCSC{Int64, Int64}
@test MI(sA) isa SparseMatrixCSC{MI, Int64}

sB = sparse([MI(1) MI(2); MI(3) MI(4)])
@test sB isa SparseMatrixCSC{MI, Int64}

sC = sparse([MI(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MI, Int64}

sD = sparse([MI(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MI, Int64}

sE = sparse([MI(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MI, Int64}

sF = sparse([1 MI(2.0); 3 4])
@test sF isa SparseMatrixCSC{MI, Int64}

# Dense/Sparse matrix comparaison

@test (A == sA) == (sA == A)
@test (B == sB) == (sB == B)
@test (C == sC) == (sC == C)
@test (D == sD) == (sD == D)
@test (E == sE) == (sE == E)
@test (F == sF) == (sF == F)

# ==============================================================================
# Max-Plus scalar to class algebra scalar conversion
# ==============================================================================

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

@test full(sC) isa Matrix{MI}
@test full(sD) isa Matrix{MI}
@test full(sE) isa Matrix{MI}
@test full(sF) isa Matrix{MI}

@test full(sC) == C
@test full(sD) == D
@test full(sE) == E
@test full(sF) == F

# dense() is a anlias for full()

@test dense(sC) isa Matrix{MI}
@test dense(sD) isa Matrix{MI}
@test dense(sE) isa Matrix{MI}
@test dense(sF) isa Matrix{MI}

@test dense(sC) == C
@test dense(sD) == D
@test dense(sE) == E
@test dense(sF) == F

# Max-Plus sparse array to Max-Plus dense array

Z = dense(spzeros(MI,2,2))
@test typeof(Z) == Matrix{MI}
@test Z == [mi0 mi0; mi0 mi0]

# Max-Plus sparse array to Max-Plus dense array

Z = full(spzeros(MI,2,2))
@test typeof(Z) == Matrix{MI}
@test Z == [mi0 mi0; mi0 mi0]

# ==============================================================================
# Bug with Julia with SparseMatrixCSC and operator== which confused zero() and 0.0
# ==============================================================================

A = MI(sparse([1, 2], [1, 2], [0.0, 0.0]))
B = spzeros(MI,2,2)
@test A.nzval == MI([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mi0, mi0])
BB = spzeros(MI,2,2)
@test AA.nzval == MI([Inf, Inf])
@test (AA == BB) == true

# ==============================================================================
# Display
# ==============================================================================

using Suppressor

set_tropical_display(0)
result = @capture_out show(stdout, mi0)
@test result == "Inf"
result = @capture_out show(stdout, mi1)
@test result == "0.0"
result = @capture_out show(stdout, MI(4.0))
@test result == "4.0"
result = @capture_out show(stdout, zero(MI))
@test result == "Inf"
result = @capture_out show(stdout, one(MI))
@test result == "0.0"

set_tropical_display(1)
result = @capture_out show(stdout, mi0)
@test result == "."
result = @capture_out show(stdout, mi1)
@test result == "0"
result = @capture_out show(stdout, MI(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MI))
@test result == "."
result = @capture_out show(stdout, one(MI))
@test result == "0"

set_tropical_display(2)
result = @capture_out show(stdout, mi0)
@test result == "."
result = @capture_out show(stdout, mi1)
@test result == "e"
result = @capture_out show(stdout, MI(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MI))
@test result == "."
result = @capture_out show(stdout, one(MI))
@test result == "e"

set_tropical_display(3)
result = @capture_out show(stdout, mi0)
@test result == "ε"
result = @capture_out show(stdout, mi1)
@test result == "0"
result = @capture_out show(stdout, MI(4))
@test result == "4"
result = @capture_out show(stdout, zero(MI))
@test result == "ε"
result = @capture_out show(stdout, one(MI))
@test result == "0"

set_tropical_display(4)
result = @capture_out show(stdout, mi0)
@test result == "ε"
result = @capture_out show(stdout, mi1)
@test result == "e"
result = @capture_out show(stdout, MI(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MI))
@test result == "ε"
result = @capture_out show(stdout, one(MI))
@test result == "e"

A = MI([4.5 0.0; 7.0 Inf])
result = @capture_out show(stdout, A)
@test result == "MI[4.5 e; 7 ε]"
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & e \\\\\n7 & \\varepsilon \\\\\n\\end{array}\n\\right]\n"
set_tropical_display(0)
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & 0 \\\\\n7 & +\\infty \\\\\n\\end{array}\n\\right]\n"
