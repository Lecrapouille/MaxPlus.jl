# ==============================================================================
# Unit tests on Max-Plus type: type promotions, scalars, matrices and operations
# ==============================================================================

# ==============================================================================
# Scalar constructor
# ==============================================================================

a = MP(1.0)
@test typeof(a) == MP
@test a isa MP
@test a.λ isa Float64
@test a.λ == 1.0
@test isnan(a) == false
@test isinf(a) == false

b = MP(1)
@test typeof(b) == MP
@test b isa MP
@test b.λ isa Float64
@test b.λ == 1.0
@test isnan(b) == false
@test isinf(b) == false

c = MP(-1)
@test c isa MP
@test c.λ isa Float64
@test c.λ == -1.0
@test isnan(c) == false
@test isinf(c) == false

d = -MP(1)
@test d isa MP
@test d.λ isa Float64
@test d.λ == -1.0
@test isnan(d) == false
@test isinf(d) == false

e = MP(NaN)
@test e isa MP
@test e.λ isa Float64
@test e.λ == -Inf
@test isnan(e) == false
@test isinf(e) == true

f = Tropical(MP, NaN)
@test f isa MP
@test f.λ isa Float64
@test f.λ == -Inf
@test isnan(f) == false
@test isinf(f) == true

g = Tropical(MP, 42.0)
@test g isa MP
@test g.λ isa Float64
@test g.λ == 42.0
@test isnan(g) == false
@test isinf(g) == false

# ==============================================================================
# Scalar copy constructor
# ==============================================================================

a = MP(MP(1))
@test a isa MP
@test a.λ == 1.0

b = MP(MI(3))
@test b isa MP
@test b.λ == 3.0

# ==============================================================================
# Neutral / absorbing constants
# ==============================================================================

@test mp0 isa MP
@test mp0.λ isa Float64
@test mp0.λ == -Inf
@test iszero(mp0) == true
@test isone(mp0) == false

@test mp1 isa MP
@test mp1.λ isa Float64
@test mp1.λ == 0.0
@test iszero(mp1) == false
@test isone(mp1) == true

@test mptop isa MP
@test mptop.λ isa Float64
@test mptop.λ == Inf
@test iszero(mptop) == false
@test isone(mptop) == false

# ==============================================================================
# Boolean constructor (used for identity matrices)
# ==============================================================================

@test MP(true) == mp1
@test MP(false) == mp0

# ==============================================================================
# Signs
# ==============================================================================

@test sign(MP(0)) == 0
@test sign(MP(5)) == 1
@test sign(MP(-5)) == -1
@test sign(mp1) == 0
@test sign(mp0) == -1

# ==============================================================================
# Scalar comparaisons
# ==============================================================================

# Max-Plus scalars

a = MP(3.0)
@test (a == a) == true
@test (a != a) == (a ≠ a) == false
@test (a >= a) == (a ≥ a) == true
@test (a <= a) == (a ≤ a) == true
@test (a > a) == false
@test (a < a) == false

# Max-Plus vs. non Max-Plus comparaison

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
# Max-Plus to classic algebra conversion
# ==============================================================================

# Scalar

@test plustimes(MP(2.0)) isa Float64
@test plustimes(MP(2.0)) == 2.0
@test plustimes(MP(-2.0)) isa Float64
@test plustimes(MP(-2.0)) == -2.0

@test float(MP(2.0)) isa Float64
@test float(MP(2.0)) == 2.0
@test float(MP(-2.0)) isa Float64
@test float(MP(-2.0)) == -2.0

@test Float64(MP(2.0)) isa Float64
@test Float64(MP(2.0)) == 2.0
@test Float64(MP(-2.0)) isa Float64
@test Float64(MP(-2.0)) == -2.0

# Array

@test plustimes(MP([2 3; 4 5])) isa Matrix{Float64}
@test plustimes(MP([2 3; 4 5])) == [2.0 3.0; 4.0 5.0]
@test float(MP([2 3; 4 5])) isa Matrix{Float64}
@test float(MP([2 3; 4 5])) == [2.0 3.0; 4.0 5.0]

# ==============================================================================
# zero / one elements
# ==============================================================================

@test zero(MP) == MP(-Inf) == -Inf
@test zero(MP) == mp0 == -Inf
@test zero(MP(42.0)) == mp0 == -Inf
@test zero(MP(42)) == mp0 == -Inf
@test zero(MP) == MP(-Inf) == -Inf
@test zero(eltype(MP([1 2]))) == mp0

@test one(MP) == MP(0.0) == 0.0
@test one(MP) == mp1 == 0.0
@test one(MP(42.0)) == mp1 == 0.0
@test one(MP) == mp1 == 0.0
@test one(MP(42)) == mp1 == 0
@test one(eltype(MP([1 2]))) == mp1

# ==============================================================================
# Operations on scalars and neutral elements and commutativity
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
# Operations on scalars of different types
# ==============================================================================

b = MP(3.0)
@test b + 3.0 == 3.0 + b == b == MP(max(3.0, 3.0)) == MP(3.0)
@test b + 3  == 3 + b == b == MP(max(3.0, 3.0)) == MP(3.0)
@test b * 3.0 == 3.0 * b == MP(3.0 + 3.0) == MP(6.0)
@test b * 3  == 3 * b == MP(3.0 + 3.0) == MP(6.0)

# ==============================================================================
#  Distributivity of operations
# ==============================================================================

a = MP(5.0)
b = MP(3.0)
c = MP(1.0)
@test (a + b) + c == a + (b + c) == a + b + c == MP(max(1.0, 3.0, 5.0)) == MP(5.0)
@test (a * b) * c == a * (b * c) == a * b * c == MP(1.0 + 3.0 + 5.0) == MP(9.0)
@test (a + b) * c == MP(max(5.0, 3.0) + 1.0) == MP(6.0)
@test (a * c) + (b * c) == MP(max(5.0 + 1.0, 3.0 + 1.0)) == MP(6.0)

# ==============================================================================
# Residuation
# ==============================================================================

@test MP(3) / MP(3) == -(MP(3) * MP(-3)) == MP(0.0)
@test MP(3) / MP(5) == -(MP(5.0) * MP(-3.0)) == MP(-2.0)
@test MP(3) \ MP(3) == -(MP(3) * MP(-3)) == MP(0.0)
@test MP(3) \ MP(6) == -(MP(3.0) * MP(-6.0)) == MP(3.0)
@test MP(3) \ mp0 == mp0
@test MP(3) \ mp1 == MP(-3)

# ==============================================================================
# Forbiden minus operator
# ==============================================================================

@test_throws ErrorException("Minus operator does not exist in (max,+) algebra") MP(3) - MP(3)
@test_throws ErrorException("Minus operator does not exist in (max,+) algebra") MP(3) - 3
@test_throws ErrorException("Minus operator does not exist in (max,+) algebra") 3 - MP(3)

# ==============================================================================
# Mixture of Max-Plus and Min-Plus is forbidden
# ==============================================================================

@test_throws ErrorException("Cannot promote (max,+) to (min,+)") MP(3) + MI(3)
@test_throws ErrorException("Cannot promote (max,+) to (min,+)") MP(3) * MI(3)

# ==============================================================================
# Scalar power operator
# ==============================================================================

@test MP(2)^4  == MP(2 * 4) == MP(8)
@test MP(0)^0  == MP(0 * 0) == MP(0)
@test MP(2)^0  == MP(2 * 0) == MP(0)
@test MP(2)^-3 == MP(2)^(-3) == MP(2 * -3) == MP(-6)
@test MP(2)^0.5 == MP(2 * 0.5) == MP(1.0)
@test MP(2)^(-0.5) == MP(2 * -0.5) == MP(-1.0)
@test mp0^0 == MP(0.0)
@test mp0^2 == mp0
@test mp0^(-2) == mptop
@test mp0^(0.5) == mp0
@test mp0^(-0.5) == mptop
@test mp1^0 == mp1
@test mp1^2 == mp1
@test mp1^(-2) == mp1
@test mp1^(0.5) == mp1
@test mp1^(-0.5) == mp1
@test inv(MP(5)) == MP(5)^-1 == MP(-5)

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

@test mp0 + mp0 == mp0
@test mp0 + mp1 == 0
@test mp0 + mptop == mptop
@test mp1 + mp0 == mp1
@test mp1 + mp1 == mp1
@test mp1 + mptop == mptop
@test mptop + mp0 == mptop
@test mptop + mp1 == mptop
@test mptop + mptop == mptop

@test mp0 * mp0 == mp0
@test mp0 * mp1 == mp0
@test mp0 * mptop == mp0
@test mp1 * mp0 == mp0
@test mp1 * mp1 == mp1
@test mp1 * mptop == mptop
@test mptop * mp0 == mp0
@test mptop * mp1 == mptop
@test mptop * mptop == mptop

@test mp0 / mp0 == mptop
@test mp0 / mp1 == mp0
@test mp0 / mptop == mp0
@test mp1 / mp0 == mptop
@test mp1 / mp1 == mp1
@test mp1 / mptop == mp0
@test mptop / mp0 == mptop
@test mptop / mp1 == mptop
@test mptop / mptop == mptop

@test mp0 \ mp0 == mptop
@test mp0 \ mp1 == mptop
@test mp0 \ mptop == mptop
@test mp1 \ mp0 == mp0
@test mp1 \ mp1 == mp1
@test mp1 \ mptop == mptop
@test mptop \ mp0 == mp0
@test mptop \ mp1 == mp0
@test mptop \ mptop == mptop

# ==============================================================================
# Other operations
# ==============================================================================

@test abs2(MP(3.0)) == MP(6.0)
@test abs(MP(-3.0)) == MP(3.0)
@test abs(mp0) == mptop
@test float(MP(2)) == MP(2)
@test round(MP(1.7)) == MP(round(1.7)) == round(1.7)
@test floor(MP(1.7)) == MP(floor(1.7)) == floor(1.7)
@test ceil(MP(1.7)) == MP(ceil(1.7)) == ceil(1.7)
@test trunc(MP(1.7)) == MP(trunc(1.7)) == trunc(1.7)

# ==============================================================================
# Matrix comparaisons
# ==============================================================================

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
# Vector comparaisons
# ==============================================================================

# Dense Max-Plus vector comparaison

B = [MP(1); MP(3)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true; true]
@test (B .<= B) == (B .≤ B) == [true; true]
@test (B .> B) == [false; false]
@test (B .< B) == [false; false]

# Sparse Max-Plus vector comparaison

S = sparse([MP(1); MP(3)])
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true; true]
@test (B .<= B) == (B .≤ B) == [true; true]
@test (B .> B) == [false; false]
@test (B .< B) == [false; false]

# Sparse/Dense Max-Plus matrix comparaison

@test (S == B) == (B == S) == true

# Sort is using isless()
v = MP([3, 1, 2]);
@test sort(v) == MP([1, 2, 3])

# ==============================================================================
# Range construction
# ==============================================================================

A = MP(1:3)
@test A == [MP(1); MP(2); MP(3)]

B = MP(1.0:0.5:3.0)
@test B == [MP(1.0); MP(1.5); MP(2.0); MP(2.5); MP(3.0)]

# ==============================================================================
# Dense matrix constructor and type promotion/contamination
# ==============================================================================

A = [MP(1) MP(2); MP(3) MP(4)]
@test A isa Matrix{MP}

B = MP([1 2; 3 4])
@test B isa Matrix{MP}

C = [MP(1) 2; 3 4]
@test C isa Matrix{MP}

D = [MP(1.0) 2; 3 4]
@test D isa Matrix{MP}

E = [MP(1) 2.0; 3 4]
@test E isa Matrix{MP}

F = [1 MP(2.0); 3 4]
@test F isa Matrix{MP}

@test A == B == C == D == E == F

# ==============================================================================
# Sparse constructor
# ==============================================================================

# Max-Plus: Using SparseArray sparse(MP)

sA = sparse([MP(1) MP(2); MP(3) MP(4)])
@test sA isa SparseMatrixCSC{MP, Int64}

sB = MP(sparse([1 2; 3 4]))
@test sB isa SparseMatrixCSC{MP, Int64}

sC = sparse([MP(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MP, Int64}

sD = sparse([MP(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MP, Int64}

sE = sparse([MP(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MP, Int64}

sF = sparse([1 MP(2.0); 3 4])
@test sF isa SparseMatrixCSC{MP, Int64}

sG = MP([1; 2; 1; 2], [1; 1; 2; 2], [1; 3; 2; 4])
@test sG isa SparseMatrixCSC{MP, Int64}

sH = MP([1, 2, 1, 2], [1, 1, 2, 2], [1, 3, 2, 4])
@test sH isa SparseMatrixCSC{MP, Int64}

@test sA == sB == sC == sD == sE == sF == sG == sH

# Using SparseArray.sparse (I, J, D) vectors

spA = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test size(spA.nzval,1) == 3
@test spA.nzval == [mp0; 2.0; 0.0]

spB = MP(dropzeros(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0])))
@test size(spB.nzval,1) == 2
@test spB.nzval == [mp0; 2.0]

spC = MP([1, 2, 3], [1, 2, 3], [-Inf, 2, 0])
@test size(spC.nzval,1) == 2
@test spC.nzval == [2.0; 0.0]
@test spC != spB

spC = MP(spA)
@test size(spC.nzval,1) == 3
@test spC.nzval == [mp0; 2.0; 0.0]
@test spC == spA

spD = sparse(spB)
@test size(spD.nzval,1) == 2
@test spD.nzval == [mp0; 2.0]

spE = MP([1, 2, 3], [1, 2, 3], MP([-Inf, 2, 0]))
@test size(spE.nzval,1) == 2
@test spE.nzval == [2.0; 0.0]

# Using MP(sparse)

spA = sparse(MP([-Inf 0; 0 -Inf]))
@test findnz(spA) == ([2, 1], [1, 2], MP([0.0, 0.0]))
spB = MP(sparse([-Inf 0; 0 -Inf]))
@test findnz(spB) == ([1, 2], [1, 2], MP([mp0, mp0]))
@test spA != spB
spC = sparse(MP([4 0; 7 -Inf]))
@test findnz(spC) == ([1, 2, 1], [1, 1, 2], MP([4, 7, 0]))

# ==============================================================================
# Matrix ones, eye, zeros constructions
# ==============================================================================

# Max-Plus Matrix of ones

@test ones(MP, 2) isa Vector{MP}
@test ones(MP, 2) == [mp1; mp1]
@test ones(MP, 2,5) isa Matrix{MP}
@test ones(MP, 2,5) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]
@test ones(MP([1 2; 3 4])) isa Matrix{MP}
@test ones(MP([1 2; 3 4])) == [mp1 mp1; mp1 mp1]

# Max-Plus Identity matrix

@test eye(MP, 2) isa Matrix{MP}
@test eye(MP, 2) == [mp1 mp0; mp0 mp1]
@test eye(MP, 2,5) isa Matrix{MP}
@test eye(MP, 2,5) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]
@test eye(MP([1 2; 3 4])) isa Matrix{MP}
@test eye(MP([1 2; 3 4])) == [mp1 mp0; mp0 mp1]

# Max-Plus Identity sparse matrix

@test speye(MP, 2) isa SparseMatrixCSC{MP}
@test speye(MP, 2) == sparse([mp1 mp0; mp0 mp1])
@test speye(MP, 2,5) isa SparseMatrixCSC{MP}
@test speye(MP, 2,5) == sparse([mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0])
@test speye(MP([1 2; 3 4])) isa SparseMatrixCSC{MP}
@test speye(MP([1 2; 3 4])) == sparse([mp1 mp0; mp0 mp1])

# Max-Plus Matrix of zeros

@test spzeros(MP, 2) isa SparseVector{MP, Int64}
@test spzeros(MP, 2).nzval == MP([])
@test spzeros(MP, 2,3) isa SparseMatrixCSC{MP, Int64}
@test spzeros(MP, 2,3).nzval == MP([])
#FIXME broken with Julia 1.6
#@test spzeros(MP([1 2; 3 4])) isa SparseMatrixCSC{MP, Int64}
#@test spzeros(MP([1 2; 3 4])).nzval == MP([])

@test zeros(MP, 2) isa Vector{MP}
@test zeros(MP, 2) == MP([mp0; mp0])
@test zeros(MP, 2,3) isa Matrix{MP}
@test zeros(MP, 2,3) == MP([mp0 mp0 mp0; mp0 mp0 mp0])
@test zeros(MP([1 2; 3 4])) isa Matrix{MP}
@test zeros(MP([1 2; 3 4])) == [mp0 mp0; mp0 mp0]

# ==============================================================================
# Matrix dense/sparse multiplication
# ==============================================================================

A = [mp0 mp0 3; -2 mp0 mp0; mp0 mp0 0]
B = [mp0 2 mp0; mp0 mp0 0; -3 mp0 mp0]
D = [3 mp0 mp0; mp0 -2 mp0; mp0 mp0 0]
P = [mp0 mp0 0; 0 mp0 mp0; mp0 mp0 0]

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

I = eye(MP,4,4)
Z = spzeros(MP,4,4)
O = ones(MP,4,4)
A = MP(rand(4,4))

@test (A * I) == (I * A) == A
@test (A + Z) == (Z + A) == A
@test (A * Z) == (Z * A) == Z
@test (A + O) == (O + A) == A
@test (A * O) != (O * A) != A

# ==============================================================================
# Max-Plus element by element operations on matrices
# ==============================================================================

# Max-Plus Dense matrix

A = MP([1.0 2.0; 3.0 4.0])
@test 2.0 .+ A == A .+ 2.0 == [MP(2.0) + MP(1.0) MP(2.0) + MP(2.0); MP(2.0) + MP(3.0) MP(2.0) + MP(4.0)] == MP([2.0 2.0; 3.0 4.0])
@test 2.0 .* A == A .* 2.0 == [MP(2.0) * MP(1.0) MP(2.0) * MP(2.0); MP(2.0) * MP(3.0) MP(2.0) * MP(4.0)] == MP([3.0 4.0; 5.0 6.0])

# Max-Plus Sparse matrix

sA = sparse(A)
@test 2.0 .+ sA == sA .+ 2.0 == [MP(2.0) + MP(1.0) MP(2.0) + MP(2.0); MP(2.0) + MP(3.0) MP(2.0) + MP(4.0)] == MP([2.0 2.0; 3.0 4.0])
@test 2.0 .* sA == sA .* 2.0 == [MP(2.0) * MP(1.0) MP(2.0) * MP(2.0); MP(2.0) * MP(3.0) MP(2.0) * MP(4.0)] == MP([3.0 4.0; 5.0 6.0])

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

F = [1 MP(2.0); 3 4]
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

sF = sparse([1 MP(2.0); 3 4])
@test sF isa SparseMatrixCSC{MP, Int64}

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

# Max-Plus sparse array to Max-Plus dense array

Z = dense(spzeros(MP,2,2))
@test typeof(Z) == Matrix{MP}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to Max-Plus dense array

Z = full(spzeros(MP,2,2))
@test typeof(Z) == Matrix{MP}
@test Z == [mp0 mp0; mp0 mp0]

# ==============================================================================
# Bug with Julia with SparseMatrixCSC and operator== which confused zero() and 0.0
# ==============================================================================

A = MP(sparse([1, 2], [1, 2], [0.0, 0.0]))
B = spzeros(MP,2,2)
@test A.nzval == MP([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mp0, mp0])
BB = spzeros(MP,2,2)
@test AA.nzval == MP([-Inf, -Inf])
@test (AA == BB) == true

# ==============================================================================
# Sparse Max-Plus element insertion
# ==============================================================================

# Note: ScicosLab will return isempty == false and length = 2
# Note: NSP, like Julia will return isempty == false and length = 4

O = spzeros(MP,2,2)
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
# Max-Plus scalar power operator
# ==============================================================================

A = MP([4 3; 7 -Inf])
@test A^0 == eye(MP,2,2)
@test A * A == A^2 == MP([10.0 7.0; 11.0 10.0])
@test A * A * A == A^3 == MP([14.0 13.0; 17.0 14.0])
@test A^-1 == MP([-4 -7; -3 mp0])

# ==============================================================================
# Max-Plus trace
# ==============================================================================

A = [5 mp0 5; mp0 6 3; 11 12 11]

@test tr(A) == tr(sparse(A)) == MP(11.0)
@test tr(eye(MP,2,2)) == tr(speye(MP,2,2)) == 0.0
@test tr(zeros(MP,2,2)) == tr(spzeros(MP,2,2)) == mp0
@test tr(MP([1.0 2.0; 3.0 4.0])) == tr(sparse(MP([1.0 2.0; 3.0 4.0]))) == MP(1.0) + MP(4.0) == MP(4.0)

# ==============================================================================
# Max-Plus norm
# ==============================================================================

@test norm(MP([1.0 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399.0)
@test norm(MP([1 20 2; 30 400 4; 4 50 10])) == MP(400 - 1) == MP(399)
@test norm(sparse(MP([1 20 2; 30 400 4; 4 50 10]))) == MP(400 - 1) == MP(399)
@test norm([mp0 1; 10 mp0]) == MP(10.0 - -Inf) == mptop

# ==============================================================================
# Max-Plus inverse
# ==============================================================================

A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
@test inv(A) == A^-1 == [mp0 -2 mp0; -1 mp0 mp0; mp0 mp0 -3]
@test A * inv(A) == inv(A) * A == eye(MP,3,3)
@test typeof(inv(A)) == typeof(A^-1) == Matrix{MP}

A = [mp0 1 mp0; 2 mp0 mp0]
@test inv(A) == A^-1 == [mp0 -2; -1 mp0; mp0 mp0]
@test A * inv(A) == eye(MP,2,2)
@test inv(A) * A == [0 mp0 mp0; mp0 0 mp0; mp0 mp0 mp0]
@test A * inv(A) != inv(A) * A

#FIXME @test_throws ErrorException("The matrix cannot be inversed") inv(MP([1 2; 3 4]))

# ==============================================================================
# TODO: Sparse matrixe x dense vector
# ==============================================================================

# ==============================================================================
# Matrix residuation
# ==============================================================================

A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
B = [3 mp0 mp0; mp0 mp0 4; mp0 5 mp0]
x = A \ B
@test x == MP([mp0 mp0 2; 2 mp0 mp0; mp0 2 mp0])
@test A * x == B
@test (A \ A) == (A / A) == eye(MP,3,3)
@test (B \ B) == (B / B) == eye(MP,3,3)

###

A = MP([3.0 4; 5 6])
B = MP([0.0 -2; 2 0])
x = A \ B
@test x == MP([-3 -5; -4 -6])
@test A * x == B

@test (A \ A) == MP([0 1; -1 0])
@test (A / A) == (B \ B) == (B / B) == MP([0 -2; 2 0])

# ==============================================================================
# Max-Plus star
# ==============================================================================

# Scalars

@test star(MP(2)) == mptop
@test star(MP(1.0)) == mptop
@test star(MP(-1.0)) == mp1
@test star(mp0) == mp1
@test star(mp1) == mp1
@test star(mptop) == mptop

# Matrices

@test_throws ErrorException("Matrix shall be squared") star(MP([]))
@test_throws ErrorException("Matrix shall be squared") star(eye(MP,3,2))
@test star(eye(MP,2,2)) == eye(MP,2,2)
@test star(full(spzeros(MP,2,2))) == eye(MP,2,2)

#

@test star(MP([1 2; 3 4])) == [mptop mptop; mptop mptop]
A = MP([-3 -2; -1 0]); B = star(A)
@test B == eye(MP,2,2) + A
@test B == B + A^2
@test B == B + A^3

#

@test star(MP([-1 2; mp0 -3])) == MP([0 2; -Inf 0])
A = [mp0 2 3; -2 -10 -1; -5 -2 mp1]
@test star(A) == MP([0 2 3; -2 0 1; -4 -2 0])
@test star(star(A)) == star(A)
@test star(A) == (A^0 + A)^2
@test (A^0 + A)^2 == (A^0 + A)^3

# FIXME KO: Mixing +inf and -inf

@test_broken star(MP([2 3; mp0 -1])) == [mi0 mi0; mp0 mp1]

# Random large matrix

A = MP(rand(64,64))
@test star(A) == fill(mptop, 64,64)

# FIXME KO
B = (((ones(1, size(A,1)) * A * ones(size(A,2), 1))[1,1])^-1) * A
@test_broken maximum(plustimes(B)) == 0.0

# ==============================================================================
# Max-Plus A* and A+
# ==============================================================================

# Scalars

@test star(MP(2)) == mptop
@test star(MP(1.0)) == mptop
@test star(MP(-1.0)) == mp1
@test star(mp0) == mp1
@test star(mp1) == mp1
@test star(mptop) == mptop

@test plus(MP(2)) == mptop
@test plus(MP(1.0)) == mptop
@test plus(MP(-1.0)) == MP(-1)
@test plus(mp0) == mp0
@test plus(mp1) == mp1
@test plus(mptop) == mptop

# Matrices

@test_throws ErrorException("Matrix shall be squared") star(MP([]))
@test_throws ErrorException("Matrix shall be squared") star(MP([1; 2]))
@test_throws ErrorException("Matrix shall be squared") star(eye(MP,3,2))

A = [mp0 2 3; -2 -10 -1; -5 -2 mp1]
B = star(A)
@test B == MP([0 2 3; -2 0 1; -4 -2 0])
@test B == A * star(A) == plus(A)

A = [mp0 2 3; -10 -10 -1; -5 -2 mp1]
B = star(A)
@test B == MP([0 2 3; -6 0 -1; -5 -2 0])
@test B == A * star(A) == plus(A)

# What happens if a circuit has strictly positive weight ?
A = [mp0 2 3; -10 -10 -1; 6 -2 mp1]
B = star(A)
@test B == [mptop mptop mptop; mptop mptop mptop; mptop mptop mptop]
@test B == A * star(A) == plus(A)

# ==============================================================================
# Max-Plus a star b
# ==============================================================================

A = MP([-3 -2; -1 0])
b = MP([mp0; mp1])
x = astarb(A, b)
@test x == MP([-2; 0])
@test x == A * x + b

# ==============================================================================
# Display
# ==============================================================================

using Suppressor

set_tropical_display(0)
result = @capture_out show(stdout, mp0)
@test result == "-Inf"
result = @capture_out show(stdout, mp1)
@test result == "0.0"
result = @capture_out show(stdout, MP(4.0))
@test result == "4.0"
result = @capture_out show(stdout, zero(MP))
@test result == "-Inf"
result = @capture_out show(stdout, one(MP))
@test result == "0.0"

set_tropical_display(1)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MP))
@test result == "."
result = @capture_out show(stdout, one(MP))
@test result == "0"

set_tropical_display(2)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MP))
@test result == "."
result = @capture_out show(stdout, one(MP))
@test result == "e"

set_tropical_display(3)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MP(4))
@test result == "4"
result = @capture_out show(stdout, zero(MP))
@test result == "ε"
result = @capture_out show(stdout, one(MP))
@test result == "0"

set_tropical_display(4)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MP(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MP))
@test result == "ε"
result = @capture_out show(stdout, one(MP))
@test result == "e"

A = MP([4.5 0.0; 7.0 -Inf])
result = @capture_out show(stdout, A)
@test result == "MP[4.5 e; 7 ε]"
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & e \\\\\n7 & \\varepsilon \\\\\n\\end{array}\n\\right]\n"
set_tropical_display(0)
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & 0 \\\\\n7 & -\\infty \\\\\n\\end{array}\n\\right]\n"

tropshow(stdout, [MP(1); 2])
