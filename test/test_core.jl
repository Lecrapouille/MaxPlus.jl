# ==============================================================================
# Check type of Max-Plus and type promotions on scalars
# ==============================================================================

# ==============================================================================
# Scalar constructor
# ==============================================================================

# Min-Plus

a = MinPlus(1.0)
@test typeof(a) == MinPlus
@test a isa MinPlus
@test a.v isa Float64
@test a.v == 1.0
@test isnan(a) == false

b = MinPlus(1)
@test typeof(b) == MinPlus
@test b isa MinPlus
@test b.v isa Float64
@test b.v == 1
@test isnan(b) == false

c = MinPlus(-1)
@test c isa MinPlus
@test c.v isa Float64
@test c.v == -1.0
@test isnan(c) == false

d = -MinPlus(1)
@test d isa MinPlus
@test d.v isa Float64
@test d.v == -1.0
@test isnan(d) == false

e = MinPlus(NaN)
@test e isa MinPlus
@test e.v isa Float64
@test isnan(e) == true

f = Trop(MinPlus, NaN)
@test f isa MinPlus
@test f.v isa Float64
@test isnan(f) == false
@test f.v == Inf

g = Trop(MinPlus, 42.0)
@test g isa MinPlus
@test g.v isa Float64
@test isnan(g) == false
@test g.v == 42.0

# Max-Plus

a = MaxPlus(1.0)
@test typeof(a) == MaxPlus
@test a isa MaxPlus
@test a.v isa Float64
@test a.v == 1.0
@test isnan(a) == false

b = MaxPlus(1)
@test typeof(b) == MaxPlus
@test b isa MaxPlus
@test b.v isa Float64
@test b.v == 1
@test isnan(b) == false

c = MaxPlus(-1)
@test c isa MaxPlus
@test c.v isa Float64
@test c.v == -1.0
@test isnan(c) == false

d = -MaxPlus(1)
@test d isa MaxPlus
@test d.v isa Float64
@test d.v == -1.0
@test isnan(d) == false

e = MaxPlus(NaN)
@test e isa MaxPlus
@test e.v isa Float64
@test isnan(e) == true

f = Trop(MaxPlus, NaN)
@test f isa MaxPlus
@test f.v isa Float64
@test isnan(f) == false
@test f.v == -Inf

g = Trop(MaxPlus, 42.0)
@test g isa MaxPlus
@test g.v isa Float64
@test isnan(g) == false
@test g.v == 42.0

# ==============================================================================
# Scalar copy constructor
# ==============================================================================

a = MaxPlus(MaxPlus(1))
@test a isa MaxPlus
@test a.v == 1.0

b = MinPlus(MinPlus(2))
@test b isa MinPlus
@test b.v == 2.0

c = MaxPlus(MinPlus(3))
@test c isa MaxPlus
@test c.v == 3.0

d = MinPlus(MaxPlus(4))
@test d isa MinPlus
@test d.v == 4.0

# ==============================================================================
# Boolean constructor
# ==============================================================================

@test MaxPlus(true) == mp1
@test MaxPlus(false) == mp0
@test MinPlus(true) == mi1
@test MinPlus(false) == mi0

# ==============================================================================
# Signs
# ==============================================================================

@test sign(MaxPlus(0)) == 0
@test sign(MaxPlus(5)) == 1
@test sign(MaxPlus(-5)) == -1
@test sign(MinPlus(0)) == 0
@test sign(MinPlus(5)) == 1
@test sign(MinPlus(-5)) == -1
@test sign(mp1) == 0
@test sign(mp0) == -1
@test sign(mi1) == 0
@test sign(mi0) == 1

# ==============================================================================
# Scalar comparaisons
# ==============================================================================

# Min-Plus scalars

a = MinPlus(3.0)
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

# Max-Plus scalars

a = MaxPlus(3.0)
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
# Max/Min-Plus scalar to classic algebra scalar conversion
# ==============================================================================

# Min-Plus

@test plustimes(MinPlus(3.0)) isa Float64
@test plustimes(MinPlus(3.0)) == 3.0

@test plustimes(MinPlus(-3.0)) isa Float64
@test plustimes(MinPlus(-3.0)) == -3.0

# Max-Plus

@test plustimes(MaxPlus(2.0)) isa Float64
@test plustimes(MaxPlus(2.0)) == 2.0

@test plustimes(MaxPlus(-2.0)) isa Float64
@test plustimes(MaxPlus(-2.0)) == -2.0

# ==============================================================================
# Neutral / absorbing constant elements
# ==============================================================================

# Min-Plus

@test mi0 isa MinPlus
@test mi0.v isa Float64
@test mi0.v == Inf
@test iszero(mp0) == true
@test isone(mp0) == false

@test mi1 isa MinPlus
@test mi1.v isa Float64
@test mi1.v == 0.0
@test iszero(mi1) == false
@test isone(mi1) == true

@test mitop isa MinPlus
@test mitop.v isa Float64
@test mitop.v == -Inf
@test iszero(mitop) == false
@test isone(mitop) == false

# Max-Plus

@test mp0 isa MaxPlus
@test mp0.v isa Float64
@test mp0.v == -Inf
@test iszero(mp0) == true
@test isone(mp0) == false

@test mp1 isa MaxPlus
@test mp1.v isa Float64
@test mp1.v == 0.0
@test iszero(mp1) == false
@test isone(mp1) == true

@test mptop isa MaxPlus
@test mptop.v isa Float64
@test mptop.v == Inf
@test iszero(mptop) == false
@test isone(mptop) == false

# ==============================================================================
# zero / one elements
# ==============================================================================

# Min-Plus

@test zero(MinPlus) == MinPlus(Inf) == Inf
@test zero(MinPlus) == mi0 == Inf
@test zero(MinPlus(42.0)) == mi0 == Inf
@test zero(MinPlus(42)) == mi0 == Inf

@test one(MinPlus) == MinPlus(0.0) == 0.0
@test one(MinPlus) == mi1 == 0.0
@test one(MinPlus(42.0)) == mi1 == 0.0
@test one(MinPlus) == mi1 == 0.0
@test one(MinPlus(42)) == mi1 == 0

# Max-Plus

@test zero(MaxPlus) == MaxPlus(-Inf) == -Inf
@test zero(MaxPlus) == mp0 == -Inf
@test zero(MaxPlus(42.0)) == mp0 == -Inf
@test zero(MaxPlus(42)) == mp0 == -Inf

@test one(MaxPlus) == MaxPlus(0.0) == 0.0
@test one(MaxPlus) == mp1 == 0.0
@test one(MaxPlus(42.0)) == mp1 == 0.0
@test one(MaxPlus) == mp1 == 0.0
@test one(MaxPlus(42)) == mp1 == 0

# ==============================================================================
# Operations on scalars and neutral elements and commutativity
# ==============================================================================

# Min-Plus

b = MinPlus(3.0)
@test b + mi0 == mi0 + b == MinPlus(min(3.0, Inf)) == MinPlus(min(Inf, 3.0)) == b == MinPlus(3.0)
@test b * mi1 == mi1 * b == MinPlus(3.0 + 0.0) == MinPlus(0.0 + 3.0) == b == MinPlus(3.0)

@test b * b == MinPlus(3.0 + 3.0) == MinPlus(6.0)
@test b + b == MinPlus(min(3.0, 3.0)) == MinPlus(3.0)

a = MinPlus(5.0)
b = MinPlus(3.0)
@test a + b == MinPlus(min(5.0, 3.0)) == b + a == MinPlus(min(3.0, 5.0)) == MinPlus(3.0)
@test a * b == MinPlus(5.0 + 3.0) == b * a == MinPlus(3.0 + 5.0) == MinPlus(8.0)

# Max-Plus

b = MaxPlus(3.0)
@test b + mp0 == mp0 + b == MaxPlus(max(3.0, -Inf)) == MaxPlus(max(-Inf, 3.0)) == b == MaxPlus(3.0)
@test b * mp1 == mp1 * b == MaxPlus(3.0 + 0.0) == MaxPlus(0.0 + 3.0) == b == MaxPlus(3.0)

@test b * b == MaxPlus(3.0 + 3.0) == MaxPlus(6.0)
@test b + b == MaxPlus(max(3.0, 3.0)) == MaxPlus(3.0)

a = MaxPlus(5.0)
b = MaxPlus(3.0)
@test a + b == MaxPlus(max(5.0, 3.0)) == b + a == MaxPlus(max(3.0, 5.0)) == MaxPlus(5.0)
@test a * b == MaxPlus(5.0 + 3.0) == b * a == MaxPlus(3.0 + 5.0) == MaxPlus(8.0)

# ==============================================================================
# Operations on scalars of different types
# ==============================================================================

# Min-Plus

b = MinPlus(3.0)
@test b + 3.0 == 3.0 + b == b == MinPlus(min(3.0, 3.0)) == MinPlus(3.0)
@test b + 3  == 3 + b == b == MinPlus(min(3.0, 3.0)) == MinPlus(3.0)
@test b * 3.0 == 3.0 * b == MinPlus(3.0 + 3.0) == MinPlus(6.0)
@test b * 3  == 3 * b == MinPlus(3.0 + 3.0) == MinPlus(6.0)

# Max-Plus

b = MaxPlus(3.0)
@test b + 3.0 == 3.0 + b == b == MaxPlus(max(3.0, 3.0)) == MaxPlus(3.0)
@test b + 3  == 3 + b == b == MaxPlus(max(3.0, 3.0)) == MaxPlus(3.0)
@test b * 3.0 == 3.0 * b == MaxPlus(3.0 + 3.0) == MaxPlus(6.0)
@test b * 3  == 3 * b == MaxPlus(3.0 + 3.0) == MaxPlus(6.0)

# ==============================================================================
#  Distributivity of operations
# ==============================================================================

# Min-Plus

a = MinPlus(5.0)
b = MinPlus(3.0)
c = MinPlus(1.0)
@test (a + b) + c == a + (b + c) == a + b + c == MinPlus(min(1.0, 3.0, 5.0)) == MinPlus(1.0)
@test (a * b) * c == a * (b * c) == a * b * c == MinPlus(1.0 + 3.0 + 5.0) == MinPlus(9.0)
@test (a + b) * c == MinPlus(min(5.0, 3.0) + 1.0) == MinPlus(4.0)
@test (a * c) + (b * c) == MinPlus(min(5.0 + 1.0, 3.0 + 1.0)) == MinPlus(4.0)

# Max-Plus

a = MaxPlus(5.0)
b = MaxPlus(3.0)
c = MaxPlus(1.0)
@test (a + b) + c == a + (b + c) == a + b + c == MaxPlus(max(1.0, 3.0, 5.0)) == MaxPlus(5.0)
@test (a * b) * c == a * (b * c) == a * b * c == MaxPlus(1.0 + 3.0 + 5.0) == MaxPlus(9.0)
@test (a + b) * c == MaxPlus(max(5.0, 3.0) + 1.0) == MaxPlus(6.0)
@test (a * c) + (b * c) == MaxPlus(max(5.0 + 1.0, 3.0 + 1.0)) == MaxPlus(6.0)

# ==============================================================================
# Residuation
# ==============================================================================

# Min-Plus TODO

# Max-Plus

@test MaxPlus(3) / MaxPlus(3) == -(MaxPlus(3) * MaxPlus(-3)) == MaxPlus(0.0)
@test MaxPlus(3) / MaxPlus(5) == -(MaxPlus(5.0) * MaxPlus(-3.0)) == MaxPlus(-2.0)
@test MaxPlus(3) \ MaxPlus(3) == -(MaxPlus(3) * MaxPlus(-3)) == MaxPlus(0.0)
@test MaxPlus(3) \ MaxPlus(6) == -(MaxPlus(3.0) * MaxPlus(-6.0)) == MaxPlus(3.0)
@test MaxPlus(3) \ mp0 == mp0
@test MaxPlus(3) \ mp1 == MaxPlus(-3)

# ==============================================================================
# Forbiden minus operator
# ==============================================================================

@test_throws ErrorException("Minus operator does not exist in (min,+) algebra") MinPlus(3) - MinPlus(3)
@test_throws ErrorException("Minus operator does not exist in (min,+) algebra") MinPlus(3) - 3
@test_throws ErrorException("Minus operator does not exist in (min,+) algebra") 3 - MinPlus(3)
@test_throws ErrorException("Minus operator does not exist in (max,+) algebra") MaxPlus(3) - MaxPlus(3)
@test_throws ErrorException("Minus operator does not exist in (max,+) algebra") MaxPlus(3) - 3
@test_throws ErrorException("Minus operator does not exist in (max,+) algebra") 3 - MaxPlus(3)

# ==============================================================================
# Scalar power operator
# ==============================================================================

# Min-Plus TODO

# Max-Plus

@test MaxPlus(2)^4  == MaxPlus(2 * 4) == MaxPlus(8)
@test MaxPlus(0)^0  == MaxPlus(0 * 0) == MaxPlus(0)
@test MaxPlus(2)^0  == MaxPlus(2 * 0) == MaxPlus(0)
@test MaxPlus(2)^-3 == MaxPlus(2)^(-3) == MaxPlus(2 * -3) == MaxPlus(-6)
@test MaxPlus(2)^0.5 == MaxPlus(2 * 0.5) == MaxPlus(1.0)
@test MaxPlus(2)^(-0.5) == MaxPlus(2 * -0.5) == MaxPlus(-1.0)
@test mp0^0 == MaxPlus(0.0)
@test mp0^2 == mp0
@test_broken mp0^(-2) == mp0
@test inv(MaxPlus(5)) == MaxPlus(5)^-1 == MaxPlus(-5)

# ==============================================================================
# Max-Plus min operator
# ==============================================================================

@test min(MaxPlus(3.0), mp0) == min(mp0, MaxPlus(3.0)) == mp0
@test min(MaxPlus(3.0), mp1) == min(mp1, MaxPlus(3.0)) == mp1
@test min(MaxPlus(1), MaxPlus(2)) == min(1, MaxPlus(2)) == min(MaxPlus(1), 2) == MaxPlus(1)
#@test min(MaxPlus([10 1; 10 1]), MaxPlus([4 5; 6 5])) == MaxPlus([4 1; 6 1])
#@test min(sparse(MaxPlus([10.0 1; 10 1])), sparse(MaxPlus([4.0 5; 6 5]))) == sparse(MaxPlus([4.0 1; 6 1]))
#@test min(sparse(MaxPlus([10 1; 10 1])), sparse(MaxPlus([4 5; 6 5]))) == sparse(MaxPlus([4 1; 6 1]))
#@test min(sparse(eye(MaxPlus,2,2)), spzeros(MaxPlus,2,2)) == spzeros(MaxPlus,2,2)
# @test min(eye(MaxPlus,2,2), ones(2,2)) == eye(MaxPlus,2,2)

# ==============================================================================
# Border cases operations
# ==============================================================================

# Min-Plus

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
#@test mi0 * mitop == mitop
@test mi1 * mi0 == mi0
@test mi1 * mi1 == mi1
@test mi1 * mitop == mitop
#@test mitop * mi0 == mitop
@test mitop * mi1 == mitop
@test mitop * mitop == mitop

#@test mi0 / mi0 == mi0
#@test mi0 / mi1 == mi0
#@test mi0 / mitop == mi0
#@test mi1 / mi0 == mi1
#@test mi1 / mi1 == mi1
#@test mi1 / mitop == mi1
#@test mitop / mi0 == mi1
#@test mitop / mi1 == mi1
#@test mitop / mitop == mi1

#@test mi0 \ mi0 == mi0
#@test mi0 \ mi1 == mi0
#@test mi0 \ mitop == mi0
#@test mi1 \ mi0 == mi1
#@test mi1 \ mi1 == mi1
#@test mi1 \ mitop == mi1
#@test mitop \ mi0 == mi1
#@test mitop \ mi1 == mi1
#@test mitop \ mitop == mi1

# Max-Plus

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

# Min-Plus

@test abs2(MinPlus(3.0)) == MinPlus(6.0)
@test abs(MinPlus(-3.0)) == MinPlus(3.0)
@test abs(mi0) == mi0
@test float(MinPlus(2)) == MinPlus(2)
#@test log(MinPlus(2)) == MinPlus(log(2.0)) == log(2.0)
#@test round(MinPlus(1.7)) == MinPlus(round(1.7)) == round(1.7)
#@test floor(MinPlus(1.7)) == MinPlus(floor(1.7)) == floor(1.7)

# Max-Plus

@test abs2(MaxPlus(3.0)) == MaxPlus(6.0)
@test abs(MaxPlus(-3.0)) == MaxPlus(3.0)
@test abs(mp0) == mptop
@test float(MaxPlus(2)) == MaxPlus(2)
#@test log(MaxPlus(2)) == MaxPlus(log(2.0)) == log(2.0)
#@test round(MaxPlus(1.7)) == MaxPlus(round(1.7)) == round(1.7)
#@test floor(MaxPlus(1.7)) == MaxPlus(floor(1.7)) == floor(1.7)

# ==============================================================================
# Max-Plus matrix comparaisons
# ==============================================================================

# Dense Max-Plus matrix comparaison

B = [MaxPlus(1) MaxPlus(2); MaxPlus(3) MaxPlus(4)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# TODO Dense vector comparaison

# Sparse Max-Plus matrix comparaison

S = sparse([MaxPlus(1) MaxPlus(2); MaxPlus(3) MaxPlus(4)])
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# TODO Sparse vector comparaison

# Sparse/Dense Max-Plus matrix comparaison

@test (S == B) == (B == S) == true

# Sort is using isless()
v = MaxPlus([3, 1, 2]);
@test sort(v) == MaxPlus([1, 2, 3])

# ==============================================================================
# Min-Plus matrix comparaisons
# ==============================================================================

# Dense Min-Plus matrix comparaison

B = [MinPlus(1) MinPlus(2); MinPlus(3) MinPlus(4)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# Sparse Min-Plus matrix comparaison

S = sparse([MinPlus(1) MinPlus(2); MinPlus(3) MinPlus(4)])
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# Sparse/Dense Min-Plus matrix comparaison

@test (S == B) == (B == S) == true

# Sort is using isless()
v = MinPlus([3, 1, 2]);
@test sort(v) == MinPlus([1, 2, 3])

# ==============================================================================
# Range construction
# ==============================================================================

# Min-Plus

A = MinPlus(1:3)
@test A == [MinPlus(1); MinPlus(2); MinPlus(3)]

B = MinPlus(1.0:0.5:3.0)
@test B == [MinPlus(1.0); MinPlus(1.5); MinPlus(2.0); MinPlus(2.5); MinPlus(3.0)]

# Max-Plus

A = MaxPlus(1:3)
@test A == [MaxPlus(1); MaxPlus(2); MaxPlus(3)]

B = MaxPlus(1.0:0.5:3.0)
@test B == [MaxPlus(1.0); MaxPlus(1.5); MaxPlus(2.0); MaxPlus(2.5); MaxPlus(3.0)]

# ==============================================================================
# Dense matrix constructor and type promotion/contamination
# ==============================================================================

# Min-Plus

A = [MinPlus(1) MinPlus(2); MinPlus(3) MinPlus(4)]
@test A isa Matrix{MinPlus}

B = MinPlus([1 2; 3 4])
@test B isa Matrix{MinPlus}

C = [MinPlus(1) 2; 3 4]
@test C isa Matrix{MinPlus}

D = [MinPlus(1.0) 2; 3 4]
@test D isa Matrix{MinPlus}

E = [MinPlus(1) 2.0; 3 4]
@test E isa Matrix{MinPlus}

F = [1 MinPlus(2.0); 3 4]
@test F isa Matrix{MinPlus}

@test A == B == C == D == E == F

# Max-Plus

A = [MaxPlus(1) MaxPlus(2); MaxPlus(3) MaxPlus(4)]
@test A isa Matrix{MaxPlus}

B = MaxPlus([1 2; 3 4])
@test B isa Matrix{MaxPlus}

C = [MaxPlus(1) 2; 3 4]
@test C isa Matrix{MaxPlus}

D = [MaxPlus(1.0) 2; 3 4]
@test D isa Matrix{MaxPlus}

E = [MaxPlus(1) 2.0; 3 4]
@test E isa Matrix{MaxPlus}

F = [1 MaxPlus(2.0); 3 4]
@test F isa Matrix{MaxPlus}

@test A == B == C == D == E == F

# ==============================================================================
# Max-Plus sparse constructor
# ==============================================================================

# Using SparseArray.sparse matrix

sA = sparse([MinPlus(1) MinPlus(2); MinPlus(3) MinPlus(4)])
@test sA isa SparseMatrixCSC{MinPlus, Int64}

sB = MinPlus(sparse([1 2; 3 4]))
@test sB isa SparseMatrixCSC{MinPlus, Int64}

sC = sparse([MinPlus(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MinPlus, Int64}

sD = sparse([MinPlus(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MinPlus, Int64}

sE = sparse([MinPlus(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MinPlus, Int64}

sF = sparse([1 MinPlus(2.0); 3 4])
@test sF isa SparseMatrixCSC{MinPlus, Int64}

@test sA == sB == sC == sD == sE == sF

# Using SparseArray.sparse (I, J, D) vectors

spA = MaxPlus(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test size(spA.nzval,1) == 3
@test spA.nzval == [mp0; 2.0; 0.0]

spB = sparse(MaxPlus([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test size(spB.nzval,1) == 3
@test spB.nzval == [mp0; 2.0; 0.0]
@test spB == spA

spC = MaxPlus(spA)
@test size(spC.nzval,1) == 3
@test spC.nzval == [mp0; 2.0; 0.0]
@test spC == spA

# Using Max-Plus

spA = sparse(MaxPlus([-Inf 0; 0 -Inf]))
@test findnz(spA) == ([2, 1], [1, 2], MaxPlus[0.0, 0.0])
spB = MaxPlus(sparse([-Inf 0; 0 -Inf]))
@test findnz(spB) == ([1, 2], [1, 2], MaxPlus[mp0, mp0])
@test spA != spB
spC = sparse(MaxPlus([4 0; 7 -Inf]))
@test findnz(spC) == ([1, 2, 1], [1, 1, 2], MaxPlus[4, 7, 0])

# ==============================================================================
# Matrix ones, eye, zeros constructions
# ==============================================================================

# Max-Plus Matrix of ones

@test ones(MaxPlus, 2) isa Vector{MaxPlus}
@test ones(MaxPlus, 2) == [mp1; mp1]
@test ones(MaxPlus, 2,5) isa Matrix{MaxPlus}
@test ones(MaxPlus, 2,5) == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]

# Max-Plus Identity matrix

@test eye(MaxPlus, 2) isa Matrix{MaxPlus}
@test eye(MaxPlus, 2) == [mp1 mp0; mp0 mp1]
@test eye(MaxPlus, 2,5) isa Matrix{MaxPlus}
@test eye(MaxPlus, 2,5) == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]

# Max-Plus Matrix of zeros

@test spzeros(MaxPlus, 2) isa SparseVector{MaxPlus, Int64}
@test spzeros(MaxPlus, 2).nzval == MaxPlus([])
@test spzeros(MaxPlus, 2,3) isa SparseMatrixCSC{MaxPlus, Int64}
@test spzeros(MaxPlus, 2,3).nzval == MaxPlus([])

# Min-Plus Matrix of ones

@test ones(MinPlus, 2) isa Vector{MinPlus}
@test ones(MinPlus, 2) == [mi1; mi1]
@test ones(MinPlus, 2,5) isa Matrix{MinPlus}
@test ones(MinPlus, 2,5) == [mi1 mi1 mi1 mi1 mi1; mi1 mi1 mi1 mi1 mi1]

# Min-Plus Identity matrix

@test eye(MinPlus, 2) isa Matrix{MinPlus}
@test eye(MinPlus, 2) == [mi1 mi0; mi0 mi1]
@test eye(MinPlus, 2,5) isa Matrix{MinPlus}
@test eye(MinPlus, 2,5) == [mi1 mi0 mi0 mi0 mi0; mi0 mi1 mi0 mi0 mi0]

# Min-Plus Matrix of zeros

@test spzeros(MinPlus, 2) isa SparseVector{MinPlus, Int64}
@test spzeros(MinPlus, 2).nzval == MinPlus([])
@test spzeros(MinPlus, 2,3) isa SparseMatrixCSC{MinPlus, Int64}
@test spzeros(MinPlus, 2,3).nzval == MinPlus([])

# ==============================================================================
# Max-Plus matrix dense/sparse multiplication
# ==============================================================================

# Min-Plus TODO

# Max-Plus

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
# Max-Plus matrix addition
# ==============================================================================

# Min-Plus TODO

# Max-Plus

A = MaxPlus([3.0 4; 5 6])
B = MaxPlus([1.0 2; 3 4])
C = MaxPlus([3.0 4; 5 6])
@test (A + B) == (B + A) == C

# ==============================================================================
# Matrix ones, eye, zeros operations
# ==============================================================================

# Max-Plus

I = eye(MaxPlus,4,4)
Z = spzeros(MaxPlus,4,4)
O = ones(MaxPlus,4,4)
A = MaxPlus(rand(4,4))

@test (A * I) == (I * A) == A
@test (A + Z) == (Z + A) == A
@test (A * Z) == (Z * A) == Z
@test (A + O) == (O + A) == A
@test (A * O) != (O * A) != A

# Min-Plus

I = eye(MinPlus,4,4)
Z = spzeros(MinPlus,4,4)
O = ones(MinPlus,4,4)
A = MinPlus(rand(4,4))

@test (A * I) == (I * A) == A
@test (A + Z) == (Z + A) == A
@test (A * Z) == (Z * A) == Z
@test (A + O) == (O + A) == O
@test (A * O) != (O * A) != A

# ==============================================================================
# Max-Plus element by element operations on matrices
# ==============================================================================

# Min-Plus TODO

# Max-Plus Dense matrix

A = MaxPlus([1.0 2.0; 3.0 4.0])
@test 2.0 .+ A == A .+ 2.0 == [MaxPlus(2.0) + MaxPlus(1.0) MaxPlus(2.0) + MaxPlus(2.0); MaxPlus(2.0) + MaxPlus(3.0) MaxPlus(2.0) + MaxPlus(4.0)] == MaxPlus([2.0 2.0; 3.0 4.0])
@test 2.0 .* A == A .* 2.0 == [MaxPlus(2.0) * MaxPlus(1.0) MaxPlus(2.0) * MaxPlus(2.0); MaxPlus(2.0) * MaxPlus(3.0) MaxPlus(2.0) * MaxPlus(4.0)] == MaxPlus([3.0 4.0; 5.0 6.0])

# Max-Plus Sparse matrix

sA = sparse(A)
@test 2.0 .+ sA == sA .+ 2.0 == [MaxPlus(2.0) + MaxPlus(1.0) MaxPlus(2.0) + MaxPlus(2.0); MaxPlus(2.0) + MaxPlus(3.0) MaxPlus(2.0) + MaxPlus(4.0)] == MaxPlus([2.0 2.0; 3.0 4.0])
@test 2.0 .* sA == sA .* 2.0 == [MaxPlus(2.0) * MaxPlus(1.0) MaxPlus(2.0) * MaxPlus(2.0); MaxPlus(2.0) * MaxPlus(3.0) MaxPlus(2.0) * MaxPlus(4.0)] == MaxPlus([3.0 4.0; 5.0 6.0])


# ==============================================================================
# Max-Plus type promotion/contamination on dense/sparse matrix
# ==============================================================================

# Dense matrix

A = [1 2; 3 4]
@test A isa Matrix{Int64}
@test MaxPlus(A) isa Matrix{MaxPlus}

B = [MaxPlus(1) MaxPlus(2); MaxPlus(3) MaxPlus(4)]
@test B isa Matrix{MaxPlus}

C = [MaxPlus(1) 2; 3 4]
@test C isa Matrix{MaxPlus}

D = [MaxPlus(1.0) 2; 3 4]
@test D isa Matrix{MaxPlus}

E = [MaxPlus(1) 2.0; 3 4]
@test E isa Matrix{MaxPlus}

F = [1 MaxPlus(2.0); 3 4]
@test F isa Matrix{MaxPlus}

# Sparse matrix

sA = sparse([1 2; 3 4])
@test sA isa SparseMatrixCSC{Int64, Int64}
@test MaxPlus(sA) isa SparseMatrixCSC{MaxPlus, Int64}

sB = sparse([MaxPlus(1) MaxPlus(2); MaxPlus(3) MaxPlus(4)])
@test sB isa SparseMatrixCSC{MaxPlus, Int64}

sC = sparse([MaxPlus(1) 2; 3 4])
@test sC isa SparseMatrixCSC{MaxPlus, Int64}

sD = sparse([MaxPlus(1.0) 2; 3 4])
@test sD isa SparseMatrixCSC{MaxPlus, Int64}

sE = sparse([MaxPlus(1) 2.0; 3 4])
@test sE isa SparseMatrixCSC{MaxPlus, Int64}

sF = sparse([1 MaxPlus(2.0); 3 4])
@test sF isa SparseMatrixCSC{MaxPlus, Int64}

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

@test full(sC) isa Matrix{MaxPlus}
@test full(sD) isa Matrix{MaxPlus}
@test full(sE) isa Matrix{MaxPlus}
@test full(sF) isa Matrix{MaxPlus}

@test full(sC) == C
@test full(sD) == D
@test full(sE) == E
@test full(sF) == F

# dense() is a anlias for full()

@test dense(sC) isa Matrix{MaxPlus}
@test dense(sD) isa Matrix{MaxPlus}
@test dense(sE) isa Matrix{MaxPlus}
@test dense(sF) isa Matrix{MaxPlus}

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

Z = dense(spzeros(MaxPlus,2,2))
@test typeof(Z) == Matrix{MaxPlus}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to Max-Plus dense array

Z = full(spzeros(MaxPlus,2,2))
@test typeof(Z) == Matrix{MaxPlus}
@test Z == [mp0 mp0; mp0 mp0]

# Max-Plus sparse array to non Max-Plus dense array

Z = array(spzeros(MaxPlus,2,2))
@test typeof(Z) == Matrix{Float64}
@test Z == [-Inf -Inf; -Inf -Inf]

# ==============================================================================
# Bug with Julia with SparseMatrixCSC and operator== which confused zero() and 0.0
# ==============================================================================

# Min-Plus

A = MinPlus(sparse([1, 2], [1, 2], [0.0, 0.0]))
B = spzeros(MinPlus,2,2)
@test A.nzval == MinPlus([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mi0, mi0])
BB = spzeros(MinPlus,2,2)
@test AA.nzval == MinPlus([Inf, Inf])
@test (AA == BB) == true

# Max-Plus

A = MaxPlus(sparse([1, 2], [1, 2], [0.0, 0.0]))
B = spzeros(MaxPlus,2,2)
@test A.nzval == MaxPlus([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mp0, mp0])
BB = spzeros(MaxPlus,2,2)
@test AA.nzval == MaxPlus([-Inf, -Inf])
@test (AA == BB) == true

# ==============================================================================
# Sparse Max-Plus element insertion
# ==============================================================================

# Note: ScicosLab will return isempty == false and length = 2
# Note: NSP, like Julia will return isempty == false and length = 4

O = spzeros(MaxPlus,2,2)
@test length(O) == 4
@test isempty(O) == false

# Zero elements are not inserted
O[1,1] = mp0
@test length(O) == 4
@test isempty(O) == false
@test O.nzval == []
@test isempty(nonzeros(O))

# Insert fake zero element
O[1,1] = MaxPlus(0.0)
@test O[1,1] == MaxPlus(0.0)
@test O.nzval == MaxPlus([0.0])
@test !isempty(nonzeros(O))

# Allow replacing a non-zero() element by a zero() element
O[1,1] = mp0
@test O[1,1] == MaxPlus(-Inf)
@test O.nzval == MaxPlus([-Inf])
@test !isempty(nonzeros(O))

# ==============================================================================
# Max-Plus scalar power operator
# ==============================================================================

A = MaxPlus([4 3; 7 -Inf])
@test A^0 == eye(MaxPlus,2,2)
@test A * A == A^2 == MaxPlus([10.0 7.0; 11.0 10.0])
@test A * A * A == A^3 == MaxPlus([14.0 13.0; 17.0 14.0])
@test_broken A^-1

# ==============================================================================
# Max-Plus trace
# ==============================================================================

A = [5 mp0 5; mp0 6 3; 11 12 11]

@test tr(A) == MaxPlus(11.0)
@test tr(MaxPlus([])) == mp0
@test tr(eye(MaxPlus,2,2)) == tr(eye(MaxPlus,2,2)) == 0.0
@test tr(spzeros(MaxPlus,2,2)) == tr(full(spzeros(MaxPlus,2,2))) == tr(spzeros(MaxPlus,2,2)) == mp0
@test tr(MaxPlus([1.0 2.0; 3.0 4.0])) == MaxPlus(1.0) + MaxPlus(4.0) == MaxPlus(4.0)
@test tr(sparse(MaxPlus([1.0 2.0; 3.0 4.0]))) == MaxPlus(4.0)

# ==============================================================================
# Max-Plus norm
# ==============================================================================

@test norm(MaxPlus([1.0 20 2; 30 400 4; 4 50 10])) == MaxPlus(400 - 1) == MaxPlus(399.0)
@test norm(MaxPlus([1 20 2; 30 400 4; 4 50 10])) == MaxPlus(400 - 1) == MaxPlus(399)
@test norm(sparse(MaxPlus([1 20 2; 30 400 4; 4 50 10]))) == MaxPlus(400 - 1) == MaxPlus(399)
@test norm([mp0 1; 10 mp0]) == MaxPlus(10.0 - -Inf) == mptop

# ==============================================================================
# Max-Plus inverse
# ==============================================================================

A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
@test inv(A) == A^-1 == [mp0 -2 mp0; -1 mp0 mp0; mp0 mp0 -3]
@test A * inv(A) == inv(A) * A == eye(MaxPlus,3,3)
@test typeof(inv(A)) == typeof(A^-1) == Matrix{MaxPlus}

A = [mp0 1 mp0; 2 mp0 mp0]
@test inv(A) == A^-1 == [mp0 -2; -1 mp0; mp0 mp0]
@test A * inv(A) == eye(MaxPlus,2,2)
@test inv(A) * A == [0 mp0 mp0; mp0 0 mp0; mp0 mp0 mp0]
@test A * inv(A) != inv(A) * A

#FIXME @test_throws ErrorException("The matrix cannot be inversed") inv(MaxPlus([1 2; 3 4]))

# ==============================================================================
# Matrix residuation
# ==============================================================================

A = [mp0 1 mp0; 2 mp0 mp0; mp0 mp0 3]
B = [3 mp0 mp0; mp0 mp0 4; mp0 5 mp0]
x = A \ B
@test x == MaxPlus([mp0 mp0 2; 2 mp0 mp0; mp0 2 mp0])
@test A * x == B
@test (A \ A) == (A / A) == eye(MaxPlus,3,3)
@test (B \ B) == (B / B) == eye(MaxPlus,3,3)

###

A = MaxPlus([3.0 4; 5 6])
B = MaxPlus([0.0 -2; 2 0])
x = A \ B
@test x == MaxPlus([-3 -5; -4 -6])
@test A * x == B

@test (A \ A) == MaxPlus([0 1; -1 0])
@test (A / A) == (B \ B) == (B / B) == MaxPlus([0 -2; 2 0])

# ==============================================================================
# Max-Plus star
# ==============================================================================

# Scalars

@test star(MaxPlus(2)) == mptop
@test star(MaxPlus(1.0)) == mptop
@test star(MaxPlus(-1.0)) == mp1
@test star(mp0) == mp1
@test star(mp1) == mp1
@test star(mptop) == mptop

# Matrices

@test_throws ErrorException("Matrix shall be squared") star(MaxPlus([]))
@test_throws ErrorException("Matrix shall be squared") star(eye(MaxPlus,3,2))
@test star(eye(MaxPlus,2,2)) == eye(MaxPlus,2,2)
@test star(full(spzeros(MaxPlus,2,2))) == eye(MaxPlus,2,2)

#

@test star(MaxPlus([1 2; 3 4])) == [mptop mptop; mptop mptop]
A = MaxPlus([-3 -2; -1 0]); B = star(A)
@test B == eye(MaxPlus,2,2) + A
@test B == B + A^2
@test B == B + A^3

#

@test star(MaxPlus([-1 2; mp0 -3])) == MaxPlus([0 2; -Inf 0])
A = [mp0 2 3; -2 -10 -1; -5 -2 mp1]
@test star(A) == MaxPlus([0 2 3; -2 0 1; -4 -2 0])
@test star(star(A)) == star(A)
@test star(A) == (A^0 + A)^2
@test (A^0 + A)^2 == (A^0 + A)^3

# FIXME KO: Mixing +inf and -inf

# @test mpstar(MaxPlus([2 3; mp0 -1])) == [mi0 mi0; mp0 mp1]

# Random large matrix

A = MaxPlus(rand(64,64))
@test star(A) == fill(mptop, 64,64)

# FIXME KO
B = (((ones(1, size(A,1)) * A * ones(size(A,2), 1))[1,1])^-1) * A
@test_broken maximum(plustimes(B)) == 0.0

# ==============================================================================
# Max-Plus plus
# ==============================================================================

# Scalars

@test plus(MaxPlus(2)) == mptop
@test plus(MaxPlus(1.0)) == mptop
@test_broken plus(MaxPlus(-1.0)) == MaxPlus(-1)
@test_broken plus(mp0) == mp0
@test plus(mp1) == mp1
@test plus(mptop) == mptop

# Matrices

@test_throws ErrorException("Matrix shall be squared") plus(MaxPlus([]))
@test_throws ErrorException("Matrix shall be squared") plus(eye(MaxPlus,3,2))
A = [mp0 2 3; -2 -10 -1; -5 -2 mp1]
B = plus(A)
@test B == MaxPlus([0 2 3; -2 0 1; -4 -2 0])
@test B == A * star(A)

# FIXME donne le bon resultat dans REPL
A[2,1] = MaxPlus(-10)
@test_broken plus(A) == MaxPlus([-2 2 3; -6 -3 -1; -5 -2 0])

# What happens if a circuit has strictly positive weight ?
A[3,1] = 6
@test_broken plus(A) == fill(mi0, 2,3) # FIXME

# ==============================================================================
# TODO Max-Plus a star b
# ==============================================================================


# ==============================================================================
# Max-Plus Howard
# ==============================================================================

S = sparse(MaxPlus([1 2; 3 4]))
λ,v = howard(S)
@test (λ,v) == (MaxPlus[4, 4], MaxPlus[2, 4])
@test (S * v) == (λ[1] * v)

S = sparse([mp0 2 mp0; mp1 mp0 mp0; mp0 mp0 2])
λ,v = howard(S)
@test (λ,v) == (MaxPlus[1, 1, 2], MaxPlus[1, 0, 2])
B = S / Matrix(Diagonal(λ))
@test B == [mp0 1 mp0; -1 mp0 mp0; mp0 mp0 0]
@test B * v == v

S = sparse([1 1; mp0 2])
λ,v = howard(S)
@test (λ,v) == (MaxPlus[2, 2], MaxPlus[1, 2])
@test (S * v) == (λ[1] * v)
@test (S * [0; mp0]) == (MaxPlus(1) * [0; mp0])

S = sparse([2 1; mp0 mp1])
λ,v = howard(S)
@test (λ,v) == (MaxPlus[2, 0], MaxPlus[2, 0])
B = S / Matrix(Diagonal(λ))
@test B == [0 1; mp0 0]
@test B * v == v


# ==============================================================================
# Display
# ==============================================================================

using Suppressor

# Min-Plus TODO

# Max-Plus

mp_change_display(0)
result = @capture_out show(stdout, mp0)
@test result == "-Inf"
result = @capture_out show(stdout, mp1)
@test result == "0.0"
result = @capture_out show(stdout, MaxPlus(4.0))
@test result == "4.0"
result = @capture_out show(stdout, zero(MaxPlus))
@test result == "-Inf"
result = @capture_out show(stdout, one(MaxPlus))
@test result == "0.0"

mp_change_display(1)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MaxPlus(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MaxPlus))
@test result == "."
result = @capture_out show(stdout, one(MaxPlus))
@test result == "0"

mp_change_display(2)
result = @capture_out show(stdout, mp0)
@test result == "."
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MaxPlus(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MaxPlus))
@test result == "."
result = @capture_out show(stdout, one(MaxPlus))
@test result == "e"

mp_change_display(3)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "0"
result = @capture_out show(stdout, MaxPlus(4))
@test result == "4"
result = @capture_out show(stdout, zero(MaxPlus))
@test result == "ε"
result = @capture_out show(stdout, one(MaxPlus))
@test result == "0"

mp_change_display(4)
result = @capture_out show(stdout, mp0)
@test result == "ε"
result = @capture_out show(stdout, mp1)
@test result == "e"
result = @capture_out show(stdout, MaxPlus(4.0))
@test result == "4"
result = @capture_out show(stdout, zero(MaxPlus))
@test result == "ε"
result = @capture_out show(stdout, one(MaxPlus))
@test result == "e"

A = MaxPlus([4.5 0.0; 7.0 -Inf])
result = @capture_out show(stdout, A)
@test result == "MaxPlus[4.5 e; 7 ε]"
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & e \\\\\n7 & \\varepsilon \\\\\n\\end{array}\n\\right]\n"
mp_change_display(0)
result = @capture_out LaTeX(stdout, A)
@test result == "\\left[\n\\begin{array}{*{20}c}\n4.5 & 0 \\\\\n7 & -\\infty \\\\\n\\end{array}\n\\right]\n"
