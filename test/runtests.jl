using MaxPlus
using Test
using SparseArrays

# ==============================================================================
# Check Typeof and promotions
# ==============================================================================

a = MP(1.0)
@test typeof(a) == MP{Float64}

b = MP(1)
@test typeof(b) == MP{Int64}

c = MP{Float64}(1)
@test typeof(c) == MP{Float64}

A = [1 2; 3 4]
@test typeof(A) == Array{Int64,2}
@test typeof(mparray(A)) == Array{MP{Int64},2}

B = [MP(1) MP(2); MP(3) MP(4)]
@test typeof(B) == Array{MP{Int64},2}

C = [MP(1) 2; 3 4]
@test typeof(C) == Array{MP{Int64},2}
@test typeof(array(C)) == Array{Int64,2}

D = [MP(1.0) 2; 3 4]
@test typeof(D) == Array{MP{Float64},2}
@test typeof(array(D)) == Array{Float64,2}

E = [MP(1) 2.0; 3 4]
@test typeof(E) == Array{MP{Int64},2}
@test typeof(array(E)) == Array{Int64,2}

F = [2 MP(1.0); 3 4]
@test typeof(F) == Array{MP{Float64},2}
@test typeof(array(F)) == Array{Float64,2}

@test typeof(array(mparray(A))) == typeof(A)
@test typeof(mparray(array(B))) == typeof(B)

# ==============================================================================
# Max-plus comparaisons
# ==============================================================================

b = MP(3.0)
@test (b == b) == true
@test (b != b) == (b ≠ b) == false
@test (b >= b) == (b ≥ b) == true
@test (b <= b) == (b ≤ b) == true
@test (b > b) == false
@test (b < b) == false

# TODO: matrix of bool
# B = [MP(1) MP(2); MP(3) MP(4)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
# @test (B >= B) == (B ≥ B) == true
# @test (B <= B) == (B ≤ B) == true
# @test (B > B) == false
# @test (B < B) == false

@test array(mparray(A)) == A
@test mparray(array(B)) == B

# ==============================================================================
# Max-Plus zero
# ==============================================================================

@test zero(MP{Float64}) == MP(-Inf) == -Inf
@test mp0 == zero(MP{Float64})
@test zero(MP{Int64}) == -9223372036854775808
@test zero(MP{Bool}) == false

# ==============================================================================
# Max-Plus one
# ==============================================================================

@test one(MP{Float64}) == MP(0.0) == 0.0
@test mp1 == one(MP{Float64})
@test one(MP{Int64}) == 0
@test one(MP{Bool}) == false

# ==============================================================================
# Max-plus operations and neutral elements and commutativity
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

@test MP(2)^3 == MP(2 * 3) == MP(6)

# ==============================================================================
# Max-plus distributivity of operations
# ==============================================================================

c = MP(1.0)
@test (a + b) + c == a + (b + c) == a + b + c == MP(max(1.0, 3.0, 5.0)) == MP(5.0)
@test (a * b) * c == a * (b * c) == a * b * c == MP(1.0 + 3.0 + 5.0) == MP(9.0)

# ==============================================================================
# Max-plus distributivity
# ==============================================================================

@test (a + b) * c == MP(max(5.0, 3.0) + 1.0) == MP(6.0)
@test (a * c) + (b * c) == MP(max(5.0 + 1.0, 3.0 + 1.0)) == MP(6.0)

# ==============================================================================
# Max-plus other operations
# ==============================================================================

@test b / b == MP(3.0 - 3.0) == MP(0.0)
@test b - b == MP(3.0 - 3.0) == MP(0.0)

@test MP(3.0) / MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)
@test MP(3.0) - MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)

# TODO min(a,b,c)

# ==============================================================================
# Max-plus between objects of different types
# ==============================================================================

b = MP(3.0)
@test b + 3.0 == 3.0 + b == b == MP(max(3.0, 3.0)) == MP(3.0)
@test b + 3   == 3   + b == b == MP(max(3.0, 3.0)) == MP(3.0)
@test b * 3.0 == 3.0 * b == MP(3.0 + 3.0) == MP(6.0)
@test b * 3   == 3   * b == MP(3.0 + 3.0) == MP(6.0)

# ==============================================================================
# Max-plus dense matrix creation
# ==============================================================================

b = MP(3.0); C = [b 4.0; 5.0 6.0]
@test typeof(C) == Array{MP{Float64},2}

D = mparray([1.0 2.0; 3.0 4.0])
@test typeof(D) == Array{MP{Float64},2}

# Max-Plus matrix of max-plus ones
O = mpones(Float64, 2)
@test typeof(O) == Array{MP{Float64},1}
@test O == [mp1; mp1]

O = mpones(Float64, 2,5)
@test typeof(O) == Array{MP{Float64},2}
@test O == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]

# Identity matrix
Id = mpeye(Float64, 2,5)
@test typeof(Id) == Array{MP{Float64},2}
@test Id == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]

# TODO
# E = mparray([1, 2, 3], [1, 2, 3], [0, 2, 0])

# ==============================================================================
# Max-plus sparse matrix creation and conversion
# ==============================================================================

# Basic dense non max-plus matrix to max-plus sparse array
A = mpsparse([1.0 2.0; 3.0 4.0])
@test typeof(A) == SparseMatrixCSC{MP{Float64},Int64}
@test (A.m == 2) && (A.n == 2)
@test nonzeros(A) == mparray([1.0; 3.0; 2.0; 4.0])

A = mpsparse(mparray([1.0 2.0; 3.0 4.0]))
@test typeof(A) == SparseMatrixCSC{MP{Float64},Int64}
@test (A.m == 2) && (A.n == 2)
@test nonzeros(A) == mparray([1.0; 3.0; 2.0; 4.0])

# 0.0 are not eliminated in the max-plus sparse array
B = mpsparse([1.0 2.0; 0.0 4.0])
@test typeof(B) == SparseMatrixCSC{MP{Float64},Int64}
@test (B.m == 2) && (B.n == 2)
@test nonzeros(B) == mparray([1.0; 0.0; 2.0; 4.0])

B = mpsparse(mparray([1.0 2.0; 0.0 4.0]))
@test typeof(B) == SparseMatrixCSC{MP{Float64},Int64}
@test (B.m == 2) && (B.n == 2)
@test nonzeros(B) == mparray([1.0; 0.0; 2.0; 4.0])

# -Inf are eliminated in the max-plus sparse array
C = mpsparse([1.0 2.0; -Inf 4.0])
@test typeof(C) == SparseMatrixCSC{MP{Float64},Int64}
@test (C.m == 2) && (C.n == 2)
@test nonzeros(C) == mparray([1.0; 2.0; 4.0])

C = mpsparse(mparray([1.0 2.0; -Inf 4.0]))
@test typeof(C) == SparseMatrixCSC{MP{Float64},Int64}
@test (C.m == 2) && (C.n == 2)
@test nonzeros(C) == mparray([1.0; 2.0; 4.0])

# Max-Plus matrix of zeros is already a max-plus sparse array
D = mpzeros(Float64, 3,4)
@test typeof(D) == SparseMatrixCSC{MP{Float64},Int64}
@test (D.m == 3) && (D.n == 4)
@test nonzeros(D) == []

# Convert max-plus array to max-plus sparse array
E = mpsparse(mpones(Float64, 3,4))
@test typeof(E) == SparseMatrixCSC{MP{Float64},Int64}
@test (E.m == 3) && (E.n == 4)
@test nonzeros(E) == mparray(zeros(Float64, 12))

# Convert non max-plus array to max-plus sparse array
F = mpsparse(ones(Float64, 3,4))
@test typeof(F) == SparseMatrixCSC{MP{Float64},Int64}
@test (F.m == 3) && (F.n == 4)
@test nonzeros(F) == ones(12) .+ MP(1.0)

#
G = mpsparse(mparray([1.0 2.0 3.0; 1.0 2.0 3.0; 0.0 2.0 0.0]))
@test typeof(G) == SparseMatrixCSC{MP{Float64},Int64}
@test (G.m == 3) && (G.n == 3)
@test nonzeros(G) == mparray([1.0; 1.0; 0.0; 2.0; 2.0; 2.0; 3.0; 3.0; 0.0])

# max-plus sparse array to max-plus dense array
Z = dense(mpzeros(Float64, 2,2))
@test typeof(Z) == Array{MP{Float64},2}
@test Z == [mp0 mp0; mp0 mp0]

# max-plus sparse array to max-plus dense array
Z = full(mpzeros(Float64, 2,2))
@test typeof(Z) == Array{MP{Float64},2}
@test Z == [mp0 mp0; mp0 mp0]

# max-plus sparse array to non max-plus dense array
Z = array(mpzeros(Float64, 2,2))
@test typeof(Z) == SparseMatrixCSC{Float64,Int64}
@test Z == [-Inf -Inf; -Inf -Inf]

# ==============================================================================
# Max-plus operations with matrices
# ==============================================================================

A = mparray([4 3; 7 -Inf])
@test A * A == A^2 == mparray([10.0 7.0; 11.0 10.0])

C = [MP(3.0) 4.0; 5.0 6.0]
D = mparray([1.0 2; 3 4])
@test C + D == D + C == mparray([3.0 4.0; 5.0 6.0]) == C

# Trace
@test mptrace(mpeye(Float64, 2,2)) == 0.0
@test mptrace(mpeye(Float64, 2,5)) == 0.0
@test mptrace(mpzeros(Float64, 2,2)) == mptrace(full(mpzeros(Float64, 2,2))) == mp0
@test mptrace(mpzeros(Float64, 2,5)) == mptrace(full(mpzeros(Float64, 2,5))) == mp0
@test mptrace([1.0 2.0; 3.0 4.0]) == MP(1.0) + MP(4.0) == MP(4.0)

# FIXME @test C / C == mparray([0.0 -2.0; 2.0 0.0])

# ==============================================================================
#
# ==============================================================================
