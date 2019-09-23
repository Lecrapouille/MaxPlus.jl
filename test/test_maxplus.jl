# ==============================================================================
# Check type of max-plus and type promotions
# ==============================================================================

a = MP(1.0)
@test typeof(a) == MP{Float64}

b = MP(1)
@test typeof(b) == MP{Int64}

c = MP{Float64}(1)
@test typeof(c) == MP{Float64}

A = [1 2; 3 4]
@test typeof(A) == Array{Int64,2}
@test typeof(MP(A)) == ArrMP{Int64,2}

B = [MP(1) MP(2); MP(3) MP(4)]
@test typeof(B) == ArrMP{Int64,2}

C = [MP(1) 2; 3 4]
@test typeof(C) == ArrMP{Int64,2}
@test typeof(array(C)) == Array{Int64,2}

D = [MP(1.0) 2; 3 4]
@test typeof(D) == ArrMP{Float64,2}
@test typeof(array(D)) == Array{Float64,2}

E = [MP(1) 2.0; 3 4]
@test typeof(E) == ArrMP{Int64,2}
@test typeof(array(E)) == Array{Int64,2}

F = [2 MP(1.0); 3 4]
@test typeof(F) == ArrMP{Float64,2}
@test typeof(array(F)) == Array{Float64,2}

@test typeof(array(MP(A))) == typeof(A)
@test typeof(MP(array(B))) == typeof(B)

# ==============================================================================
# Max-plus constructor
# ==============================================================================

c = MP(1.0)
@test typeof(MP(c)) == MP{Float64}
d = MP(c)
@test typeof(d) == MP{Float64}

# ==============================================================================
# Max-plus comparaisons
# ==============================================================================

# Max-plus scalars
b = MP(3.0)
@test (b == b) == true
@test (b != b) == (b ≠ b) == false
@test (b >= b) == (b ≥ b) == true
@test (b <= b) == (b ≤ b) == true
@test (b > b) == false
@test (b < b) == false

# Max-plus vs. real comparaison
@test (b == 3.0) == true
@test (b == 4.0) == false
@test (b != 4.0) == true
@test (b < 4.0) == true
@test (b <= 4.0) == true
@test (b > 4.0) == false
@test (b >= 4.0) == false
@test (4.0 > b) == true
@test (4.0 >= b) == true

# Matrix comparaison
B = [MP(1) MP(2); MP(3) MP(4)]
@test (B == B) == true
@test (B != B) == (B ≠ B) == false
@test (B .>= B) == (B .≥ B) == [true true; true true]
@test (B .<= B) == (B .≤ B) == [true true; true true]
@test (B .> B) == [false false; false false]
@test (B .< B) == [false false; false false]

# ==============================================================================
# Max-Plus zero
# ==============================================================================

@test zero(MP{Float64}) == MP(-Inf) == -Inf
@test mp0 == zero(MP{Float64})
@test zero(MP{Int64}) == -9223372036854775808
@test zero(MP{Bool}) == false
@test zero(MP(42.0)) == mp0
@test iszero(mp0) == true
@test isone(mp1) == true
@test iszero(mp1) == false
@test isone(mp0) == false

# ==============================================================================
# Max-Plus one
# ==============================================================================

@test one(MP{Float64}) == MP(0.0) == 0.0
@test mp1 == one(MP{Float64})
@test one(MP{Int64}) == 0
@test one(MP{Bool}) == false
@test one(MP(42.0)) == mp1

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
# mp0, mp1 and mptop operations
# ==============================================================================

@test mp0 + mp0 == mp0
@test mp0 + mp1 == 0
@test mp0 + mptop == mptop
@test mp1 + mp0 == 0
@test mp1 + mp1 == 0
@test mp1 + mptop == mptop

@test mp0 * mp0 == mp0
@test mp0 * mp1 == mp0
#FIXME @test mp0 * mptop == mp0
@test mp1 * mp0 == mp0
@test mp1 * mp1 == 0
@test mp1 * mptop == mptop

#FIXME @test mp0 - mp0 == mp0
@test mp0 - mp1 == mp0
#FIXME @test mp1 - mp0 == 0
@test mp1 - mp1 == 0

@test mptop - mp0 == mptop
@test mptop - mp1 == mptop
@test mp0 - mptop == mp0
@test mp1 - mptop == mp0

# ==============================================================================
# Max-plus element by element operations
# ==============================================================================

A=MP([1.0 2.0; 3.0 4.0])
@test 2.0 .+ A == A .+ 2.0 == [MP(2.0) + MP(1.0) MP(2.0) + MP(2.0); MP(2.0) + MP(3.0) MP(2.0) + MP(4.0)] == MP([2.0 2.0; 3.0 4.0])
@test 2.0 .* A == A .* 2.0 == [MP(2.0) * MP(1.0) MP(2.0) * MP(2.0); MP(2.0) * MP(3.0) MP(2.0) * MP(4.0)] == MP([3.0 4.0; 5.0 6.0])

# ==============================================================================
# Max-plus other operations
# ==============================================================================

@test MP(2)^4   == MP(2 * 4) == MP(8)
@test MP(0)^0   == MP(0 * 0) == MP(0)
@test MP(2)^0   == MP(2 * 0) == MP(0)
@test MP(2)^-3  == MP(2)^(-3) == MP(2 * -3) == MP(-6)
@test MP(2)^0.5 == MP(2 * 0.5) == MP(1.0)

@test b / b == MP(3.0 - 3.0) == MP(0.0)
@test b - b == MP(3.0 - 3.0) == MP(0.0)

@test MP(3.0) / MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)
@test MP(3.0) - MP(5.0) == MP(3.0 - 5.0) == MP(-2.0)

@test abs2(MP(3.0)) == MP(6.0)

# ==============================================================================
# Max-plus min
# ==============================================================================

b=MP(3.0)
@test min(b, mp0) == mp0
@test min(b, mp1) == mp1
@test min(MP(1), MP(2)) == MP(1)
@test min(MP([10 1; 10 1]), MP([4 5; 6 5])) == MP([4 1; 6 1])
@test min(mpsparse([10 1; 10 1]), mpsparse([4 5; 6 5])) == mpsparse([4 1; 6 1])

# ==============================================================================
# Julia found bugs. Check if they are really fixed.

A = mpzeros(Float64, 2,2)

# Do not insert zero() element
A[1,1] = mp0
@test A.nzval == []
@test isempty(nonzeros(A))

# Insert fake zero element
A[1,1] = MP(0.0)
@test A[1,1] == MP(0.0)
@test A.nzval == MP([0.0])
@test !isempty(nonzeros(A))

# Allow replacing a non-zero() element by a zero() element
A[1,1] = mp0
@test A[1,1] == MP(-Inf)
@test A.nzval == MP([-Inf])
@test !isempty(nonzeros(A))

# ==============================================================================
# Max-plus length() and isempty()
# ==============================================================================
# Note: ScicosLab will return isempty == false and length = 2
# Note: NSP, like Julia will return isempty == false and length = 4

A = mpzeros(Float64, 2,2)
@test length(A) == 4
@test !isempty(A)

A[1,1] = MP(0.0)
@test length(A) == 4
@test !isempty(A)

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

b = MP(3.0); C = [b 4.0; 5.0 6.0]; D = MP([3.0 4.0; 5.0 6.0])
@test C == D
@test typeof(C) == ArrMP{Float64,2}

# Max-Plus matrix of max-plus ones
O = mpones(Float64, 2)
@test typeof(O) == Array{MP{Float64},1}
@test O == [mp1; mp1]

O = mpones(Float64, 2,5)
@test typeof(O) == ArrMP{Float64,2}
@test O == [mp1 mp1 mp1 mp1 mp1; mp1 mp1 mp1 mp1 mp1]

# Identity matrix
Id = mpeye(Float64, 2)
@test Id == mpeye(Float64, 2)
@test typeof(Id) == ArrMP{Float64,2}
@test Id == [mp1 mp0; mp0 mp1]

Id = mpeye(Float64, 2,5)
@test typeof(Id) == ArrMP{Float64,2}
@test Id == [mp1 mp0 mp0 mp0 mp0; mp0 mp1 mp0 mp0 mp0]

# ==============================================================================
# Max-plus / non max-plus conversion
# ==============================================================================

A=MP([1.0 2.0; 3.0 4.0])
n = 2.0
@test plustimes(MP(n)) == n
@test plustimes(mpzeros(Float64,2,2)) == ones(Float64, 2,2) .* mp0
@test array(MP(A)) == plustimes(MP(A)) == A
@test MP(array(B)) == MP(plustimes(B)) == B

# ==============================================================================
# Convert values usable for min-plus
# ==============================================================================
@test minplus(MP(1.0)) == MP(1.0)
@test minplus(mp0) == mptop
@test minplus(mptop) == mp0
@test minplus(MP([Inf 0.0; 7 -Inf])) == MP([-Inf 0.0; 7 Inf])
#FIXME @test minplus(mpsparse([Inf 0.0; 7 -Inf], preserve=true)) == mpsparse([-Inf 0.0; 7 Inf], preserve=true)
A = MP([0 3 Inf 1; 1 2 2 -Inf; -Inf Inf 1 0])
@test minplus(A) == MP([0 3 -Inf 1; 1 2 2 Inf; Inf -Inf 1 0])
@test minplus(minplus(A)) == A

# ==============================================================================
# Max-plus sparse matrix creation and conversion
# ==============================================================================

# Construct a sparse max-plus vector
V = MP(sparse([1.0, 0.0, 1.0]))
@test typeof(V) == SpvMP{Float64,Int64}
@test V.nzval == MP([1.0; 1.0])

# Construct a sparse max-plus vector of zeros
V = mpzeros(Float64, 2)
@test typeof(V) == SpvMP{Float64,Int64}
@test V.nzval == MP([])

# Construct a sparse max-plus matrix
A = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]))
@test typeof(A) == SpaMP{Float64,Int64}
@test A.nzval == [MP(-Inf), MP(2.0), MP(0.0)]

# Construct a sparse max-plus matrix
A = MP(sparse([1, 2, 3], [1, 2, 3], [-Inf, 2, 0]), preserve=false)
@test typeof(A) == SpaMP{Float64,Int64}
@test A.nzval == [MP(2.0)]

# Basic dense non max-plus matrix to max-plus sparse array
A = mpsparse([1.0 2.0; 3.0 4.0])
@test typeof(A) == SpaMP{Float64,Int64}
@test (A.m == 2) && (A.n == 2)
@test nonzeros(A) == MP([1.0; 3.0; 2.0; 4.0])

A = mpsparse(MP([1.0 2.0; 3.0 4.0]))
@test typeof(A) == SpaMP{Float64,Int64}
@test (A.m == 2) && (A.n == 2)
@test nonzeros(A) == MP([1.0; 3.0; 2.0; 4.0])

# 0.0 are not eliminated in the max-plus sparse array
B = mpsparse([1.0 2.0; 0.0 4.0])
@test typeof(B) == SpaMP{Float64,Int64}
@test (B.m == 2) && (B.n == 2)
@test nonzeros(B) == MP([1.0; 0.0; 2.0; 4.0])

B = mpsparse(MP([1.0 2.0; 0.0 4.0]))
@test typeof(B) == SpaMP{Float64,Int64}
@test (B.m == 2) && (B.n == 2)
@test nonzeros(B) == MP([1.0; 0.0; 2.0; 4.0])

# -Inf are eliminated in the max-plus sparse array
C = mpsparse([1.0 2.0; -Inf 4.0])
@test typeof(C) == SpaMP{Float64,Int64}
@test (C.m == 2) && (C.n == 2)
@test nonzeros(C) == MP([1.0; 2.0; 4.0])

C = mpsparse(MP([1.0 2.0; -Inf 4.0]))
@test typeof(C) == SpaMP{Float64,Int64}
@test (C.m == 2) && (C.n == 2)
@test nonzeros(C) == MP([1.0; 2.0; 4.0])

# Max-Plus matrix of zeros is already a max-plus sparse array
D = mpzeros(Float64, 3,4)
@test typeof(D) == SpaMP{Float64,Int64}
@test (D.m == 3) && (D.n == 4)
@test nonzeros(D) == []

# Convert max-plus array to max-plus sparse array
E = mpsparse(mpones(Float64, 3,4))
@test typeof(E) == SpaMP{Float64,Int64}
@test (E.m == 3) && (E.n == 4)
@test nonzeros(E) == MP(zeros(Float64, 12))

# Convert non max-plus array to max-plus sparse array
F = mpsparse(ones(Float64, 3,4))
@test typeof(F) == SpaMP{Float64,Int64}
@test (F.m == 3) && (F.n == 4)
@test nonzeros(F) == ones(12) .+ MP(1.0)

#
G = mpsparse(MP([1.0 2.0 3.0; 1.0 2.0 3.0; 0.0 2.0 0.0]))
@test typeof(G) == SpaMP{Float64,Int64}
@test (G.m == 3) && (G.n == 3)
@test nonzeros(G) == MP([1.0; 1.0; 0.0; 2.0; 2.0; 2.0; 3.0; 3.0; 0.0])

# max-plus sparse array to max-plus dense array
Z = dense(mpzeros(Float64, 2,2))
@test typeof(Z) == ArrMP{Float64,2}
@test Z == [mp0 mp0; mp0 mp0]

# max-plus sparse array to max-plus dense array
Z = full(mpzeros(Float64, 2,2))
@test typeof(Z) == ArrMP{Float64,2}
@test Z == [mp0 mp0; mp0 mp0]

# max-plus sparse array to non max-plus dense array
Z = array(mpzeros(Float64, 2,2))
@test typeof(Z) == SparseMatrixCSC{Float64,Int64}
@test Z == [-Inf -Inf; -Inf -Inf]

# ==============================================================================
# Bug with Julia with SparseMatrixCSC and operator== which confused zero() and 0.0

A = MP(sparse([1, 2], [1, 2], [0.0, 0.0]))
B = mpzeros(Float64, 2,2)
@test A.nzval == MP([0.0, 0.0])
@test (A == B) == false

AA = sparse([1, 2], [1, 2], [mp0, mp0])
BB = mpzeros(Float64, 2,2)
@test AA.nzval == MP([-Inf, -Inf])
@test (AA == BB) == true

# ==============================================================================
# Max-plus operations with matrices
# ==============================================================================

A = MP([4 3; 7 -Inf])
@test A^0 == mpeye(Float64, 2,2)
@test A * A == A^2 == MP([10.0 7.0; 11.0 10.0])
@test A * A * A == A^3 == MP([14.0 13.0; 17.0 14.0])
# TODO A^-1

C = [MP(3.0) 4.0; 5.0 6.0]
D = MP([1.0 2; 3 4])
@test C + D == D + C == MP([3.0 4.0; 5.0 6.0]) == C

# Trace
# TODO shall we forbid trace of non squared matrices ?
@test mptrace([]) == mp0
@test mptrace(mpeye(Float64, 2,2)) == tr(mpeye(Float64, 2,2)) == 0.0
@test mptrace(mpeye(Float64, 2,5)) == 0.0
@test mptrace(mpzeros(Float64, 2,2)) == mptrace(full(mpzeros(Float64, 2,2))) == tr(mpzeros(Float64, 2,2)) == mp0
@test mptrace(mpzeros(Float64, 2,5)) == mptrace(full(mpzeros(Float64, 2,5))) == mp0
@test mptrace([1.0 2.0; 3.0 4.0]) == tr(MP([1.0 2.0; 3.0 4.0])) == MP(1.0) + MP(4.0) == MP(4.0)
@test mptrace(sparse([1.0 2.0; 3.0 4.0])) == mptrace(mpsparse([1.0 2.0; 3.0 4.0])) == MP(4.0)

# ==============================================================================
# Max-plus norm
# ==============================================================================

@test mpnorm(MP([1 20 2;30 400 4;4 50 10])) == MP(400 - 1) == MP(399)
@test mpnorm(mpsparse([1 20 2;30 400 4;4 50 10])) == MP(400 - 1) == MP(399)
@test mpnorm([mp0 1; 10 mp0]) == MP(10.0) - mp0 == mptop

# ==============================================================================
# Max-plus star
# ==============================================================================
@test mpstar(MP(1.0)) == mpstar(1.0) == mptop
@test mpstar(MP(-1.0)) == mpstar(-1.0) == mp1
@test mpstar(MP([1.0 2; 3 4])) == mpstar([1.0 2; 3 4]) == [mptop mptop; mptop mptop]
A=MP([-3.0 -2; -1 0]); B = mpstar(A)
@test B == mpeye(Float64, 2,2) + A
@test B == B + A^2
@test B == B + A^3

@test mpstar(mpeye(Float64, 2,2)) == mpeye(Float64, 2,2)

#TODO sparse
#sA = mpsparse(A); sB = mpstar(sA)
#@test B == mpeye(Float64, 2,2) + A
#@test B == B + A^2
#@test B == B + A^3

# ==============================================================================
# TODO @test C / C == MP([0.0 -2.0; 2.0 0.0])
# ==============================================================================


# ==============================================================================
# Max-plus display
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
