# ==============================================================================
# Max-plus system linear
# ==============================================================================

# ==============================================================================
# Creation of max-plus system linear

### First system - using dense matrices
S1 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]),
              MP([1.0; 1.0; 1.0]),
              MP([1.0 1.0 1.0]),
              eye(MP, 3,3),
              zeros(MP, 3,1))

@test S1.A == [1 2 3; 4 5 6; 7 8 9]
@test S1.B == reshape(MP([1.0;1.0;1.0]), 3, 1)
@test S1.C == MP([1.0 1.0 1.0])
@test S1.D == eye(MP, 3,3)
@test S1.x0 == zeros(MP, 3,1)

# Test equivalent constructors
S4 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]), MP([1.0;1.0;1.0]), MP([1.0 1.0 1.0]), eye(MP, 3,3), MP([1.0; 1.0; 1.0]))
S5 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]), MP([1.0;1.0;1.0]), MP([1.0 1.0 1.0]), eye(MP, 3,3))
S6 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]), MP([1.0;1.0;1.0]), MP([1.0 1.0 1.0]))

@test S1.A == S4.A == S5.A == S6.A
@test S1.B == S4.B == S5.B == S6.B
@test S1.C == S4.C == S5.C == S6.C
@test S1.x0 == S5.x0 == S6.x0
@test S4.x0 == reshape(MP([1.0; 1.0; 1.0]), 3, 1)  # x0 is passed as [1.0; 1.0; 1.0]

### Second system
S7 = MPSysLin(MP([1 2; 3 4]), MP([1.0; 1.0]), MP([1.0 1.0]), eye(MP, 2,2), zeros(MP, 2,1))

@test S7.A == [1 2; 3 4]
@test S7.B == reshape(MP([1.0;1.0]), 2, 1)
@test S7.C == MP([1.0 1.0])
@test S7.D == eye(MP, 2,2)
@test S7.x0 == zeros(MP, 2,1)

S10 = MPSysLin(MP([1 2; 3 4]), MP([1.0; 1.0]), MP([1.0 1.0]), eye(MP, 2,2), MP([1.0; 1.0]))
S11 = MPSysLin(MP([1 2; 3 4]), MP([1.0; 1.0]), MP([1.0 1.0]), eye(MP, 2,2))
S12 = MPSysLin(MP([1 2; 3 4]), MP([1.0; 1.0]), MP([1.0 1.0]))

@test S7.A == S10.A == S11.A == S12.A
@test S7.B == S10.B == S11.B == S12.B
@test S7.C == S10.C == S11.C == S12.C
@test S7.x0 == S11.x0 == S12.x0
@test S10.x0 == reshape(MP([1.0; 1.0]), 2, 1)  # x0 is passed as [1.0; 1.0]

### Third system
A = MP([1 2; 3 4]); B = MP([5 6; 7 8]); C = MP([9 10; 11 12]); D = MP([13 14; 15 16])
x0 = MP([17; 18])

S13 = MPSysLin(A,B,C,D,x0)
S14 = MPSysLin(A,B,C,D)
S15 = MPSysLin(A,B,C)
@test S13.A == S14.A == S15.A == A
@test S13.B == S14.B == S15.B == B
@test S13.C == S14.C == S15.C == C
@test S13.D == S14.D == D
@test S13.x0 == reshape(x0, length(x0), 1)
@test S14.x0 == S15.x0 == zeros(MP, 2, 1)

# ==============================================================================
# Simulation

# Create a simple system for simulation tests
Ssim = MPSysLin(MP([1.0 2; 3 4]), MP([1.0; 1.0]), MP([1.0 1.0]), eye(MP, 2,2))
res = simul(Ssim, MP(1:10), false)
@test length(res) > 0

res = simul(Ssim, MP(1:10), true)
@test size(res, 1) == 10

# ==============================================================================
# Composition of systems

S1 = MPSysLin(MP([1.0 2; 3 4]), MP([1.0; 1.0]), MP([1.0 1.0]), eye(MP, 2,2))
S2 = MPSysLin(MP([1.0 2 3; 4 5 6; 7 8 9]), MP([1.0;1.0;1.0]), MP([1.0 1.0 1.0]), eye(MP, 3,3))

# Diagonal composition
S = S1 | S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test size(S.B) == (5, 2)
@test size(S.C) == (2, 5)
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5, 1)

# Parallel composition
S = S1 + S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test size(S.B) == (5, 1)
@test size(S.C) == (1, 5)
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5, 1)

# Series composition
S = S1 * S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test size(S.B) == (5, 1)
@test size(S.C) == (1, 5)
@test S.D == [eye(MP, 2, 2) spzeros(MP, 2, 3); S2.B * S1.C eye(MP, 3, 3)]
@test S.x0 == zeros(MP, 5, 1)

# Input in common (vcat)
S = [S1; S2]

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test size(S.B) == (5, 1)
@test size(S.C) == (2, 5)
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5,1)

# Output addition (hcat) - note: [S1, S2] creates a vector, [S1 S2] uses hcat
S = [S1 S2]

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test size(S.B) == (5, 2)
@test size(S.C) == (1, 5)
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5,1)

# Feedback composition
S = S1 / S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test size(S.B) == (5, 1)
@test size(S.C) == (1, 5)
@test S.x0 == zeros(MP, 5,1)

# ==============================================================================
# Scalar multiplication

S1 = MPSysLin(MP([1.0 2; 3 4]), MP([2.0; 2.0]), MP([3.0 3.0]), eye(MP, 2,2))

S = 5 * S1
@test S.A == [1.0 2; 3 4]
@test S.B == reshape(MP([2.0; 2.0]), 2, 1)
@test S.C == MP([8.0 8.0])  # 3 + 5 = 8 in max-plus
@test S.D == eye(MP, 2,2)
@test S.x0 == zeros(MP, 2,1)
@test S == MP(5.0) * S1

S = S1 * 5
@test S.A == [1.0 2; 3 4]
@test S.B == reshape(MP([7.0; 7.0]), 2, 1)  # 2 + 5 = 7 in max-plus
@test S.C == MP([3.0 3.0])
@test S.D == eye(MP, 2,2)
@test S.x0 == zeros(MP, 2,1)
@test S1 * 5 == S1 * MP(5.0)

# ==============================================================================
# Implicit to Explicit
# explicit converts implicit max-plus systems to explicit state-space form.
# If D is the identity (I), the system is already explicit and A, B, C are unchanged.

S1 = MPSysLin(MP([1.0 2; 3 4]), MP([2.0; 2.0]), MP([3.0 3.0]), eye(MP, 2, 2))
# The explicit form should eliminate the implicit part: Sex.D should be zero (spzeros)
Sex = explicit(S1)
@test size(Sex.A, 1) == size(Sex.A, 2)   # A should be square
@test Sex.D == spzeros(MP, size(Sex.A, 1), size(Sex.A, 1)) # D should be zero (explicit)

Si1 = implicit(S1)
@test Si1.A == S1.A && Si1.B == S1.B && Si1.C == S1.C
@test Si1.D == eye(MP, 2, 2)
@test implicit(S1) == Si1
u1 = MP([1.0, 2.0])
@test simul(S1, u1, false) == simul(Si1, u1, false)

# Conversion between dense and sparse system forms
Sf = mpfull(Sex)
@test Sf.A isa Matrix{MP}
Ss = mpsparse(Sf)
@test Ss.A isa SparseMatrixCSC{MP}

# Smoke tests for shift, flowshop, and flowshop_simu functions
sh = mpshift(2, 0)
@test size(sh.A, 1) == 3
@test sh.x0 isa SparseMatrixCSC{MP}
Efs = MP.([2.0 1; 3.0 2])
sfs = flowshop(Efs)
@test size(sfs.B, 2) == size(Efs, 1) + size(Efs, 2)
@test sfs.A isa SparseMatrixCSC{MP}
@test sfs.B isa SparseMatrixCSC{MP}
@test sfs.C isa SparseMatrixCSC{MP}
@test sfs.D isa SparseMatrixCSC{MP}
@test sfs.x0 isa SparseMatrixCSC{MP}

Tg, Ng = flowshop_graph(Efs, ones(2), ones(2))
@test Tg isa SparseMatrixCSC{MP}
@test Ng isa SparseMatrixCSC{MP}
@test size(Tg) == size(Ng)
@test size(Tg, 1) ≥ 1

# Simulate the flowshop system
chi, yfs = flowshop_simu(sfs, [1, 1], [1, 1], zeros(Float64, 4, 3))
@test size(yfs, 1) == 3
@test size(yfs, 2) == 4
