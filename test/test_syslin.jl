# ==============================================================================
# Max-plus system linear
# ==============================================================================

mpsparse(S) = MP(sparse(S))

# ==============================================================================
# Creation of max-plus system linear

### First system
S1 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]),
              mpsparse([0;0;0]),
              mpsparse([0 0 0]),
              sparse(eye(MP, 3,3)),
              zeros(MP, 3,1))

@test S1.A == [1 2 3; 4 5 6; 7 8 9]
@test_broken S1.B == reshape([0;0;0], length([0;0;0]), 1)
@test S1.C == [0 0 0]
@test S1.D == eye(MP, 3,3)
@test S1.x0 == zeros(MP, 3,1)

S2 = MPSysLin(mpsparse([1 2 3; 4 5 6; 7 8 9]),
              mpsparse([0;0;0]),
              mpsparse([0 0 0]),
              mpsparse(eye(MP, 3,3)))
S3 = MPSysLin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]))
S4 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), eye(MP, 3,3), MP([0; 0; 0]))
S5 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), eye(MP, 3,3))
S6 = MPSysLin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]))

@test S1.A == S2.A == S3.A == S4.A == S5.A == S6.A
@test S1.B == S2.B == S3.B == S4.B == S5.B == S6.B
@test S1.C == S2.C == S3.C == S4.C == S5.C == S6.C
#@test S1.D == S2.D == S3.D == S4.D == S5.D == S6.D
@test S1.x0 == S2.x0 == S3.x0 == S5.x0 == S6.x0
@test S4.x0 == ones(MP, 3)

### Second system
S7 = MPSysLin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(eye(MP, 2,2)), zeros(MP, 2,1))

@test S7.A == [1 2; 3 4]
@test_broken S7.B == reshape([0;0], length([0;0]), 1)
@test S7.C == [0 0]
@test S7.D == eye(MP, 2,2)
@test S7.x0 == zeros(MP, 2,1)

S8 = MPSysLin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(eye(MP, 2,2)))
S9 = MPSysLin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]))
S10 = MPSysLin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), eye(MP, 2,2), MP([0; 0]))
S11 = MPSysLin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), eye(MP, 2,2))
S12 = MPSysLin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]))

@test S7.A == S8.A == S9.A == S10.A == S11.A == S12.A
@test S7.B == S8.B == S9.B == S10.B == S11.B == S12.B
@test S7.C == S8.C == S9.C == S10.C == S11.C == S12.C
@test S7.D == S8.D
@test_broken S9.D == S10.D
@test S10.D == S11.D
@test_broken S11.D == S12.D
@test S7.x0 == S8.x0 == S9.x0 == S11.x0 == S12.x0
@test S10.x0 == ones(MP, 2)

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
@test_broken S13.x0 == reshape(x0, length(x0), 1)
@test S14.x0 == S15.x0 == zeros(MP, 2, 1)

# ==============================================================================
# Simulation. Note here we should call the explicit systems but here
# they are the same than implicit.

res = mpsimul(S7, MP(1:10), false)
@test res == MP([37])

res = mpsimul(S7, MP(1:10), true)
@test res == MP([1 5 9 13 17 21 25 29 33 37])
@test res == mpsimul(S8, MP(1:10), true)
@test res == mpsimul(S9, MP(1:10), true)
@test res == mpsimul(S11, MP(1:10), true)
@test res == mpsimul(S12, MP(1:10), true)

@test mpsimul(S10, MP(1:10), true) == MP([4  8  12  16  20  24  28  32  36  40])

# ==============================================================================
# Composition of systems

S1 = MPSysLin(MP([1.0 2; 3 4]), MP([0.0; 0]), MP([0.0 0]), eye(MP, 2,2))
S2 = MPSysLin(MP([1.0 2 3; 4 5 6; 7 8 9]), MP([0.0;0;0]), MP([0.0 0 0]), eye(MP, 3,3))

# Diagonal composition
S = S1 | S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == [mp1 mp0; mp1 mp0; mp0 mp1; mp0 mp1; mp0 mp1]
@test S.C == [mp1 mp1 mp0 mp0 mp0; mp0 mp0 mp1 mp1 mp1]
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5)

# Parallel composition
S = S1 + S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == ones(MP, 5,1)
@test S.C == ones(MP, 1,5)
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5)

# Series composition
S = S1 * S2

@test S.A == [1.0 2 3 mp0 mp0; 4 5 6 mp0 mp0; 7 8 9 mp0 mp0; mp0 mp0 mp0 1 2; mp0 mp0 mp0 3 4]
@test S.B == [ones(MP, 3,1); dense(zeros(MP, 2,1))]
@test S.C == [dense(zeros(MP, 1,3)) ones(MP, 1,2)]
D = eye(MP, 5,5); D[4:5, 1:3] .= mp1;
@test S.D == D
@test S.x0 == zeros(MP, 5,1)

# Input in common
S = [S1; S2]

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == ones(MP, 5,1)
@test S.C == [mp1 mp1 mp0 mp0 mp0; mp0 mp0 mp1 mp1 mp1]
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5,1)

# Output addition FIXME syntax S = [S1, S2] is not working
S = [S1 S2]

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == [mp1 mp0; mp1 mp0; mp0 mp1; mp0 mp1; mp0 mp1]
@test S.C == ones(MP, 1,5)
@test S.D == eye(MP, 5,5)
@test S.x0 == zeros(MP, 5,1)

# Feedback composition FIXME syntax /.
S = S1 / S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == [ones(MP, 2,1); dense(zeros(MP, 3,1))]
@test S.C == [ones(MP, 1,2) dense(zeros(MP, 1,3))]
D = eye(MP, 5,5); D[3:5, 1:2] .= mp1; D[1:2, 3:5] .= mp1;
@test S.D == D
@test S.x0 == zeros(MP, 5,1)

# ==============================================================================
#

S1 = MPSysLin(MP([1.0 2; 3 4]), MP([0.0; 0]), MP([0.0 0]), eye(MP, 2,2))

S = 5 * S1
@test S.A == [1.0 2; 3 4]
@test_broken S.B == reshape([0.0; 0], length([0;0]), 1)
@test S.C == [5.0 5]
@test S.D == eye(MP, 2,2)
@test S.x0 == zeros(MP, 2,1)
@test S == MP(5.0) * S1

S = S1 * 5
@test S.A == [1.0 2; 3 4]
@test_broken S.B == reshape([5.0; 5], length([0;0]), 1)
@test S.C == [0.0 0]
@test S.D == eye(MP, 2,2)
@test S.x0 == zeros(MP, 2,1)
@test S1 * 5 == S1 * MP(5.0)

# ==============================================================================
# Implicit to Explicit

S1 = MPSysLin(MP([1.0 2; 3 4]), MP([0.0; 0]), MP([0.0 0]), eye(MP, 2,2))
S = mpexplicit(S1)

@test S.A == S1.A
@test S.B == S1.B
@test S.C == S1.C
@test S.D == eye(MP, 2,2)
@test S.x0 == S1.x0

S2 = MPSysLin(MP([1.0 2 3; 4 5 6; 7 8 9]), MP([0.0;0;0]), MP([0.0 0 0]), eye(MP, 3,3))
S = mpexplicit(S2)

@test S.A == S2.A
@test S.B == S2.B
@test S.C == S2.C
@test S.D == eye(MP, 3,3)
@test S.x0 == S2.x0
