# ==============================================================================
# max-plus system linear
# ==============================================================================
A = MP([1 2; 3 4]); B = MP([5 6; 7 8]); C = MP([9 10; 11 12]); D = MP([13 14; 15 16])
x0 = MP([17; 18])

S1 = mpsyslin(A,B,C,D,x0)
S2 = mpsyslin(A,B,C,D)
S3 = mpsyslin(A,B,C)
@test S1.A == S2.A == S3.A == A
@test S1.B == S2.B == S3.B == B
@test S1.C == S2.C == S3.C == C
@test S1.D == S2.D == D
@test S3.D == full(mpzeros(Int64, 2, 2))
#@test S1.x0 == x0
#@test S2.x0 == S3.x0 == full(mpzeros(Int64, 1, 2))

S4 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]), mpsparse(mpeye(Int64, 3,3)), mpzeros(Int64, 3,1))
S1 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(mpeye(Int64, 2,2)), mpzeros(Int64, 2,1))

#S4 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]), mpsparse(mpeye(Int64, 3,3)))
#S1 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(mpeye(Int64, 2,2)))

#S4 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]))
#S1 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]))

#S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), mpeye(Int64, 3,3), MP([0; 0; 0]))
#S1 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), mpeye(Int64, 2,2), MP([0; 0]))

#S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), mpeye(Int64, 3,3))
#S1 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), mpeye(Int64, 2,2))

#S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]))
#S1 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]))

# S1+S4 S4+S1
# S1*S4 S1\S4 S1|S4 [S1 S4] [S1; S4]

# ==============================================================================
# TODO deprecated
# ==============================================================================
#A = MP([1.0 2.0; 3.0 4.0])
#x0 = MP([1.0; 2.0])
#x1 = mpsyslin(A, x0, 1)
#x2 = mpsyslin(A, x1, 1)
#x3 = mpsyslin(A, x2, 1)
#X = mpsimul(A, x0, 3, history=true)
#@test [x0 x1 x2 x3] == X
