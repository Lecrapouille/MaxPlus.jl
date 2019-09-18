# ==============================================================================
# Max-plus system linear
# ==============================================================================

# ==============================================================================
# Creation of max-plus system linear

### First system
S1 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]), mpsparse(mpeye(Int64, 3,3)), mpzeros(Int64, 3,1))

@test S1.A == [1 2 3; 4 5 6; 7 8 9]
@test S1.B[:,1] == [0;0;0]
@test S1.C == [0 0 0]
@test S1.D == mpeye(Int64, 3,3)
@test S1.x0 == mpzeros(Int64, 3,1)

S2 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]), mpsparse(mpeye(Int64, 3,3)))
S3 = mpsyslin(mpsparse([1 2 3; 4 5 6; 7 8 9]), mpsparse([0;0;0]), mpsparse([0 0 0]))
S4 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), mpeye(Int64, 3,3), MP([0; 0; 0]))
S5 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]), mpeye(Int64, 3,3))
S6 = mpsyslin(MP([1 2 3; 4 5 6; 7 8 9]), MP([0;0;0]), MP([0 0 0]))

@test S1.A == S2.A == S3.A == S4.A == S5.A == S6.A
@test S1.B == S2.B == S3.B == S4.B == S5.B == S6.B
@test S1.C == S2.C == S3.C == S4.C == S5.C == S6.C
@test S1.D == S2.D == S3.D == S4.D == S5.D == S6.D
@test S1.x0 == S2.x0 == S3.x0 == S4.x0 == S5.x0 == S6.x0

### Second system
S7 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(mpeye(Int64, 2,2)), mpzeros(Int64, 2,1))

@test S7.A == [1 2; 3 4]
@test S7.B[:,1] == [0;0]
@test S7.C == [0 0]
@test S7.D == mpeye(Int64, 2,2)
@test S7.x0 == mpzeros(Int64, 2,1)

S8 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]), mpsparse(mpeye(Int64, 2,2)))
S9 = mpsyslin(mpsparse([1 2; 3 4]), mpsparse([0; 0]), mpsparse([0 0]))
S10 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), mpeye(Int64, 2,2), MP([0; 0]))
S11 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]), mpeye(Int64, 2,2))
S12 = mpsyslin(MP([1 2; 3 4]), MP([0; 0]), MP([0 0]))

@test S7.A == S8.A == S9.A == S10.A == S11.A == S12.A
@test S7.B == S8.B == S9.B == S10.B == S11.B == S12.B
@test S7.C == S8.C == S9.C == S10.C == S11.C == S12.C
@test S7.D == S8.D == S9.D == S10.D == S11.D == S12.D
@test S7.x0 == S8.x0 == S9.x0 == S10.x0 == S11.x0 == S12.x0

### Third system
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
@test S1.x0[:,1] == x0
@test S2.x0 == S3.x0 == mpzeros(Int64, 2, 1)

# ==============================================================================
# Composition of systems

S1 = mpsyslin(MP([1.0 2; 3 4]), MP([0.0; 0]), MP([0.0 0]), mpeye(Float64, 2,2))
S2 = mpsyslin(MP([1.0 2 3; 4 5 6; 7 8 9]), MP([0.0;0;0]), MP([0.0 0 0]), mpeye(Float64, 3,3))

# Diagonal composition
S = S1 | S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == [mp1 mp0; mp1 mp0; mp0 mp1; mp0 mp1; mp0 mp1]
@test S.C == [mp1 mp1 mp0 mp0 mp0; mp0 mp0 mp1 mp1 mp1]
@test S.D == mpeye(Float64, 5,5)
@test S.x0 == mpzeros(Float64, 5,1)

# Parallel composition
S = S1 + S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == mpones(Float64, 5,1)
@test S.C == mpones(Float64, 1,5)
@test S.D == mpeye(Float64, 5,5)
@test S.x0 == mpzeros(Float64, 5,1)

# Series composition
S = S1 * S2

@test S.A == [1.0 2 3 mp0 mp0; 4 5 6 mp0 mp0; 7 8 9 mp0 mp0; mp0 mp0 mp0 1 2; mp0 mp0 mp0 3 4]
@test S.B == [mpones(Float64, 3,1); dense(mpzeros(Float64, 2,1))]
@test S.C == [dense(mpzeros(Float64, 1,3)) mpones(Float64, 1,2)]
D = mpeye(Float64, 5,5); D[4:5, 1:3] .= mp1;
@test S.D == D
@test S.x0 == mpzeros(Float64, 5,1)

# Input in common
S = [S1; S2]

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == mpones(Float64, 5,1)
@test S.C == [mp1 mp1 mp0 mp0 mp0; mp0 mp0 mp1 mp1 mp1]
@test S.D == mpeye(Float64, 5,5)
@test S.x0 == mpzeros(Float64, 5,1)

# Output addition FIXME syntax S = [S1, S2] is not working
S = [S1 S2]

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == [mp1 mp0; mp1 mp0; mp0 mp1; mp0 mp1; mp0 mp1]
@test S.C == mpones(Float64, 1,5)
@test S.D == mpeye(Float64, 5,5)
@test S.x0 == mpzeros(Float64, 5,1)

# Feedback composition FIXME syntax /.
S = S1 / S2

@test S.A == [1.0 2 mp0 mp0 mp0; 3 4 mp0 mp0 mp0; mp0 mp0 1 2 3; mp0 mp0 4 5 6; mp0 mp0 7 8 9]
@test S.B == [mpones(Float64, 2,1); dense(mpzeros(Float64, 3,1))]
@test S.C == [mpones(Float64, 1,2) dense(mpzeros(Float64, 1,3))]
D = mpeye(Float64, 5,5); D[3:5, 1:2] .= mp1; D[1:2, 3:5] .= mp1;
@test S.D == D
@test S.x0 == mpzeros(Float64, 5,1)

# ==============================================================================
# Extraction

(A,B,C,D,x0) = full(S1 / S2)
@test A == S.A
@test B == S.B
@test C == S.C
@test D == S.D
@test x0 == S.x0

(A,B,C,D,x0) = dense(S1 / S2)
@test A == S.A
@test B == S.B
@test C == S.C
@test D == S.D
@test x0 == S.x0

# ==============================================================================
# Simulation

S1 = mpsyslin(MP([1.0 2; 3 4]), MP([0.0; 0]), MP([0.0 0]), mpeye(Float64, 2,2))
res = mpsimul(S1, MP(1:10), true)
@test res == MP([1.0 5 9 13 17 21 25 29 33 37])

res = mpsimul(S1, MP(1:10), false)
@test res == MP{Float64}[37.0]
