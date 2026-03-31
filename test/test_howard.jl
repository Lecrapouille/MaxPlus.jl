# ==============================================================================
# Max-Plus Howard algorithm
# ==============================================================================

# Irreducible matrix, only 1 eigenvalue. l is constant.
S = sparse(MP([1 2; 3 4]))
λ,v = mpeigen(S)
@test (λ,v) == (MP[4, 4], MP[2, 4])
@test S * v == λ[1] * v
@test S * v == λ[2] * v

# Two blocks diagonal matrix. The entries of l take two values.
S = sparse([mp0 2 mp0; mp1 mp0 mp0; mp0 mp0 2])
λ,v = mpeigen(S)
@test (λ,v) == (MP[1, 1, 2], MP[1, 0, 2])
@test diagm(λ) == [1 mp0 mp0; mp0 1 mp0; mp0 mp0 2]
@test S / diagm(λ) == [mp0 1 mp0; -1 mp0 mp0; mp0 mp0 0]
@test (S / diagm(λ)) * v == v

# Block triangular matrix with 2 eigenvalues.
# l is constant only one eigen value is found.
S = sparse([1 1; mp0 2])
λ,v = mpeigen(S)
@test (λ,v) == (MP[2, 2], MP[1, 2])
@test S * v == λ[1] * v
@test S * v == λ[2] * v
@test (S * [0; mp0]) == (MP(1) * [0; mp0])
# But 1 also en eigen value
@test (S * [0; mp0]) == (MP(1) * [0; mp0])

# Block triangular matrix with 1 eigenvalue
# l is not constant λ(1) is eigen value
# with eigen vector [v(1);0]
S = sparse([2 1; mp0 mp1])
λ,v = mpeigen(S)
@test (λ,v) == (MP[2, 0], MP[2, 0])
@test diagm(λ) == [2 mp0; mp0 0]
@test S / diagm(λ) == [0 1; mp0 0]
@test (S / diagm(λ)) * v == v
@test (S * [v[1]; mp0]) == (λ[1] * [v[1]; mp0])

# Book "Synchronization and Linearity An Algebra for Discrete Event Systems"
S = MP([3.0 7; 2 4])
λ,v = mpeigen(sparse(S))
@test (λ,v) == (MP[4.5, 4.5], MP[6.5, 4])
@test S * v == λ[1] * v
@test S * v == λ[2] * v

# previous matrix modified
S = MP([3.0 7; 2 -Inf])
λ,v = mpeigen(sparse(S))
@test (λ,v) == (MP[4.5, 4.5], MP[4.5, 2])
@test S * v == λ[1] * v
@test S * v == λ[2] * v

# previous matrix modified
S = MP([3.0 -Inf; 2 4])
λ,v = mpeigen(sparse(S))
@test (λ,v) == (MP[3, 4], MP[3, 4])
@test S * v != λ[1] * v
@test S * v != λ[2] * v

# Demo Scilab with fixed random seed for reproducibility
Random.seed!(42)
S = MP(sprand(10,10,0.0005))+0.001*sparse(eye(MP, 10,10))
λ,v = mpeigen(S)
@test S * v == λ[1] * v
@test S * v == λ[2] * v

# Test dense matrix interface
A = MP([3.0 7; 2 4])
λ,v = mpeigen(A)
@test (λ,v) == (MP[4.5, 4.5], MP[6.5, 4])

# ==============================================================================
# Semi-Howard algorithm tests
# ==============================================================================

# Simple test: when all tau = 1, semihoward should give same result as howard
S = sparse(MP([1.0 2; 3 4]))
Tau = sparse(MP([1.0 1.0; 1.0 1.0]))  # All ones = standard howard
r = semihoward(S, Tau)
λh, vh = mpeigen(S)
@test r.eigenvalues == λh
@test r.eigenvectors == vh

# Test with varying tau values
# Cycle time = sum(weights) / sum(tau) instead of sum(weights) / length
S = sparse(MP([2.0 mp0; mp0 3.0]))  # Two diagonal 1-cycles
Tau = sparse(MP([2.0 mp0; mp0 1.0]))  # tau[1]=2, tau[2]=1
r = semihoward(S, Tau)
# Node 1: cycle (1->1), weight=2, tau=2, lambda = 2/2 = 1
# Node 2: cycle (2->2), weight=3, tau=1, lambda = 3/1 = 3
@test r.eigenvalues == MP[1.0, 3.0]
