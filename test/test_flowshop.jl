# ==============================================================================
# Flowshop unit tests.
#
# Reference workload: the `complex.flowshop` example shipped with the
# TimedPetriNetEditor companion project. It describes a shop with 8 machines
# and 6 pieces, each machine class and each piece class holding exactly one
# token (single resource, single piece in flight per class).
#
# The matrix below is the verbatim transcription of the processing times
# in `complex.flowshop` (-Inf = no task). We inline it so that the test
# stays self-contained and does not depend on an external file outside this
# repository.
# ==============================================================================

using MaxPlus, Test, SparseArrays, LinearAlgebra

const _NA = -Inf

# Processing times: rows = machines, columns = pieces.
const COMPLEX_E = MP.([
    2.0  3.9  0.95 1.1  0.7  1.4 ;
    _NA  _NA  2.0  1.2  _NA  1.7 ;
    3.7  _NA  2.2  _NA  6.4  _NA ;
    _NA  _NA  2.0  _NA  1.0  1.0 ;
    1.7  3.1  3.0  _NA  1.3  _NA ;
    0.5  3.2  4.3  1.9  1.6  0.4 ;
    1.0  1.0  1.0  1.0  1.0  1.0 ;
    1.5  1.5  1.5  1.2  1.2  1.2 ;
])

# One token per machine class and per piece class (cf. `nm:` and `np:`
# lines in the original .flowshop file).
const COMPLEX_NM = ones(Int, 8)
const COMPLEX_NP = ones(Int, 6)

# Reference critical cycle, obtained by closing the feedback loops with
# `mpshift`, switching to explicit state-space form and running Howard.
# This is the (max,+) eigenvalue, i.e. the throughput period of the shop.
const COMPLEX_CRITICAL_CYCLE = 16.95

# Helper: build the explicit closed-loop state matrix of an E-flowshop
# parametrised by (nm, np). Mirrors what `flowshop_simu` does internally
# but exposes the intermediate Howard input so we can assert on it.
function _closed_loop_state(E, nm, np)
    s = flowshop(E)
    fbm = mpshift(nm[1], 0)
    for i in 2:length(nm)
        fbm = fbm | mpshift(nm[i], 0)
    end
    fbp = mpshift(np[1], 0)
    for i in 2:length(np)
        fbp = fbp | mpshift(np[i], 0)
    end
    sb  = s / (fbp | fbm)
    sbs = explicit(sb)
    return sparse(sbs.A)
end

# ------------------------------------------------------------------------------
# Shape checks: flowshop must produce a well-formed MPSysLin whose dimensions
# match the (nmach + npiece) input/output count and the (nmach * npiece)
# internal state count.
# ------------------------------------------------------------------------------
@testset "flowshop: shape on complex.flowshop" begin
    nmach, npiece = size(COMPLEX_E)
    s = flowshop(COMPLEX_E)
    @test size(s.A) == (nmach * npiece, nmach * npiece)
    @test size(s.B) == (nmach * npiece, nmach + npiece)
    @test size(s.C) == (nmach + npiece, nmach * npiece)
    @test s.A isa SparseMatrixCSC{<:MP}
    @test s.B isa SparseMatrixCSC{<:MP}
    @test s.C isa SparseMatrixCSC{<:MP}
end

# ------------------------------------------------------------------------------
# Critical cycle: Howard converges to a single strongly connected component
# whose eigenvalue is the shop throughput period. With these specific
# processing times the period is 16.95 time units.
# ------------------------------------------------------------------------------
@testset "flowshop: critical cycle on complex.flowshop" begin
    Sa  = _closed_loop_state(COMPLEX_E, COMPLEX_NM, COMPLEX_NP)
    chi = howard(Sa)

    # Howard returns one eigenvalue per node; on a strongly connected
    # graph they must all be equal to the critical cycle time.
    @test !isempty(chi.eigenvalues)
    eigs = plustimes.(chi.eigenvalues)
    @test all(isapprox.(eigs, COMPLEX_CRITICAL_CYCLE; atol = 1e-9))
    @test isapprox(maximum(eigs), COMPLEX_CRITICAL_CYCLE; atol = 1e-9)

    # Single strongly connected component (the closed-loop graph is one
    # big circuit), and Howard converges in very few iterations.
    @test chi.components == 1
    @test chi.iterations ≥ 1
end

# ------------------------------------------------------------------------------
# End-to-end simulation: `flowshop_simu` must run on the same data and
# return outputs shaped (ntimes, nmach + npiece) — see `flowshop_simu`.
# We only check the dimensions and the cycle time it reports (same chi as
# the standalone Howard run above).
# ------------------------------------------------------------------------------
@testset "flowshop_simu: complex.flowshop end-to-end" begin
    s = flowshop(COMPLEX_E)
    nmach, npiece = size(COMPLEX_E)
    ntimes = 5
    u = zeros(Float64, nmach + npiece, ntimes)
    chi, y = flowshop_simu(s, COMPLEX_NM, COMPLEX_NP, u)

    @test size(y, 1) == ntimes
    @test size(y, 2) == nmach + npiece
    @test !isempty(chi.eigenvalues)
    @test isapprox(plustimes(chi.eigenvalues[1]), COMPLEX_CRITICAL_CYCLE; atol = 1e-9)
end
