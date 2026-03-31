################################################################################
# Flowshop (max,+)-linear systems — porting from flowshop.sci, shift.sci,
# flowshop_simu.sci (ScicosLab).
################################################################################

# ==============================================================================
# shift.sci — delay n on event indices, t on dates.
function mpshift(n::Integer, t::Real)
    n >= 1 || error("mpshift: n must be >= 1")
    na = n + 1
    Ii = 1:n
    Jj = 2:na
    Vv = fill(one(MP), n)
    A = SparseArrays.sparse(Ii, Jj, Vv, na, na)
    B = spzeros(MP, na, 1)
    B[na, 1] = one(MP)
    C = spzeros(MP, 1, na)
    C[1, 1] = MP(t)
    D = spzeros(MP, na, na)
    x0 = spzeros(MP, na, 1)
    return MPSysLin(A, B, C, D, x0)
end

# ==============================================================================
# flowshop.sci — E[machine, piece]: (max-plus) duration; ε = no task.
function flowshop(E::AbstractMatrix{MP})
    (nmach, npiece) = size(E)
    nstate = nmach * npiece
    nio = nmach + npiece
    A = spzeros(MP, nstate, nstate)
    D = spzeros(MP, nstate, nstate)
    B = spzeros(MP, nstate, nio)
    C = spzeros(MP, nio, nstate)
    x0 = spzeros(MP, nstate, 1)
    l = zeros(Int, nmach)
    d = zeros(Int, npiece)
    for i in 1:npiece
        for j in 1:nmach
            ij = i + (j - 1) * npiece
            if E[j, i] == mp0
                if l[j] != 0
                    l[j] += 1
                end
                if d[i] != 0
                    d[i] += 1
                end
            else
                if l[j] != 0
                    D[ij, ij - l[j]] = E[j, i - l[j]]
                end
                if d[i] != 0
                    D[ij, ij - d[i] * npiece] = E[j - d[i], i]
                end
                if d[i] == 0
                    B[ij, i] = one(MP)
                end
                d[i] = 1
                if l[j] == 0
                    B[ij, j + npiece] = one(MP)
                end
                l[j] = 1
            end
        end
        ij = i + (nmach - 1) * npiece
        ijc = ij - (d[i] - 1) * npiece
        C[i, ijc] = E[nmach - d[i] + 1, i]
    end
    for j in 1:nmach
        jp = j * npiece - l[j] + 1
        C[j + npiece, jp] = E[j, npiece - l[j] + 1]
    end
    return MPSysLin(A, B, C, D, x0)
end

flowshop(E::AbstractMatrix{<:Real}) = flowshop(map(MP, E))

# ==============================================================================
# flowshop_graph.sci — uniquement les matrices T et N (pas de graphe ni GUI).
# Retourne les mêmes objets qu’en Scilab après `T = sparse(T'); N = sparse(N')`.

function flowshop_graph(
    E::AbstractMatrix{MP},
    m::AbstractVector{<:Real},
    p::AbstractVector{<:Real}
)
    (nmach, npiece) = size(E)
    length(m) == nmach || error("flowshop_graph: length(m) must equal size(E,1)")
    length(p) == npiece || error("flowshop_graph: length(p) must equal size(E,2)")
    A = zeros(Int, nmach, npiece)
    l = zeros(Int, nmach)
    d = zeros(Int, npiece)
    bp = zeros(Int, npiece)
    bm = zeros(Int, nmach)
    nd = 0
    Tacc = Dict{Tuple{Int,Int}, MP}()
    Nacc = Dict{Tuple{Int,Int}, MP}()
    # Scilab T(r,c)=v ; ici on construit R = T' donc R[c,r]=v
    function tset!(r::Int, c::Int, v::MP)
        Tacc[(c, r)] = v
    end
    function nset!(r::Int, c::Int, v::MP)
        Nacc[(c, r)] = v
    end

    for i in 1:npiece
        for j in 1:nmach
            if E[j, i] == mp0
                A[j, i] = 0
                if l[j] != 0
                    l[j] += 1
                end
                if d[i] != 0
                    d[i] += 1
                end
            else
                nd += 1
                A[j, i] = nd
                if l[j] != 0
                    hc = A[j, i - l[j]]
                    tset!(hc, nd, E[j, i - l[j]])
                    nset!(hc, nd, one(MP))
                end
                if d[i] != 0
                    hc = A[j - d[i], i]
                    tset!(hc, nd, E[j - d[i], i])
                    nset!(hc, nd, one(MP))
                end
                if d[i] == 0
                    nd += 1
                    bp[i] = nd
                    nset!(nd, A[j, i], MP(Float64(p[i])))
                    tset!(nd, A[j, i], one(MP))
                end
                d[i] = 1
                if l[j] == 0
                    nd += 1
                    bm[j] = nd
                    nset!(nd, A[j, i], MP(Float64(m[j])))
                    tset!(nd, A[j, i], one(MP))
                end
                l[j] = 1
            end
        end
        nd += 1
        hc = A[nmach - d[i] + 1, i]
        tset!(hc, nd, E[nmach - d[i] + 1, i])
        nset!(hc, nd, one(MP))
        tset!(nd, bp[i], one(MP))
        nset!(nd, bp[i], one(MP))
    end
    for j in 1:nmach
        nd += 1
        hc = A[j, npiece - l[j] + 1]
        tset!(hc, nd, E[j, npiece - l[j] + 1])
        nset!(hc, nd, one(MP))
        tset!(nd, bm[j], one(MP))
        nset!(nd, bm[j], one(MP))
    end

    function mat_from_acc(acc::Dict{Tuple{Int,Int}, MP})
        I = Int[]
        J = Int[]
        V = MP[]
        for ((i, j), v) in acc
            push!(I, i)
            push!(J, j)
            push!(V, v)
        end
        return SparseArrays.sparse(I, J, V, nd, nd)
    end
    return mat_from_acc(Tacc), mat_from_acc(Nacc)
end

function flowshop_graph(
    E::AbstractMatrix{<:Real},
    m::AbstractVector{<:Real},
    p::AbstractVector{<:Real}
)
    flowshop_graph(map(MP, E), m, p)
end

# ==============================================================================
# flowshop_simu.sci — u: each column is a time step (as in Scilab);
# rows of u are the inputs (size = columns of B).
function flowshop_simu(
    s::MPSysLin,
    nm::AbstractVector{<:Integer},
    np::AbstractVector{<:Integer},
    u::AbstractMatrix{MP}
)
    nmach = length(nm)
    npiece = length(np)
    nt = size(u, 2)
    nt >= 1 || error("flowshop_simu: u must have at least one column (time)")
    fbm = mpshift(nm[1], 0)
    for i in 2:nmach
        fbm = fbm | mpshift(nm[i], 0)
    end
    fbp = mpshift(np[1], 0)
    for i in 2:npiece
        fbp = fbp | mpshift(np[i], 0)
    end
    sb = s / (fbp | fbm)
    sbs = explicit(sb)
    Sa = SparseArrays.sparse(sbs.A)
    chi = howard(Sa)
    y = simul(sbs, Matrix(transpose(u)), true)
    nout = nmach + npiece
    size(y, 2) == nout || error("flowshop_simu: inconsistent output dimensions")
    if isempty(chi.eigenvalues)
        return chi, plustimes(y)
    end
    λ = plustimes(chi.eigenvalues[1])
    # simul(..., history=true) returns nt × nout (rows = time); same shape as plustimes(y)
    adj = Float64[λ * t for t in 1:nt, _ in 1:nout]
    ypl = plustimes(y)
    ycorr = ypl .- adj
    return chi, ycorr
end

function flowshop_simu(
    s::MPSysLin,
    nm::AbstractVector{<:Integer},
    np::AbstractVector{<:Integer},
    u::AbstractMatrix{<:Real}
)
    flowshop_simu(s, nm, np, map(MP, u))
end
