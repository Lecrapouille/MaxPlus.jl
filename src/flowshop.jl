################################################################################
# Flowshop (max,+)-linear systems — porting from flowshop.sci, shift.sci,
# flowshop_simu.sci (ScicosLab).
################################################################################

# ==============================================================================
# flowshop.sci — E[machine, piece]: (max-plus) duration; ε = no task.
# (docstring: docstrings/flowshop.jl)
function flowshop(E::AbstractMatrix{<:MP})
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
# (docstring: docstrings/flowshop.jl)
function flowshop_graph(
    E::AbstractMatrix{<:MP},
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
    Tacc = Dict{Tuple{Int,Int}, MP{Float64}}()
    Nacc = Dict{Tuple{Int,Int}, MP{Float64}}()
    # Scilab T(r,c)=v ; ici on construit R = T' donc R[c,r]=v
    function tset!(r::Int, c::Int, v::MP{Float64})
        Tacc[(c, r)] = v
    end
    function nset!(r::Int, c::Int, v::MP{Float64})
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

    function mat_from_acc(acc::Dict{Tuple{Int,Int}, MP{Float64}})
        I = Int[]
        J = Int[]
        V = MP{Float64}[]
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
# Pure-Julia exporter (no external dependency): write the flowshop description
# as a `.flowshop` file (the format understood by the TimedPetriNetEditor
# companion project). `G` is simply the path of the produced file, so MaxPlus.jl
# never has to call any C++ binding. (docstring: docstrings/flowshop.jl)
function save_flowshop(
    E::AbstractMatrix{<:MP},
    m::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    path::AbstractString;
    machine_names::AbstractVector{<:AbstractString} = ["M$(j)" for j in 1:size(E, 1)],
    piece_names::AbstractVector{<:AbstractString} = ["P$(i)" for i in 1:size(E, 2)],
)
    (nmach, npiece) = size(E)
    length(m) == nmach || error("save_flowshop: length(m) must equal size(E,1)")
    length(p) == npiece || error("save_flowshop: length(p) must equal size(E,2)")
    length(machine_names) == nmach || error("save_flowshop: length(machine_names) must equal size(E,1)")
    length(piece_names) == npiece || error("save_flowshop: length(piece_names) must equal size(E,2)")
    reserved = ("npieces", "nmachines", "nm", "np", "pieces")
    for name in machine_names
        (name in reserved) && error("save_flowshop: machine name '$name' collides with a reserved keyword")
        occursin(':', name) && error("save_flowshop: machine name '$name' must not contain ':'")
    end
    open(path, "w") do io
        println(io, "# Flowshop generated by MaxPlus.save_flowshop")
        println(io, "npieces: ", npiece)
        println(io, "nmachines: ", nmach)
        println(io)
        println(io, "# Tokens per machine class (nm) - horizontal cycle")
        println(io, "nm: ", join(Int.(round.(m)), " "))
        println(io)
        println(io, "# Tokens per piece class (np) - vertical cycle")
        println(io, "np: ", join(Int.(round.(p)), " "))
        println(io)
        println(io, "# Piece names (columns)")
        println(io, "pieces: ", join(piece_names, " "))
        println(io)
        println(io, "# Processing times: machine_name: t1 t2 ... (nan = no task)")
        for j in 1:nmach
            vals = String[E[j, i] == mp0 ? "nan" : string(plustimes(E[j, i])) for i in 1:npiece]
            println(io, machine_names[j], ": ", join(vals, " "))
        end
    end
    return path
end

function save_flowshop(
    E::AbstractMatrix{<:Real},
    m::AbstractVector{<:Real},
    p::AbstractVector{<:Real},
    path::AbstractString;
    kwargs...
)
    save_flowshop(map(MP, E), m, p, path; kwargs...)
end

# ==============================================================================
# flowshop_simu.sci — u: each column is a time step (as in Scilab);
# rows of u are the inputs (size = columns of B). (docstring: docstrings/flowshop.jl)
function flowshop_simu(
    s::MPSysLin,
    nm::AbstractVector{<:Integer},
    np::AbstractVector{<:Integer},
    u::AbstractMatrix{<:MP}
)
    nmach = length(nm)
    npiece = length(np)
    nt = size(u, 2)
    nt >= 1 || error("flowshop_simu: u must have at least one column (time)")

    # The machine controller
    fbm = mpshift(nm[1], 0)
    for i in 2:nmach
        fbm = fbm | mpshift(nm[i], 0)
    end
    # The piece controller
    fbp = mpshift(np[1], 0)
    for i in 2:npiece
        fbp = fbp | mpshift(np[i], 0)
    end
    # Complete feedback system
    sb = s / (fbp | fbm)
    # Reducing a system and putting it in explicit form
    sbs = explicit(sb)
    # Spectral analysis
    Sa = SparseArrays.sparse(sbs.A)
    chi = howard(Sa)
    # Simulating the system
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
