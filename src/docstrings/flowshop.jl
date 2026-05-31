# ==============================================================================
# Docstrings for the Flowshop (max,+)-linear systems.
# Ported from ScicosLab: shift.sci, flowshop.sci, flowshop_graph.sci,
# flowshop_simu.sci. The function bodies live in src/flowshop.jl.
# ==============================================================================

################################################################################
###
### Flowshop
###
################################################################################

# ==============================================================================
"""
    mpshift(n::Integer, t::Real)

Port of ScicosLab `shift.sci`. Build the [`MPSysLin`](@ref) of an `n`-event
shift (delay) with a date offset `t`. Used to model the holding of `n`
machines/pallets in a flowshop feedback. `n` must be `>= 1`.
"""
mpshift(n::Integer, t::Real)

# ==============================================================================
"""
    flowshop(E::AbstractMatrix)

Port of ScicosLab `flowshop.sci`. From the processing-time matrix `E`
(rows = machine classes, columns = part classes, entry = task duration, `mp0`
= no task), build the (max,+) linear system [`MPSysLin`](@ref) of the cyclic
flowshop. `E` may be given as a (max,+) matrix or as a classic real matrix
(it is then converted with `mp0` for missing tasks).
"""
flowshop(E::AbstractMatrix{<:MP})

# ==============================================================================
"""
    T, N = flowshop_graph(E, m, p)

Port of ScicosLab `flowshop_graph.sci` (matrices only; no graph drawing nor
`show_graph`). From the processing-time matrix `E`, the vector `m` of machine
counts per class and the vector `p` of pallet counts per part class, build:

- `T`: sparse (max,+) matrix of dates (durations on the arcs);
- `N`: sparse (max,+) matrix of counters (token weights on the arcs).

`T` and `N` share the same sparsity pattern and are returned already
transposed, as in ScicosLab (`T = sparse(T'); N = sparse(N')`). Feed them to
[`semihoward`](@ref) / [`mpeigen`](@ref) to get the cycle time.

!!! note "No graph is returned"
    Unlike the ScicosLab `flowshop_graph`, this port does **not** draw nor
    return the timed event graph. Export with [`save_flowshop`](@ref) and call
    `show_cr_graph` in [TimedPetriNetEditor.jl](https://github.com/Lecrapouille/TimedPetriNetEditor.jl).
"""
flowshop_graph(E::AbstractMatrix{<:MP}, m::AbstractVector{<:Real}, p::AbstractVector{<:Real})

# ==============================================================================
"""
    chi, y = flowshop_simu(s, nm, np, u)

Port of ScicosLab `flowshop_simu.sci`. Close the flowshop system `s`
([`MPSysLin`](@ref) built by [`flowshop`](@ref)) with the shift feedbacks of
`nm` machines and `np` pallets, make it [`explicit`](@ref) and simulate it on
the input `u` (rows = inputs, columns = time steps).

Returns:
- `chi`: the [`HowardResult`](@ref) of [`howard`](@ref) on the closed-loop `A`;
- `y`: the classic (`Float64`) outputs `nt × (nmach + npiece)` with the cyclic
  drift `λ * t` already subtracted (`λ = chi.eigenvalues[1]`).
"""
flowshop_simu(s::MPSysLin, nm::AbstractVector{<:Integer}, np::AbstractVector{<:Integer}, u::AbstractMatrix{<:MP})

# ==============================================================================
"""
    save_flowshop(E, m, p, path; machine_names, piece_names)

Export the flowshop description to a `.flowshop` text file at `path` and return
that path (so it can be used as the graph handle `G`). The format is the one
read by [TimedPetriNetEditor.jl](https://github.com/Lecrapouille/TimedPetriNetEditor.jl)
(C++ backend: [TimedPetriNetEditor](https://github.com/Lecrapouille/TimedPetriNetEditor)):
`npieces`, `nmachines`, `nm`, `np`, `pieces`, then one line per machine with the
processing times (`mp0` / no task is written as `nan`).

This is **pure Julia** with no external dependency: `MaxPlus.jl` only *writes*
the description, it never calls the C++ binding. Load the produced file in
`TimedPetriNetEditor.jl` (`show_cr_graph(path)`) to build the network and show
its critical cycle.

# Arguments
- `E`: processing-time matrix (rows = machines, columns = parts, `mp0` = no task).
- `m`: machine counts per class (`length(m) == size(E,1)`).
- `p`: pallet counts per part class (`length(p) == size(E,2)`).
- `machine_names` / `piece_names`: optional labels (default `M1..`, `P1..`).
"""
save_flowshop(E::AbstractMatrix{<:MP}, m::AbstractVector{<:Real}, p::AbstractVector{<:Real}, path::AbstractString)
