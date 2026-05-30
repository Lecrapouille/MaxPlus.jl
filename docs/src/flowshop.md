# API: Flowshop

Functions porting the ScicosLab flowshop toolbox (`flowshop.sci`, `shift.sci`,
`flowshop_graph.sci`, `flowshop_simu.sci`) to build and simulate cyclic
flowshop `(max,+)` linear systems.

A worked example (with the same processing-time matrix as the ScicosLab demo)
is available in the [tutorial notebooks](tutorial.md) `flowshop-fr.ipynb` and
`flowshop-en.ipynb`.

!!! note "No graph drawing"
    [`flowshop_graph`](@ref) returns only the `(T, N)` sparse `(max,+)`
    matrices; it does **not** draw nor return the timed event graph. To build
    the network from the same matrix and display its critical cycles, use the
    companion project
    [TimedPetriNetEditor](https://github.com/Lecrapouille/TimedPetriNetEditor).

## Flowshop system construction

```@docs
MaxPlus.flowshop(E::AbstractMatrix{<:MP})
```

```@docs
MaxPlus.mpshift(n::Integer, t::Real)
```

## Flowshop graph and spectral analysis

```@docs
MaxPlus.flowshop_graph(E::AbstractMatrix{<:MP}, m::AbstractVector{<:Real}, p::AbstractVector{<:Real})
```

## Flowshop simulation

```@docs
MaxPlus.flowshop_simu(s::MPSysLin, nm::AbstractVector{<:Integer}, np::AbstractVector{<:Integer}, u::AbstractMatrix{<:MP})
```
