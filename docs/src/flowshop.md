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
    the network from the same matrix and display its critical cycles, use
    [`MaxPlus.save_flowshop`](@ref) then
    [`TimedPetriNetEditor.show_cr_graph`](https://github.com/Lecrapouille/TimedPetriNetEditor.jl)
    from the companion Julia package
    [TimedPetriNetEditor.jl](https://github.com/Lecrapouille/TimedPetriNetEditor.jl).

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

## Flowshop graph export (for TimedPetriNetEditor)

```@docs
MaxPlus.save_flowshop(E::AbstractMatrix{<:MP}, m::AbstractVector{<:Real}, p::AbstractVector{<:Real}, path::AbstractString)
```

## Flowshop simulation

```@docs
MaxPlus.flowshop_simu(s::MPSysLin, nm::AbstractVector{<:Integer}, np::AbstractVector{<:Integer}, u::AbstractMatrix{<:MP})
```
