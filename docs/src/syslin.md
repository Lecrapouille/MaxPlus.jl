# Toolbox (max,+) Linear Systems

## (max,+) Linear System constructor

```@docs
MaxPlus.MPSysLin
```

## (max,+) Linear System operations

```@docs
MaxPlus.mpexplicit(S::MPSysLin)
```

```@docs
MaxPlus.mpsimul(S::MPSysLin, u::MPAbstractVecOrMat, history::Bool)
```

```@docs
Base.:(+)(x::MPSysLin, y::MPSysLin)
```

```@docs
Base.:(*)(y::MPSysLin, x::MPSysLin)
```

```@docs
Base.:(|)(x::MPSysLin, y::MPSysLin)
```

```@docs
Base.:vcat(x::MPSysLin, y::MPSysLin)
```

```@docs
Base.:hcat(x::MPSysLin, y::MPSysLin)
```

```@docs
Base.:(/)(x::MPSysLin, y::MPSysLin)
```

```@docs
Base.:(==)(x::MPSysLin, y::MPSysLin)
```

```@docs
MaxPlus.mpshow(io::IO, S::MPSysLin)
```

```@docs
Base.show(io::IO, S::MPSysLin)
Base.show(io::IO, ::MIME"text/plain", S::MPSysLin)
```

```@docs
MaxPlus.LaTeX(S::MPSysLin)
MaxPlus.LaTeX(io::IO, S::MPSysLin)
```
