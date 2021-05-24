# Max-Plus Algebra Toolbox Julia API

## Main Max-Plus structure

### Scalar

```@docs
MaxPlus.MP
```

### Dense Matrix and Dense Vector

```@docs
MP(::Array)
```

### Sparse Matrix and Sparse Vector

```@docs
MP(::SparseVector)
```

```@docs
MP(::SparseMatrixCSC)
```

```@docs
MP(I::AbstractVector{Ti}, J::AbstractVector{Ti}, V::AbstractVector{Tv}) where {Tv,Ti<:Integer}
```

### Unit Range

```@docs
MP(::UnitRange)
MP(::StepRangeLen)
```

## Overriden operators

```@docs
Base.:(+)(::MP, ::MP)
```

```@docs
Base.:(*)(::MP, ::MP)
```

```@docs
Base.:(\)(::ArrMP, ::ArrMP)
```

```@docs
Base.:(^)(::MP, ::Number)
```

## Max-Plus constants

```@docs
mp0
```

```@docs
ε
```

```@docs
mp1
```

```@docs
mpe
```

```@docs
mptop
```

```@docs
mpI
```

## Max-Plus Zeros

```@docs
zero(::MP)
zero(::Type{MP})
```

```@docs
mpzero
```

```@docs
mpzeros
```

## Max-Plus Ones

```@docs
one(::MP)
one(::Type{MP})
```

```@docs
mpone
```

```@docs
mpones
```

## Max-Plus Identity Matrix

```@docs
mpeye
```

## Max-Plus Convertion

```@docs
plustimes
```

```@docs
full
```

```@docs
dense
```

```@docs
array
```

## Power

```@docs
Base.inv(::ArrMP)
```

## Max-Plus Matrix operations

```@docs
mptrace
```

```@docs
mpnorm
```

```@docs
mpplus
```

```@docs
mpstar
```

```@docs
mpastarb
```

```@docs
howard
```

## Display control of Max-Plus scalar and Matrices

```@docs
mp_change_display
LaTeX
mpshow
Base.show(::IO, ::MP)
Base.show(::IO, ::ArrMP)
Base.show(::IO, ::SpaMP)

Base.show(::IO, ::MIME"text/plain", ::MP)
Base.show(::IO, ::MIME"text/plain", ::ArrMP)
Base.show(::IO, ::MIME"text/plain", ::SpaMP)

Base.show(::IO, ::MIME"text/latex", ::MP)
Base.show(::IO, ::MIME"text/latex", ::ArrMP)
Base.show(::IO, ::MIME"text/latex", ::SpaMP)
```

## Index

```@index
```
