# (max,+) Algebra Toolbox Julia API

## Main (max,+) structure

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
MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)
```

```@docs
MP(I::AbstractVector, V::AbstractVector)
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
Base.:(\)(::Matrix{MP}, ::Matrix{MP})
```

```@docs
Base.:(^)(::MP, ::Number)
```

## (max,+) constants

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

## (max,+) Zeros

```@docs
zero(::MP)
zero(::Type{MP})
```

```@docs
zero
```

```@docs
spzeros
```

```@docs
zeros
```

## (max,+) Ones

```@docs
one(::MP)
one(::Type{MP})
```

```@docs
one
```

```@docs
ones
```

## (max,+) Identity Matrix

```@docs
eye
```

```@docs
mpI
```

## (max,+) Convertion

```@docs
plustimes
```

```@docs
full
```

```@docs
dense
```

## Power

```@docs
Base.inv(::Matrix{MP})
```

## (max,+) Matrix operations

```@docs
tr
```

```@docs
norm
```

```@docs
plus
```

```@docs
star
```

```@docs
astarb
```

```@docs
howard
```

```@docs
semihoward
```

## Display control of (max,+) scalar and Matrices

```@docs
set_tropical_display
LaTeX
mpshow
Base.show(::IO, ::MP)
Base.show(::IO, ::Matrix{MP})
Base.show(::IO, ::SparseMatrixCSC{MP})

Base.show(::IO, ::MIME"text/plain", ::MP)
Base.show(::IO, ::MIME"text/plain", ::Matrix{MP})
Base.show(::IO, ::MIME"text/plain", ::SparseMatrixCSC{MP})

Base.show(::IO, ::MIME"text/latex", ::MP)
Base.show(::IO, ::MIME"text/latex", ::Matrix{MP})
Base.show(::IO, ::MIME"text/latex", ::SparseMatrixCSC{MP})
```

## Index

```@index
```
