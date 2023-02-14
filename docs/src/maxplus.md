# API: (max,+) Algebra

## (max,+) constructor

### Scalar Constructors

```@docs
MaxPlus.MP(::Float64)
```

```@docs
MaxPlus.MP(::Bool)
```

### Dense Matrix and Dense Vector Constructors

```@docs
MaxPlus.MP(::Array)
```

### Sparse Matrix Constructors

```@docs
MaxPlus.MP(::SparseMatrixCSC)
```

```@docs
MaxPlus.MP(I::AbstractVector, J::AbstractVector, V::AbstractVector)
```

### Sparse Vector Constructors

```@docs
MaxPlus.MP(::SparseVector)
```

```@docs
MaxPlus.MP(I::AbstractVector, V::AbstractVector)
```

### Unit Range Constructors

```@docs
MaxPlus.MP(::UnitRange)
```

```@docs
MaxPlus.MP(::StepRangeLen)
```

## Overloaded Algebraic Operators

```@docs
Base.zero(::MP)
Base.zero(::Type{MP})
```

```@docs
Base.one(::MP)
Base.one(::Type{MP})
```

```@docs
Base.:(+)(::MP, ::MP)
```

```@docs
Base.:(*)(::MP, ::MP)
```

```@docs
Base.:(^)(::MP, ::Number)
```

```@docs
Base.:(/)(::MP, ::MP)
```

```@docs
Base.:(\)(::MP, ::MP)
```

```@docs
Base.:(-)(::MP, ::MP)
```

```@docs
Base.:min(x::MP, y::MP)
```

## (max,+) algebra to classic algebra conversion

```@docs
MaxPlus.plustimes(::MP)
```

```@docs
MaxPlus.star(::MP)
```

```@docs
MaxPlus.plus(::MP)
```

## (max,+) Constants

```@docs
MaxPlus.mp0
```

```@docs
MaxPlus.ε
```

```@docs
MaxPlus.mp1
```

```@docs
MaxPlus.mpe
```

```@docs
MaxPlus.mptop
```

```@docs
MaxPlus.mpI
```

## (max,+) Dense matrices constructions

```@docs
Base.ones(MP, m::Int64, n::Int64)
Base.ones(MP, n::Int64)
Base.ones(MP, A::Array)
```

```@docs
Base.zeros(MP, m::Int64, n::Int64)
Base.zeros(MP, n::Int64)
Base.zeros(MP, A::Array)
```

```@docs
MaxPlus.eye(MP, m::Int64, n::Int64)
MaxPlus.eye(MP, n::Int64)
MaxPlus.eye(MP, A::Array)
```

## (max,+) Sparse matrices constructions

```@docs
SparseArrays.spzeros(MP, m::Int64, n::Int64)
SparseArrays.spzeros(MP, n::Int64)
SparseArrays.spzeros(MP, A::Array)
```

```@docs
MaxPlus.speye(MP, m::Int64, n::Int64)
MaxPlus.speye(MP, n::Int64)
```

## (max,+) matrices Conversion

```@docs
MaxPlus.plustimes(A::Array{MP})
```

```@docs
MaxPlus.plustimes(S::SparseMatrixCSC{MP})
```

```@docs
MaxPlus.full(S::SparseMatrixCSC{MP})
```

```@docs
MaxPlus.dense(S::SparseMatrixCSC{MP})
```

## (max,+) Matrix operations

```@docs
Base.:(\)(::AbstractMatrix{MP}, ::AbstractMatrix{MP})
Base.:(\)(::AbstractMatrix{MP}, ::MP)
Base.:(\)(::MP, ::AbstractMatrix{MP})
```

```@docs
Base.:(/)(::AbstractMatrix{MP}, ::AbstractMatrix{MP})
Base.:(/)(::AbstractMatrix{MP}, ::MP)
Base.:(/)(::MP, ::AbstractMatrix{MP})
```

```@docs
Base.inv(::Matrix{MP})
```

```@docs
MaxPlus.star(A::Array{MP})
```

```@docs
MaxPlus.plus(A::Array{MP})
```

```@docs
MaxPlus.astarb(A::Array{MP}, b::Array{MP})
```

```@docs
MaxPlus.mpeigen
```

```@docs
MaxPlus.howard(S::SparseMatrixCSC{MP})
```

```@docs
MaxPlus.tr(A::Array{MP})
```

```@docs
MaxPlus.norm(A::Array{MP})
```

```@docs
MaxPlus.spget(S::SparseMatrixCSC{MP})
```

```@docs
MaxPlus.sparse_map(f, S::SparseMatrixCSC{MP})
```

## Display control of (max,+) scalar and Matrices

```@docs
MaxPlus.set_tropical_display
MaxPlus.LaTeX(io::IO, A::Matrix{MP})

#Base.show(::IO, ::Matrix{MP})
#Base.show(::IO, ::SparseMatrixCSC{MP})

Base.show(::IO, ::MIME"text/plain", A::Matrix{MP})
Base.show(::IO, ::MIME"text/latex", A::Matrix{MP})
```

## Index

```@index
```
