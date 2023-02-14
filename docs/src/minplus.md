# API: (min,+) Algebra

## (min,+) constructor

### Scalar Constructors

```@docs
MaxPlus.MI(::Float64)
```

```@docs
MaxPlus.MI(::Bool)
```

### Dense Matrix and Dense Vector Constructors

```@docs
MaxPlus.MI(::Array)
```

### Sparse Matrix Constructors

```@docs
MaxPlus.MI(::SparseMatrixCSC)
```

```@docs
MaxPlus.MI(I::AbstractVector, J::AbstractVector, V::AbstractVector)
```

### Sparse Vector Constructors

```@docs
MaxPlus.MI(::SparseVector)
```

```@docs
MaxPlus.MI(I::AbstractVector, V::AbstractVector)
```

### Unit Range Constructors

```@docs
MaxPlus.MI(::UnitRange)
```

```@docs
MaxPlus.MI(::StepRangeLen)
```

## Overloaded Algebraic Operators

```@docs
Base.zero(::MI)
Base.zero(::Type{MI})
```

```@docs
Base.one(::MI)
Base.one(::Type{MI})
```

```@docs
Base.:(+)(::MI, ::MI)
```

```@docs
Base.:(*)(::MI, ::MI)
```

```@docs
Base.:(^)(::MI, ::Number)
```

```@docs
Base.:(/)(::MI, ::MI)
```

```@docs
Base.:(\)(::MI, ::MI)
```

```@docs
Base.:(-)(::MI, ::MI)
```

```@docs
Base.:min(x::MI, y::MI)
```

## (min,+) algebra to classic algebra conversion

```@docs
MaxPlus.plustimes(::MI)
```

## (min,+) Constants

```@docs
MaxPlus.mi0
```

```@docs
MaxPlus.mi1
```

```@docs
MaxPlus.mie
```

```@docs
MaxPlus.mitop
```

```@docs
MaxPlus.miI
```

## (min,+) Dense matrices constructions

```@docs
Base.ones(MI, m::Int64, n::Int64)
Base.ones(MI, n::Int64)
```

```@docs
Base.zeros(MI, m::Int64, n::Int64)
Base.zeros(MI, n::Int64)
```

```@docs
MaxPlus.eye(MI, m::Int64, n::Int64)
MaxPlus.eye(MI, n::Int64)
```

## (min,+) Sparse matrices constructions

```@docs
SparseArrays.spzeros(MI, m::Int64, n::Int64)
SparseArrays.spzeros(MI, n::Int64)
```

```@docs
MaxPlus.speye(MI, m::Int64, n::Int64)
MaxPlus.speye(MI, n::Int64)
```

## (min,+) matrices Conversion

```@docs
MaxPlus.plustimes(A::Array{MI})
```

```@docs
MaxPlus.plustimes(S::SparseMatrixCSC{MI})
```

```@docs
MaxPlus.full(S::SparseMatrixCSC{MI})
```

```@docs
MaxPlus.dense(S::SparseMatrixCSC{MI})
```

## (min,+) Matrix operations

```@docs
Base.:(\)(::AbstractMatrix{MI}, ::AbstractMatrix{MI})
Base.:(\)(::AbstractMatrix{MI}, ::MI)
Base.:(\)(::MI, ::AbstractMatrix{MI})
```

```@docs
Base.:(/)(::AbstractMatrix{MI}, ::AbstractMatrix{MI})
Base.:(/)(::AbstractMatrix{MI}, ::MI)
Base.:(/)(::MI, ::AbstractMatrix{MI})
```

```@docs
Base.inv(::Matrix{MI})
```

```@docs
MaxPlus.star(A::Array{MI})
```

```@docs
MaxPlus.plus(A::Array{MI})
```

```@docs
MaxPlus.astarb(A::Array{MI}, b::Array{MI})
```

```@docs
MaxPlus.MIeigen
```

```@docs
MaxPlus.howard(S::SparseMatrixCSC{MI})
```

```@docs
MaxPlus.tr(A::Array{MI})
```

```@docs
MaxPlus.norm(A::Array{MI})
```

## Display control of (min,+) scalar and Matrices

```@docs
MaxPlus.set_tropical_display
MaxPlus.LaTeX(io::IO, A::Matrix{MI})

#Base.show(::IO, ::Matrix{MI})
#Base.show(::IO, ::SparseMatrixCSC{MI})

Base.show(::IO, ::MIME"text/plain", A::Matrix{MI})
Base.show(::IO, ::MIME"text/latex", A::Matrix{MI})
```

## Index

```@index
```
