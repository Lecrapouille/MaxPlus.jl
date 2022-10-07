# ==============================================================================
# Trop{T} types

# Fake templates to make difference between Min-Plus and Max-Plus numbers
struct Min end
struct Max end
const MM = Union{Min, Max}

# Base class for Max-Plus and Min-Plus structures
struct Trop{T <: MM} <: Number
    λ::Float64
end

# Max-Plus structure
const MP = Trop{Max}

# Min-Plus structure
const MI = Trop{Min}

# Sparse matrix (classic algebra) shorter name
const SparseMatrix{T,U} = SparseMatrixCSC{T,U}

# Sparse matrix (max-plus algebra) shorter name
const SpaTrop{T,U} = SparseMatrix{Trop{T},U}
const SpaMP{U} = SparseMatrix{MP,U}
const SpaMI{U} = SparseMatrix{MI,U}

# Sparse vector (max-plus algebra) shorter name
const SpvTrop{T,U} = SparseVector{Trop{T},U}
const SpvMP{U} = SparseVector{MP,U}
const SpvMI{U} = SparseVector{MI,U}

# Dense matrix (max-plus algebra) shorter name
const ArrTrop{T,N} = Array{Trop{T},N}
const ArrMP{N} = Array{MP,N}
const ArrMI{N} = Array{MI,N}
const AbsMatMP = AbstractMatrix{MP}

# Dense vector (max-plus algebra) shorter name
const VecTrop{T} = Vector{Trop{T}}
const VecMP = Vector{MP}
const VecMI = Vector{MI}
