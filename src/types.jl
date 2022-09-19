# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# A portage of the ScicosLab Max-Plus toolbox http://www.scicoslab.org/
# License: public domain
#
# Note: the documentation of functions for the REPL are placed in docstrings.jl
# ==============================================================================
# This file defines types for (max,+) and (min,+) type and some of their
# containers.
# ==============================================================================

using SparseArrays

export Tropical
export MaxPlus, SpmMaxPlus, SpvMaxPlus, ArrMaxPlus, VecMaxPlus
export MinPlus, SpmMinPlus, SpvMinPlus, ArrMinPlus, VecMinPlus

# ==============================================================================
# Dummy templates to make the distinction between Min-Plus and Max-Plus numbers
struct Min end
struct Max end

# ==============================================================================
# Base class for Max-Plus and Min-Plus structures
struct Tropical{T <: Union{Min, Max}} <: Number
    v::Float64
end

# ==============================================================================
# Container of Tropical structure
const TropicalSparseMatrix{T, U} = SparseMatrixCSC{Tropical{T}, U}
const TropicalSparseVector{T, U} = SparseVector{Tropical{T}, U}
const TropicalArray{T, U} = Array{Tropical{T}, U}
const TropicalVector{T} = Vector{Tropical{T}}
const AbstractTropicalVecOrMat{T} = AbstractVecOrMat{Tropical{T}}

# ==============================================================================
# Max-Plus short structure names
const MaxPlus = Tropical{Max}
const SpmMaxPlus{U} = TropicalSparseMatrix{Max, U}
const SpvMaxPlus{U} = TropicalSparseVector{Max, U}
const ArrMaxPlus = TropicalArray{Max}
const VecMaxPlus = TropicalVector{Max}

# ==============================================================================
# Min-Plus short structure names
const MinPlus = Tropical{Min}
const SpmMinPlus{U} = TropicalSparseMatrix{Min, U}
const SpvMinPlus{U} = TropicalSparseVector{Min, U}
const ArrMinPlus = TropicalArray{Min}
const VecMinPlus = TropicalVector{Min}
