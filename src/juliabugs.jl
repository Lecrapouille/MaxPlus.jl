# ==============================================================================
# This file fixes a some issues found in Julia.
# ==============================================================================

# ==============================================================================
# https://github.com/JuliaLang/julia/issues/33332
# Fix bug in Sparse Matrix in Julia. They used vals1[j1] != 0 instead of using
# iszero(vals1[j1]) idem for vals2[j2]

function Base.:(==)(A1::SparseMatrixCSC, A2::SparseMatrixCSC)
    size(A1) != size(A2) && return false
    vals1, vals2 = nonzeros(A1), nonzeros(A2)
    rows1, rows2 = rowvals(A1), rowvals(A2)
    m, n = size(A1)
    @inbounds for i = 1:n
        nz1,nz2 = nzrange(A1,i), nzrange(A2,i)
        j1,j2 = first(nz1), first(nz2)
        # step through the rows of both matrices at once:
        while j1 <= last(nz1) && j2 <= last(nz2)
            r1,r2 = rows1[j1], rows2[j2]
            if r1==r2
                vals1[j1]!=vals2[j2] && return false
                j1+=1
                j2+=1
            else
                if r1<r2
                    !iszero(vals1[j1]) && return false
                    j1+=1
                else
                    !iszero(vals2[j2]) && return false
                    j2+=1
                end
            end
        end
        # finish off any left-overs:
        for j = j1:last(nz1)
            !iszero(vals1[j]) && return false
        end
        for j = j2:last(nz2)
            !iszero(vals2[j]) && return false
        end
    end
    return true
end

# ==============================================================================
# https://github.com/JuliaLang/julia/issues/33332 fixed in Julia 1.3
# Fix bug of sparse[i,j] = mp0 confusing with sparse[i,j] = 0

getcolptr(S::SparseMatrixCSC) = getfield(S, :colptr)

function _insert!(v::Vector, pos::Integer, item, nz::Integer)
    if nz > length(v)
        insert!(v, pos, item)
    else # nz < length(v)
        Base.unsafe_copyto!(v, pos+1, v, pos, nz - pos)
        v[pos] = item
        v
    end
end

function _setindex_scalar!(A::SparseMatrixCSC{MP{Tv},Ti}, _v, _i::Integer, _j::Integer) where {Tv,Ti<:Integer}
    v = convert(MP{Tv}, _v)
    i = convert(Ti, _i)
    j = convert(Ti, _j)
    if !((1 <= i <= size(A, 1)) & (1 <= j <= size(A, 2)))
        throw(BoundsError(A, (i,j)))
    end
    coljfirstk = Int(getcolptr(A)[j])
    coljlastk = Int(getcolptr(A)[j+1] - 1)
    searchk = searchsortedfirst(rowvals(A), i, coljfirstk, coljlastk, Base.Order.Forward)
    if searchk <= coljlastk && rowvals(A)[searchk] == i
        nonzeros(A)[searchk] = v
        return A
    end
    if !iszero(v)
        nz = getcolptr(A)[size(A, 2)+1]
        !isbitstype(Ti) || nz < typemax(Ti) ||
            throw(ArgumentError("nnz(A) going to exceed typemax(Ti) = $(typemax(Ti))"))
        _insert!(rowvals(A), searchk, i, nz)
        _insert!(nonzeros(A), searchk, v, nz)
        @simd for m in (j + 1):(size(A, 2) + 1)
            @inbounds getcolptr(A)[m] += Ti(1)
        end
    end
    return A
end

@inline Base.:setindex!(A::SparseMatrixCSC{MP{Tv},Ti}, _v, _i::Integer, _j::Integer) where {Tv,Ti<:Integer} =
    _setindex_scalar!(A, _v, _i, _j)

# Used by SimpleWeightedDiGraph
function SparseMatrixCSC{MP{Tv},Ti}(M::StridedMatrix) where {Tv,Ti}
    nz = count(!iszero, M)
    colptr = zeros(Ti, size(M, 2) + 1)
    nzval = Vector{MP{Tv}}(undef, nz)
    rowval = Vector{Ti}(undef, nz)
    colptr[1] = 1
    cnt = 1
    @inbounds for j in 1:size(M, 2)
        for i in 1:size(M, 1)
            v = M[i, j]
            if !iszero(v)
                rowval[cnt] = i
                nzval[cnt] = v
                cnt += 1
            end
        end
        colptr[j+1] = cnt
    end
    return SparseMatrixCSC(size(M, 1), size(M, 2), colptr, rowval, nzval)
end

# ==============================================================================
# Bug https://github.com/JuliaLang/julia/issues/33036

"""
    mpI

Fix an algebra conception in Julia official LinearAlgebra (uniformscaling.jl)
`I` representes an identity matrix of any size is defined with a booleen instead
of the function `one()`. As consequence, in Julia 0.4 the `eye(T,m,n)` could
created a max-plus identity matrix. Since Julia 0.7 `eye()` has been deprecated
and replaced by the buggy `Matrix{T}(I, m, n)`. This function uses `zero()` but
not `one()` and as consequence the max-plus identity matrix is not well formed.

This const allows to be more algebra compliant by calling `one()` and fixing the
fucntion `Matrix{T}(I, m, n)`.
"""
const global mpI = UniformScaling(one(MP{Float64}).λ)

# ==============================================================================
# Because Julia will create the ill-formed identity matrix mixing zero() and true
# instead of zero() and one()

@inline Base.literal_pow(::typeof(^), A::ArrMP{T}, ::Val{0}) where T =
    mpeye(T, size(A,1), size(A,2))
