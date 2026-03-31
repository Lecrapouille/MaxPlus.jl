# ==============================================================================
# This file contains workarounds and display customizations for Max-Plus types.
# Many historical Julia bugs have been fixed; this file is cleaned for Julia 1.10+.
# ==============================================================================

# ==============================================================================
# Sparse * Dense multiplication workaround
# Convert to dense to avoid dispatch issues with tropical types
Base.:(*)(A::Array{Tropical{T}}, S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax} = A * full(S)
Base.:(*)(S::SparseMatrixCSC{Tropical{T}}, A::Array{Tropical{T}}) where {T<:MinOrMax} = full(S) * A

# ==============================================================================
# Fix A^0 to return proper identity matrix with one() elements
# Julia's default creates ill-formed identity mixing zero() and true
@inline Base.literal_pow(::typeof(^), A::Matrix{Tropical{T}}, ::Val{0}) where {T<:MinOrMax} =
    eye(Tropical{T}, size(A,1), size(A,2))

# ==============================================================================
# Custom sparse matrix display (compact Julia 1.5 style)
# We prefer this over the Julia 1.6+ style which displays sparse as dense with dots

function fallback_show(io::IOContext, S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax}
    nnz(S) == 0 && return

    ioc = IOContext(io, :compact => true)
    function _format_line(r, col, padr, padc)
        print(ioc, "\n  [", rpad(rowvals(S)[r], padr), ", ", lpad(col, padc), "]  =  ")
        if isassigned(nonzeros(S), Int(r))
            show(ioc, nonzeros(S)[r])
        else
            print(ioc, Base.undef_ref_str)
        end
    end

    function _get_cols(from, to)
        idx = eltype(SparseArrays.getcolptr(S))[]
        c = searchsortedlast(SparseArrays.getcolptr(S), from)
        for i = from:to
            while i == SparseArrays.getcolptr(S)[c+1]
                c += 1
            end
            push!(idx, c)
        end
        idx
    end

    rows = displaysize(io)[1] - 4
    if !get(io, :limit, false) || rows >= nnz(S)
        cols = _get_cols(1, nnz(S))
        padr, padc = ndigits.((maximum(rowvals(S)[1:nnz(S)]), cols[end]))
        _format_line.(1:nnz(S), cols, padr, padc)
    else
        if rows <= 2
            return
        end
        s1, e1 = 1, div(rows - 1, 2)
        s2, e2 = nnz(S) - (rows - 1 - e1) + 1, nnz(S)
        cols1, cols2 = _get_cols(s1, e1), _get_cols(s2, e2)
        padr = ndigits(max(maximum(rowvals(S)[s1:e1]), maximum(rowvals(S)[s2:e2])))
        padc = ndigits(cols2[end])
        _format_line.(s1:e1, cols1, padr, padc)
        _format_line.(s2:e2, cols2, padr, padc)
    end
    return
end

function Base.show(io::IO, ::MIME"text/plain", S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax}
    xnnz = nnz(S)
    m, n = size(S)
    print(io, m, "×", n, " ", name(Tropical{T}), "sparse matrix with ", xnnz, " stored ",
              xnnz == 1 ? "entry" : "entries")
    if xnnz != 0
        print(io, ":")
        fallback_show(IOContext(io, :typeinfo => eltype(S)), S)
    end
end

Base.show(io::IO, ::Type{MIME{Symbol("text/plain")}}, S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax} =
    Base.show(io, MIME{Symbol("text/plain")}(), S)

function Base.show(io::IO, S::SparseMatrixCSC{Tropical{T}}) where {T<:MinOrMax}
    show(io, MIME{Symbol("text/plain")}(), S)
end

# ==============================================================================
mpshow(io::IO, S::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax} = show(io, S)
