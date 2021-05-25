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

function _setindex_scalar!(A::SparseMatrixCSC{MP,Ti}, _v, _i::Integer, _j::Integer) where {Tv,Ti<:Integer}
    v = convert(MP, _v)
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

@inline Base.:setindex!(A::SparseMatrixCSC{MP,Ti}, _v, _i::Integer, _j::Integer) where {Tv,Ti<:Integer} =
    _setindex_scalar!(A, _v, _i, _j)

# Used by SimpleWeightedDiGraph
function SparseMatrixCSC{MP,Ti}(M::StridedMatrix) where {Tv,Ti}
    nz = count(!iszero, M)
    colptr = zeros(Ti, size(M, 2) + 1)
    nzval = Vector{MP}(undef, nz)
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
# Since Julia 1.4.x the matrix product Max-Plus sparse * full or full * sparse
# is no longer working (while sparse * sparse keep working)
Base.:(*)(A::ArrMP, S::SpaMP) = A * full(S)
Base.:(*)(S::SpaMP, A::ArrMP) = full(S) * A

# ==============================================================================
# Since Julia 1.6.x sparse matrices are displayed like dense matrices but with
# dots instead of zero elements. We prefer the older Julia 1.5 display because
# more compact, so we force older display but only for sparse Max-Plus only and
# for sparse matrices of other type because our purpose is not to change the
# Julia behavior but we cherry pick good behavior for our toolbox.

function Base.show(io::IO, ::MIME"text/plain", S::SpaMP)
    xnnz = nnz(S)
    m, n = size(S)
    print(io, m, "×", n, " Max-Plus sparse matrix with ", xnnz, " stored ",
              xnnz == 1 ? "entry" : "entries")
    if xnnz != 0
        print(io, ":")
        show(IOContext(io, :typeinfo => eltype(S)), S)
    end
end

function Base.show(io::IO, S::SpaMP)
    Base.show(convert(IOContext, io), S::SpaMP)
end

function Base.show(io::IOContext, S::SpaMP)
    nnz(S) == 0 && return show(io, MIME("text/plain"), S)

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
        idx = eltype(getcolptr(S))[]
        c = searchsortedlast(getcolptr(S), from)
        for i = from:to
            while i == getcolptr(S)[c+1]
                c +=1
            end
            push!(idx, c)
        end
        idx
    end

    rows = displaysize(io)[1] - 4 # -4 from [Prompt, header, newline after elements, new prompt]
    if !get(io, :limit, false) || rows >= nnz(S) # Will the whole matrix fit when printed?
        cols = _get_cols(1, nnz(S))
        padr, padc = ndigits.((maximum(rowvals(S)[1:nnz(S)]), cols[end]))
        _format_line.(1:nnz(S), cols, padr, padc)
    else
        if rows <= 2
            print(io, "\n  \u22ee")
            return
        end
        s1, e1 = 1, div(rows - 1, 2) # -1 accounts for \vdots
        s2, e2 = nnz(S) - (rows - 1 - e1) + 1, nnz(S)
        cols1, cols2 = _get_cols(s1, e1), _get_cols(s2, e2)
        padr = ndigits(max(maximum(rowvals(S)[s1:e1]), maximum(rowvals(S)[s2:e2])))
        padc = ndigits(cols2[end])
        _format_line.(s1:e1, cols1, padr, padc)
        print(io, "\n  \u22ee")
        _format_line.(s2:e2, cols2, padr, padc)
    end
    return
end
