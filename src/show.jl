# ##############################################################################
# Display Max-Plus scalars and matrices
# ##############################################################################

# ==============================================================================
name(::Type{Tropical{Max}}) = "(max,+) "
name(::Type{Tropical{Min}}) = "(min,+) "
name(::Tropical{Max}) = "(max,+) "
name(::Tropical{Min}) = "(min,+) "

# ==============================================================================

# Memory for saving the style of display for neutral and absorbing tropical
# elements. By default the display of ScicosLab will be used.
const DISPLAY_STYLE = Ref{Int}(1)

# ==============================================================================
# Change the style of behavior of functions Base.show()
function set_tropical_display(style::Int)
    if style < 0 || style > 4
        error("style shall be 0 .. 4")
    end
    DISPLAY_STYLE[] = style
    if style == 0
        (@printf stdout "I will show -Inf, +Inf and 0.0")
    elseif style == 1
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as .")
    elseif style == 2
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as . and 0.0 as e")
    elseif style == 3
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as ε")
    elseif style == 4
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as ε and 0.0 as e")
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) scalar.
function tropshow(io::IO, x::Tropical{T}) where {T<:MinOrMax}
    style = DISPLAY_STYLE[]
    if style == 0
        show(io, x.λ)
    elseif x == zero(x)
        (style == 1 || style == 2) ? (@printf io ".") : (@printf io "ε")
    elseif x == one(x)
        (style == 1 || style == 3) ? (@printf io "0") : (@printf io "e")
    elseif x.λ == trunc(x.λ)
        (@printf io "%d" x.λ)
    else
        show(io, x.λ)
    end
end

# ==============================================================================
# Internal function to display a matrix with aligned columns (replaces PrettyTables)
function _matrix_show(io::IO, A::AbstractArray{<:Tropical})
    buf = IOBuffer()
    rows, cols = size(A, 1), size(A, 2)
    strs = Matrix{String}(undef, rows, cols)
    for i in 1:rows
        for j in 1:cols
            tropshow(buf, A[i,j])
            strs[i,j] = String(take!(buf))
        end
    end
    widths = [maximum(length(strs[i,j]) for i in 1:rows) for j in 1:cols]
    for i in 1:rows
        print(io, "  ")
        for j in 1:cols
            print(io, lpad(strs[i,j], widths[j]))
            j < cols && print(io, "   ")
        end
        println(io)
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) dense array.
function tropshow(io::IO, A::AbstractArray{Tropical{T},N}) where {T<:MinOrMax,N}
    if N == 1
        print(io, size(A,1), "-element ", name(Tropical{T}), "vector:\n")
        _matrix_show(io, reshape(A, :, 1))
    elseif N == 2
        print(io, size(A,1), '×', size(A,2), " ", name(Tropical{T}), "dense matrix:\n")
        _matrix_show(io, A)
    else
        show(io, A)
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) transposed matrix.
function tropshow(io::IO, A::LinearAlgebra.Transpose{Tropical{T}, Array{Tropical{T},N}}) where {T<:MinOrMax,N}
    if N == 1
        print(io, size(A,1), '×', size(A,2), " ", name(Tropical{T}), "transposed vector:\n")
        _matrix_show(io, A)
    elseif N == 2
        print(io, size(A,1), '×', size(A,2), " ", name(Tropical{T}), "transposed dense matrix:\n")
        _matrix_show(io, A)
    else
        show(io, A)
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) transposed vector.
function tropshow(io::IO, V::LinearAlgebra.Transpose{Tropical{T}, Vector{Tropical{T}}}) where {T<:MinOrMax}
    print(io, size(V,1), '×', size(V,2), " ", name(Tropical{T}), "transposed vector:\n")
    _matrix_show(io, V)
end

# ==============================================================================
function tropshow(io::IO, S::MPSysLin)
    (@printf io "Implicit dynamic linear Max-Plus system:\n")
    (@printf io "  x(n) = D*x(n) + A*x(n-1) + B*u(n)\n  y(n) = C*x(n)\n  x(0) = x0\n\nwith:")
    (@printf io "\nD = "); tropshow(io, S.D)
    (@printf io "\nA = "); tropshow(io, S.A)
    (@printf io "\nB = "); tropshow(io, S.B)
    (@printf io "\nC = "); tropshow(io, S.C)
    (@printf io "\nx0 = "); tropshow(io, S.x0)
end

# ==============================================================================
# Julia 1.12+: `MIME"text/plain"` in *value* position is the DataType
# `MIME{Symbol("text/plain")}` (singleton type object), not a `MIME` instance.
# repr/IJulia then call `show(io, ::Type{MIME{Symbol("text/plain")}}, x)`.
# Normalize by forwarding to the instance-based 3-arg show.
function _show_mime_plain(io::IO, x::Tropical{T}) where {T<:MinOrMax}
    print(io, name(x))
    tropshow(io, x)
end

Base.show(io::IO, x::Tropical) = tropshow(io, x)
Base.show(io::IO, S::MPSysLin) = tropshow(io, S)
Base.show(io::IO, ::MIME"text/plain", x::Tropical{T}) where {T<:MinOrMax} = _show_mime_plain(io, x)
Base.show(io::IO, ::Type{MIME{Symbol("text/plain")}}, x::Tropical{T}) where {T<:MinOrMax} =
    Base.show(io, MIME{Symbol("text/plain")}(), x)

# ==============================================================================
# Called by the REPL through the display() method. These functions fix
# misaliged columns made by the default show() Julia. We use the package
# PrettyTables.

Base.show(io::IO, ::MIME"text/plain", A::VecOrMat{Tropical{T}}) where {T<:MinOrMax} = tropshow(io, A)
Base.show(io::IO, ::Type{MIME{Symbol("text/plain")}}, A::VecOrMat{Tropical{T}}) where {T<:MinOrMax} =
    Base.show(io, MIME{Symbol("text/plain")}(), A)

Base.show(io::IO, ::MIME"text/plain", S::MPSysLin) = tropshow(io, S)
Base.show(io::IO, ::Type{MIME{Symbol("text/plain")}}, S::MPSysLin) =
    Base.show(io, MIME{Symbol("text/plain")}(), S)

Base.show(io::IO, ::MIME"text/plain", A::LinearAlgebra.Transpose{Tropical{T}, AbstractArray{Tropical{T},N}}) where {T<:MinOrMax,N} = tropshow(io, A)
Base.show(io::IO, ::Type{MIME{Symbol("text/plain")}}, A::LinearAlgebra.Transpose{Tropical{T}, AbstractArray{Tropical{T},N}}) where {T<:MinOrMax,N} =
    Base.show(io, MIME{Symbol("text/plain")}(), A)

# ==============================================================================
#
function LaTeX(A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax}
    style = DISPLAY_STYLE[]
    s = "\\left[\n\\begin{array}{*{20}c}\n"
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == zero(Tropical{T})
                if style == 0
                    if A[i,j] < 0
                        s = s * "-\\infty"
                    else
                        s = s * "+\\infty"
                    end
                elseif style == 3 || style == 4
                    s = s * "\\varepsilon"
                else
                    s = s * "."
                end
            elseif A[i,j] == one(Tropical{T})
                if style == 2 || style == 4
                    s = s * "e"
                else
                    s = s * string(Int64(A[i,j].λ))
                end
            elseif A[i,j].λ == trunc(A[i,j].λ)
                s = s * string(Int64(A[i,j].λ))
            else
                s = s * string(A[i,j].λ)
            end
            if j < size(A, 2)
                s = s * " & "
            end
        end
        s = s * " \\\\\n"
    end
    return s * "\\end{array}\n\\right]\n"
end

function _latex_matrix_print(io::IO, A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax}
    # Use print (not @printf) for static LaTeX: Printf can mis-parse `%` / escapes in
    # `{*{20}c}` and corrupt `\left[...]` (e.g. into `\lef[...]`).
    style = DISPLAY_STYLE[]
    print(io, "\\left[\n\\begin{array}{*{20}c}\n")
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == zero(Tropical{T})
                if style == 0
                    if A[i,j] < 0
                        print(io, "-\\infty")
                    else
                        print(io, "+\\infty")
                    end
                elseif style == 3 || style == 4
                    print(io, "\\varepsilon")
                else
                    print(io, ".")
                end
            elseif A[i,j] == one(Tropical{T})
                if style == 2 || style == 4
                    print(io, "e")
                else
                    print(io, Int64(A[i,j].λ))
                end
            elseif A[i,j].λ == trunc(A[i,j].λ)
                print(io, Int64(A[i,j].λ))
            else
                show(io, A[i,j].λ)
            end
            if j < size(A, 2)
                print(io, " & ")
            end
        end
        print(io, " \\\\\n")
    end
    print(io, "\\end{array}\n\\right]\n")
    return nothing
end

LaTeX(io::IO, A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax} = _latex_matrix_print(io, A)

# ==============================================================================
function LaTeX(S::MPSysLin)
    "\\left\\{\\begin{array}{lcl}\nx_n & = & " *
    LaTeX(S.D) * " x_n \\oplus " *
    LaTeX(S.A) * " x_{n-1} \\oplus " *
    LaTeX(S.B) * " u_n\\\\ y_n & = & " *
    LaTeX(S.C) * " x_n\\\\ x_0 & = & " *
    LaTeX(S.x0) * "\\end{array}\\right."
end

# ==============================================================================
function LaTeX(io::IO, S::MPSysLin)
    # Do not use show(io, LaTeX(S)): MIME/IO can escape or quote strings and break `\left`.
    print(io, LaTeX(S))
end

# ==============================================================================
# Scalars: Jupyter often requests `text/latex`; default to plain tropical text
# (same as notebook workaround: avoid raw LaTeX for single elements).
Base.show(io::IO, ::MIME"text/latex", x::Tropical{T}) where {T<:MinOrMax} =
    Base.show(io, MIME{Symbol("text/plain")}(), x)
Base.show(io::IO, ::Type{MIME{Symbol("text/latex")}}, x::Tropical{T}) where {T<:MinOrMax} =
    Base.show(io, MIME{Symbol("text/plain")}(), x)

# ==============================================================================
# Convert a Max-Plus dense/sparse matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on set_tropical_display(style).

function Base.show(io::IO, ::MIME"text/latex", A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax}
    print(io, "\$\$\n")
    _latex_matrix_print(io, A)
    print(io, "\$\$")
end
Base.show(io::IO, ::Type{MIME{Symbol("text/latex")}}, A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax} =
    Base.show(io, MIME{Symbol("text/latex")}(), A)

function Base.show(io::IO, ::MIME"text/latex", S::MPSysLin)
    print(io, "\$\$\n", LaTeX(S), "\n\$\$")
end
Base.show(io::IO, ::Type{MIME{Symbol("text/latex")}}, S::MPSysLin) =
    Base.show(io, MIME{Symbol("text/latex")}(), S)

function Base.show(io::IO, ::MIME"text/latex", A::LinearAlgebra.Transpose{Tropical{T}, AbstractArray{Tropical{T},N}}) where {T<:MinOrMax,N}
    print(io, "\$\$\n")
    _latex_matrix_print(io, A)
    print(io, "\$\$")
end
Base.show(io::IO, ::Type{MIME{Symbol("text/latex")}}, A::LinearAlgebra.Transpose{Tropical{T}, AbstractArray{Tropical{T},N}}) where {T<:MinOrMax,N} =
    Base.show(io, MIME{Symbol("text/latex")}(), A)
