# ==============================================================================
# Display Max-Plus scalars and matrices

"""
    mpstyle

Memory for saving the style of display for neutral and absorbing elements. By
default the display of ScicosLab will be used.
"""
global mpstyle = 1; # ScicosLab style

"""
    mp_change_display(style::Int)

Change the style of behavior of functions `Base.show()`:
- `-Inf` are displayed either with `ε` (style 2 or 3) or `.` symbols (style 1).
- `0` are displayed either with `e `(style 3) or '0' symbols (style 1 or 2).
- else: `-Inf` and `0` are displayed in Julia default sytle (style 0).

If this function is not called, by default the ScicosLab style will be used
(style 1).
"""
function mp_change_display(style::Int)
    global mpstyle = min(max(style, 0), 4)
end

#
name(::Type{MP}) = "Max-Plus "
name(::Type{MI}) = "Min-Plus "

# Base function for display a Max-Plus number.
function mpshow(io::IO, x::Trop{T}) where T
    if (mpstyle == 0)
        show(io, x.v)
    elseif x == zero(Trop{T})
        (mpstyle == 1 || mpstyle == 2) ? (@printf io ".") : (@printf io "ε")
    elseif x == one(Trop{T})
        (mpstyle == 1 || mpstyle == 3) ? (@printf io "0") : (@printf io "e")
    elseif x.v == trunc(x.v)
        (@printf io "%d" x.v)
    else
        show(io, x.v)
    end
end

function mpshow(io::IO, A::ArrTrop)
    if (size(A,2) == 1)
        print(io, size(A,1), "-element Max-Plus vector:\n")
    else
        print(io, size(A,1), '×', size(A,2), " Max-Plus dense matrix:\n")
    end
    pretty_table(io, A, tf = tf_borderless, noheader = true)
end

function mpshow(io::IO, A::LinearAlgebra.Transpose{Trop, ArrTrop})
    print(io, size(A,1), '×', size(A,2), " Max-Plus transposed dense matrix:\n")
    pretty_table(io, A, tf = tf_borderless, noheader = true)
end

function mpshow(io::IO, V::LinearAlgebra.Transpose{Trop, VecTrop})
    print(io, size(A,1), "-element Max-Plus transposed vector:\n")
    pretty_table(io, V, tf = tf_borderless, noheader = true)
end

# Note: we force the sparse display like done in Julia 1.5 because since Julia
# 1.6 sparse matrices are displayed like dense matrices with dots for zeros.
# This sounds weird since displayed huge sparse matrices take the same space
# than dense matrix.
mpshow(io::IO, S::SpaTrop) = show(io, S)

"""
    LaTeX(io::IO, A::Array{MP})

Base function for convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
neutral and absorbing elements depends on mp_change_display(style).

# Examples
```julia-repl
julia> LaTeX(stdout, MP([4 3; 7 -Inf]))
\\left[
\\begin{array}{*{20}c}
4 & 3 \\\\
7 & . \\\\
\\end{array}
\\right]
```
"""
# FIXME: not a column-major traversal
function LaTeX(io::IO, A::ArrTrop)
    (@printf io "\\left[\n\\begin{array}{*{20}c}\n")
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == zero(MP)
                if mpstyle == 0
                    (@printf io "-\\infty")
                elseif mpstyle == 3 || mpstyle == 4
                    (@printf io "\\varepsilon")
                else
                    (@printf io ".")
                end
            elseif A[i,j] == one(MP)
                if mpstyle == 2 || mpstyle == 4
                    (@printf io "e")
                else
                    (@printf io "%d" A[i,j].v)
                end
            elseif A[i,j].v == trunc(A[i,j].v)
                (@printf io "%d" A[i,j].v)
            else
                show(io, A[i,j].v)
            end
            if j < size(A, 2)
                (@printf io " & ")
            end
        end
        (@printf io " \\\\\n")
    end
    (@printf io "\\end{array}\n\\right]\n")
end

# Base function for displaying a Max-Plus sparse matrix.
LaTeX(io::IO, S::SpaTrop) = LaTeX(io, full(S))

function LaTeX(A::ArrTrop)
    s = "\\left[\n\\begin{array}{*{20}c}\n"
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == zero(MP)
                if mpstyle == 0
                    s = s * "-\\infty"
                elseif mpstyle == 3 || mpstyle == 4
                    s = s * "\\varepsilon"
                else
                    s = s * "."
                end
            elseif A[i,j] == one(MP)
                if mpstyle == 2 || mpstyle == 4
                    s = s * "e"
                else
                    s = s * string(Int64(A[i,j].v))
                end
            elseif A[i,j].v == trunc(A[i,j].v)
                s = s * string(Int64(A[i,j].v))
            else
                s = s * string(A[i,j].v)
            end
            if j < size(A, 2)
                s = s * " & "
            end
        end
        s = s * " \\\\\n"
    end
    return s * "\\end{array}\n\\right]\n"
end

"""
    show(io::IO, x::MP)

Display a Max-Plus number depending on the currently set style:

# Examples
```julia-repl
julia> mp_change_display(0); mpeye(2,2)
2×2 Max-Plus dense array:
  0.0  -Inf
 -Inf   0.0

julia> mp_change_display(1); mpeye(2,2)
2×2 Max-Plus dense array:
 0  .
 .  0

julia> mp_change_display(2); mpeye(2,2)
2×2 Max-Plus dense array:
 e  .
 .  e

julia> mp_change_display(3); mpeye(2,2)
2×2 Max-Plus dense array:
 0  ε
 ε  0

julia> mp_change_display(4); mpeye(2,2)
2×2 Max-Plus dense array:
 e  ε
 ε  e
```
"""
# Called by pretty_table() when REPL shows a MP matrix.
# julia> [MP(1) MP(2); MP(3) MP(4)]
Base.show(io::IO, x::Trop) = mpshow(io, x)

# Called by REPL when showing a MP scalar.
# julia> MP(1)
function Base.show(io::IO, ::MIME"text/plain", x::MP)
    print(io, name(MP))
    mpshow(io, x)
end

function Base.show(io::IO, ::MIME"text/plain", x::MI)
    print(io, name(MI))
    mpshow(io, x)
end

# Called by the REPL through the display() method. This function fixes
# misaliged columns made by the default show() Julia. We use the package
# PrettyTables.
# julia> [MP(1) MP(2); MP(3) MP(4)]
Base.show(io::IO, ::MIME"text/plain", A::ArrTrop) = mpshow(io, A)

# Convert a Max-Plus dense matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", A::ArrTrop) = LaTeX(io, A)

# Convert a Max-Plus sparse matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on mp_change_display(style).
Base.show(io::IO, ::MIME"text/latex", x::SpaTrop) = LaTeX(io, A)

Base.show(io::IO, ::MIME"text/plain", A::LinearAlgebra.Transpose{Trop, Matrix{Trop}}) = mpshow(io, A)
Base.show(io::IO, ::MIME"text/plain", V::LinearAlgebra.Transpose{Trop, Vector{Trop}}) = mpshow(io, V)
