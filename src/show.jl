# ##############################################################################
# Display Max-Plus scalars and matrices
# ##############################################################################

# ==============================================================================
name(::Type{MP}) = "(max,+) "
name(::Type{MI}) = "(min,+) "
name(::MP) = "(max,+) "
name(::MI) = "(min,+) "

# ==============================================================================

# Memory for saving the style of display for neutral and absorbing tropical
# elements. By default the display of ScicosLab will be used.
global display_style = 1; # Default display style: ScicosLab display style (+Inf as .)

# ==============================================================================
# Change the style of behavior of functions Base.show()
function set_tropical_display(style::Int)
    global display_style = style
    if (display_style == 0)
        (@printf stdout "I will show -Inf and 0.0")
    elseif (display_style == 1)
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as .")
    elseif (display_style == 2)
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as . and 0.0 as e")
    elseif (display_style == 3)
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as ε")
    elseif (display_style == 4)
        (@printf stdout "I will show (max,+) -Inf and (min,+) +Inf as ε and 0.0 as e")
    else
        error("style shall be 0 .. 4")
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) scalar.
function tropshow(io::IO, x::Tropical{T}) where {T<:MinOrMax}
    if (display_style == 0)
        show(io, x.λ)
    elseif x == zero(x)
        (display_style == 1 || display_style == 2) ? (@printf io ".") : (@printf io "ε")
    elseif x == one(x)
        (display_style == 1 || display_style == 3) ? (@printf io "0") : (@printf io "e")
    elseif x.λ == trunc(x.λ)
        (@printf io "%d" x.λ)
    else
        show(io, x.λ)
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) dense array. We use pretty_table
# to fix misaligned elements when mixing Inf, ε, e, numbers.
function tropshow(io::IO, A::AbstractArray{Tropical{T},N}) where {T<:MinOrMax,N}
    if (N == 1)
        print(io, size(A,1), "-element ", name(Tropical{T}), "vector:\n")
        pretty_table(io, A, tf = tf_borderless, noheader = true)
    elseif (N == 2)
        print(io, size(A,1), '×', size(A,2), " ", name(Tropical{T}), "dense matrix:\n")
        pretty_table(io, A, tf = tf_borderless, noheader = true)
    else
        show(io, A)
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) transposed matrix.
function tropshow(io::IO, A::LinearAlgebra.Transpose{Tropical{T}, Array{Tropical{T},N}}) where {T<:MinOrMax,N}
    if (N == 1)
        print(io, size(A,1), "-element ", name(A), "transposed vector:\n")
        pretty_table(io, A, tf = tf_borderless, noheader = true)
    elseif (N == 2)
        print(io, size(A,1), '×', size(A,2), " ", name(A), "transposed dense matrix:\n")
        pretty_table(io, A, tf = tf_borderless, noheader = true)
    else
        show(io, A)
    end
end

# ==============================================================================
# Base function for display a (max,+) or (min,+) transposed vector.
function tropshow(io::IO, V::LinearAlgebra.Transpose{Tropical{T}, Vector{Tropical{T}}}) where {T<:MinOrMax}
    print(io, size(A,1), "-element ", name(A), "transposed vector:\n")
    pretty_table(io, V, tf = tf_borderless, noheader = true)
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
Base.show(io::IO, x::Tropical) = tropshow(io, x)
Base.show(io::IO, S::MPSysLin) = tropshow(io, S)
Base.show(io::IO, ::MIME"text/plain", x::Tropical{T}) where {T<:MinOrMax} = (print(io, name(x)); tropshow(io, x))

# ==============================================================================
# Called by the REPL through the display() method. These functions fix
# misaliged columns made by the default show() Julia. We use the package
# PrettyTables.

Base.show(io::IO, ::MIME"text/plain", A::VecOrMat{Tropical{T}}) where {T<:MinOrMax} = tropshow(io, A)
Base.show(io::IO, ::MIME"text/plain", S::MPSysLin) = tropshow(io, S)
Base.show(io::IO, ::MIME"text/plain", A::LinearAlgebra.Transpose{Tropical{T}, AbstractArray{Tropical{T},N}}) where {T<:MinOrMax,N} = tropshow(io, A)

# ==============================================================================
#
function LaTeX(A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax}
    s = "\\left[\n\\begin{array}{*{20}c}\n"
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == zero(Tropical{T})
                if display_style == 0
                    if A[i,j] < 0
                        s = s * "-\\infty"
                    else
                        s = s * "+\\infty"
                    end
                elseif display_style == 3 || display_style == 4
                    s = s * "\\varepsilon"
                else
                    s = s * "."
                end
            elseif A[i,j] == one(Tropical{T})
                if display_style == 2 || display_style == 4
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

#FIXME Does not pass unit test. Why ? LaTeX(io::IO, A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax} = show(io, LaTeX(A))
#FIXME: temporary
function LaTeX(io::IO, A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax}
    (@printf io "\\left[\n\\begin{array}{*{20}c}\n")
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            if A[i,j] == zero(Tropical{T})
                if display_style == 0
                    if A[i,j] < 0
                        (@printf io "-\\infty")
                    else
                        (@printf io "+\\infty")
                    end
                elseif display_style == 3 || display_style == 4
                    (@printf io "\\varepsilon")
                else
                    (@printf io ".")
                end
            elseif A[i,j] == one(Tropical{T})
                if display_style == 2 || display_style == 4
                    (@printf io "e")
                else
                    (@printf io "%d" A[i,j].λ)
                end
            elseif A[i,j].λ == trunc(A[i,j].λ)
                (@printf io "%d" A[i,j].λ)
            else
                show(io, A[i,j].λ)
            end
            if j < size(A, 2)
                (@printf io " & ")
            end
        end
        (@printf io " \\\\\n")
    end
    (@printf io "\\end{array}\n\\right]\n")
end

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
LaTeX(io::IO, S::MPSysLin) = show(io, LaTeX(S))

# ==============================================================================
# Convert a Max-Plus dense/sparse matrix to a LaTeX formula. Symbols of
# neutral and absorbing elements depends on set_tropical_display(style).

Base.show(io::IO, ::MIME"text/latex", A::AbstractVecOrMat{Tropical{T}}) where {T<:MinOrMax} = LaTeX(io, A)
Base.show(io::IO, ::MIME"text/latex", S::MPSysLin) = LaTeX(io, S)
Base.show(io::IO, ::MIME"text/latex", A::LinearAlgebra.Transpose{Tropical{T}, AbstractArray{Tropical{T},N}}) where {T<:MinOrMax,N} = LaTeX(io, A)
