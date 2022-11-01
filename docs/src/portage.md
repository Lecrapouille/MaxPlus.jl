# Rosetta stone for converting Scilab function names to Julia

## (min,+) core functions

Scilab/ScicosLab/NSP does not implement (min,+) algebra natively.

## (max,+) core functions

In the Julia REPL, type `?` and the name of the function to show its documentation and examples. Functions are given in the alphabteical order of ScicosLab functions.

| ScicosLab            | Julia                         | Comment                                                                                                          |
|----------------------|-------------------------------|------------------------------------------------------------------------------------------------------------------|
| astarb               | [`astarb`](@ref)              | Return the x of the equation x = A x + b                                                                         |
| / (division)         | Base.:(/)                     | (max,+) division operator on scalar, dense or sparse matrix.                                                     |
| full                 | [`full`](@ref), [`dense`](@ref), Array            | Transform a (max,+) sparse matrix to a (max,+) dense matrix.                                 |
| howard               | [`howard`](@ref) [`mpeigen`](@ref)   | Return eigenvalues and eigenvectors from a sparse matrix using Howard algorithm.                          |
| inv                  | [`Base.inv`](@ref), Base.:(^,-1)   | Compute the inverse of the (max,+) matrix                                                                   |
| karp                 | [`howard`](@ref) [`mpeigen`](@ref)   | Karp algorithm is not implemented since Howard algorithm is faster.                                       |
|                      | [`LaTeX`](@ref), [`show`](@ref) MIME"text/latex"  | Output a (max,+) matrix to LaTeX code.                                                       |
| # (maxplus)          | [`MP`](@ref)                  | Create a (max,+) number or dense array or sparse array.                                                          |
| & (min)              | Base.min                      | Return the min of two (max,+) scalars or matrix.                                                                 |
| - (minus)            | Base.:(-)                     | (max,+) minus unary and binary operators on scalar, dense or sparse matrix.                                      |
| mnorm                | [`norm`](@ref)                | Compute the norm of a (max,+) dense or sparse matrix.                                                            |
| mptrace              | [`tr`](@ref)                  | Compute the trace of a (max,+) dense or sparse matrix.                                                           |
| percent0, %0         | [`Base.zero`](@ref), [`zero`](@ref), [`mp0`](@ref), [`ϵ`](@ref)     | Return the (max,+) element zero.                                           |
| percent1, %1         | [`Base.one`](@ref), [`one`](@ref), [`mp1`](@ref),[` mpe`](@ref)     | Return the (max,+) element one.                                            |
| percenteye, %eye     | [`eye`](@ref)                 | Return the (max,+) identity matrix.                                                                              |
|                      | [`mpI`](@ref)                 | The diagonal elements used for creating (max,+) identity matrix.                                                 |
| percentones, %ones   | [`ones`](@ref)                | Return the once (max,+) matrix or column vector.                                                                 |
| percenttop, %top     | [`mptop`](@ref)               | Return the min+ element zero.                                                                                    |
| percentzeros, %zeros | [`spzeros`](@ref)             | Return the zero (max,+) sparse matrix or column vector (Julia builtin).                                          |
|                      | [`zeros`](@ref)               | Return the zero (max,+) dense matrix or column vector.                                                           |
| + (plus)             | [`Base.:(+)`](@ref)           | (max,+) plus operator on scalar, dense or sparse matrix.                                                         |
| plus                 | [`plus`](@ref)                | Compute the (max,+) matrix A^+                                                                                   |
| plustimes            | [`plustimes`](@ref)           | Convert a (max,+) number or matrix to a standard number or matrix (where plus and times operators are the classic operators). |
| ^ (power)            | [`Base.:(^)`](@ref)           | Compute the power of a (max,+) number or matrix.                                                                 |
|                      | [`Base.show`](@ref)           | Display a (max,+) number or sparse or dense matrix.                                                              |
|                      | [`SparseArrays.sparse`](@ref) | Create a (max,+) sparse matrix or convert a (max,+) dense matrix.                                                |
| \ (residu)           | [`Base.(\)`](@ref)            | Return the x of the equation A x = b                                                                             |
| semihoward           | [`semihoward`](@ref)          | Critical cycles                                                                                                  |
| star                 | [`star`](@ref)                | Compute the (max,+) matrix A^*                                                                                   |
| * (times)            | [`Base.:(*)`](@ref)           | (max,+) times operator on scalar, dense or sparse matrix.                                                        |
| typeof               | [`typeof`](@ref)              | Return the type of the (max,+) number.                                                                           |

## Dynamic linear (max,+) system

| ScicosLab| Julia                 | Comment                                                                                            |
|----------|-----------------------|----------------------------------------------------------------------------------------------------|
|          | [`MPSysLin`](@ref)    | Structure holding state matrices of the syslin.                                                    |
| mpsyslin | [`mpsyslin`](@ref)    | Function building a MPSysLin structure.                                                            |
|          | [`LaTeX`](@ref)       | Output a (max,+) syslin to a LaTeX formula.                                                        |
| simul    | [`mpsimul`](@ref)     | Simulation of the (max,+) linear system.                                                           |
| explicit | [`mpexplicit`](@ref)  | Conversion of an implicit dynamic linear (max,+) system to an explicit form (where D state matrix is zero'ed).  |
|          | Base.show             | Display the (max,+) syslin.                                                                        |
| plus, +  | Base.:(+)             | (max,+) plus operator, or product with a (max,+) scalar or (max,+) matrix.                         |
| times, * | Base.:(*)             | Series composition.                                                                                |
| |        | Base.:(|)             | Diagonal composition.                                                                              |
| /.       | Base.:(/)             | Feedback composition.                                                                              |
| [;]      | Base.:vcat, [;]       | Inputs in common and concatenation of outputs.                                                     |
| [,]      | Base.:hcat, [,]       |  Concatenation of inputs and addition of outputs.                                                  |
| typeof   | typeof                | Return the type of the (max,+) syslin.                                                             |

## Flowshop

| ScicosLab           | Julia | Comment                                         |
|---------------------|-------|-------------------------------------------------|
| flowshop            |       | Flowshop maxplus linear system builder.         |
|                     | LaTeX | Generate the flowshop graph as a LaTeX formula. |
| flowshop_graph      |       | Flowshop graph builder.                         |
| flowshop_simu       |       | Maxplus linear system simulation.               |
| saturation_graph    |       | Saturation graph.                               |
| show_cr_graph       |       | Flowshop critical graph.                        |
| strong_connex_graph |       | Strong connected component extracting.          |

## Incompatibility

| Function             | ScicosLab       | NSP             | Julia       | Comment                 |
|----------------------|-----------------|-----------------|-------------|-------------------------|
| length(mpzeros(2,2)) | 2               | 4               | 4           | 4 is the correct answer |
| norm, trace, min     | not implemented | not implemented | implemented |                         |
| MP^(-x)              | not implemented | not implemented | implemented | Negative power (x > 0)  |
