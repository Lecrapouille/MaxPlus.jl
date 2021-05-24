# Scilab and Julia function names

## Max+ core functions

In the Julia REPL, type `?` and the name of the function to show its documentation and examples. Functions are given in the alphabteical order of ScicosLab functions.

| ScicosLab            | Julia                         | Comment                                                                                                                    |
|----------------------|-------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| astarb               | mpastarb                      | Return the x of the equation x = A x + b                                                                                   |
| / (division)         | Base.:(/)                     | max+ division operator on scalar, dense or sparse matrix.                                                                  |
| full                 | full, dense, Array            | Transform a max+ sparse matrix to a max+ dense matrix.                                                                     |
| howard               | howard                        | Return eigenvalues and eigenvectors from, ideally, a sparse matrix or dense matrix or a scalar.                            |
| inv                  | Base.inv, Base.:(^,-1)        | Compute the inverse of the max+ matrix                                                                                     |
| karp                 |                               | Not implemented since howard is faster than karp.                                                                          |
|                      | LaTeX, show MIME"text/latex"  | Output a max+ matrix to LaTeX code.                                                                                        |
| # (maxplus)          | MP                            | Create a max+ number or dense array or sparse array.                                                                       |
| & (min)              | Base.min                      | Return the min of two max+ scalars or matrix.                                                                              |
| - (minus)            | Base.:(-)                     | max+ minus unary and binary operators on scalar, dense or sparse matrix.                                                   |
| mnorm                | mpnorm                        | Compute the norm of a max+ dense or sparse matrix.                                                                         |
| mptrace              | mptrace                       | Compute the trace of a max+ dense or sparse matrix.                                                                        |
| percent0, %0         | Base.zero, mpzero, mp0, ϵ     | Return the max+ element zero.                                                                                              |
| percent1, %1         | Base.one, mpone, mp1, mpe     | Return the max+ element one.                                                                                               |
| percenteye, %eye     | mpeye                         | Return the max+ identity matrix.                                                                                           |
|                      | mpI                           | The diagonal elements used for creating max+ identity matrix.                                                              |
| percentones, %ones   | mpones                        | Return the once max+ matrix or column vector.                                                                              |
| percenttop, %top     | mptop                         | Return the min+ element zero.                                                                                              |
| percentzeros, %zeros | mpzeros                       | Return the zero max+ sparse matrix or column vector.                                                                       |
| + (plus)             | Base.:(+)                     | max+ plus operator on scalar, dense or sparse matrix.                                                                      |
| plus                 | mpplus                        | Compute the max+ matrix A^+                                                                                                |
| plustimes            | plustimes                     | Convert a max+ number or matrix to a standard number or matrix (where plus and times operators are the classic operators). |
| ^ (power)            | Base.:(^)                     | Compute the power of a max+ number or matrix.                                                                              |
|                      | Base.show                     | Display a max+ number or sparse or dense matrix.                                                                           |
|                      | SparseArrays.sparse           | Create a max+ sparse matrix or convert a max+ dense matrix.                                                                |
| \ (residu)           | Base.(\)                      | Return the x of the equation A x = b                                                                                       |
| semihoward           | semihoward                    | TODO                                                                                                                       |
| star                 | mpstar                        | Compute the max+ matrix A^*                                                                                                |
| * (times)            | Base.:(*)                     | max+ times operator on scalar, dense or sparse matrix.                                                                     |
| typeof               | typeof                        | Return the type of the max+ number.                                                                                        |

## Dynamic linear maxplus system

| ScicosLab| Julia            | Comment                                                                                                      |
|----------|------------------|--------------------------------------------------------------------------------------------------------------|
|          | MPSysLin         | Structure holding state matrices of the syslin.                                                              |
| mpsyslin | mpsyslin         | Function building a MPSysLin structure.                                                                      |
|          | LaTeX            | Output a max+ syslin to a LaTeX formula.                                                                     |
| simul    | mpsimul          | Simulation of the max+ linear system.                                                                        |
| explicit | mpexplicit       | Conversion of an implicit dynamic linear max+ system to an explicit form (where D state matrix is zero'ed).  |
|          | Base.show        | Display the max+ syslin.                                                                                     |
| plus, +  | Base.:(+)        | max+ plus operator, or product with a max+ scalar or max+ matrix.                                            |
| times, * | Base.:(*)        | Series composition.                                                                                          |
| |        | Base.:(|)        | Diagonal composition.                                                                                        |
| /.       | Base.:(/)        | Feedback composition.                                                                                        |
| [;]      | Base.:vcat, [;]  | Inputs in common and concatenation of outputs.                                                               |
| [,]      | Base.:hcat, [,]  |  Concatenation of inputs and addition of outputs.                                                            |
| typeof   | typeof           | Return the type of the max+ syslin.                                                                          |

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

| Function             | ScicosLab | NSP | Julia   | Comment                 |
|----------------------|-----------|-----|---------|-------------------------|
| length(mpzeros(2,2)) | 2         | 4   | 4       | 4 is the correct answer |
