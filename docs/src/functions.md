# Scilab and Julia function names

## Max+ core functions

| Scilab               | Julia                         | Comment                                                                                                                    |
|----------------------|-------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| full                 | full, dense, Array            | Transform a sparse max+ matrix to a dense max+ matrix.                                                                     |
| howard               | todo                          | Get eigenvalues and eigenvectors.                                                                                          |
| karp                 |                               | howard is more complete and faster than karp.                                                                              |
|                      | LaTeX                         | Output a max+ matrix to LaTeX formula.                                                                                     |
| maxplus, #           | MP                            | Create a max+ number or dense array or sparse array.                                                                       |
| min, &               | Base.min                      | Return the min of two max+ numbers.                                                                                        |
| minus, -             | Base.:(-)                     |                                                                                                                            |
|                      | minplus                       | Transform the max+ number into a min+ number.                                                                              |
| mnorm                | mpnorm                        | Compute the norm of a max+ matrix.                                                                                         |
| mptrace              | mptrace                       | Compute the trace of a max+ matrix.                                                                                        |
| percent0, %0         | Base.zero, mpzero, mp0        | Return the max+ element zero.                                                                                              |
| percent1, %1         | Base.one, mpone, mp1          | Return the max+ element one.                                                                                               |
| percenteye, %eye     | mpeye                         | Return the max+ identity matrix.                                                                                           |
|                      | mpI                           | The diagonal elements used for creating max+ identity matrix.                                                              |
| percentones, %ones   | mpones                        | Return the zero'ed max+ matrix.                                                                                            |
| percenttop, %top     | mptop                         | Return the min+ top value.                                                                                                 |
| percentzeros, %zeros | mpzeros                       | Return the one'd max+ matrix.                                                                                              |
| plus, +              | Base.:(+)                     | max+ plus operator.                                                                                                        |
| plustimes            | plustimes                     | Convert a max+ number or matrix to a standard number or matrix (where plus and times operators are the classic operators). |
| power, ^             | Base.:(^)                     | Compute the power of a max+ number or matrix.                                                                              |
|                      | Base.show                     | Display a max+ number or matrix.                                                                                           |
|                      | SparseArrays.sparse, mpsparse | Create a max+ sparse matrix.                                                                                               |
| star                 | mpstar                        | Compute the A* of the max+ matrix.                                                                                         |
| times, *             | Base.:(*)                     | max+ times operator.                                                                                                       |
| typeof               | typeof                        | Return the type of the max+ number.                                                                                        |

## Dynamic linear maxplus system

| Scilab   | Julia            | Comment                                                                                                      |
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

| Scilab              | Julia | Comment                                         |
|---------------------|-------|-------------------------------------------------|
| flowshop            |       | Flowshop maxplus linear system builder.         |
|                     | LaTeX | Generate the flowshop graph as a LaTeX formula. |
| flowshop_graph      |       | Flowshop graph builder.                         |
| flowshop_simu       |       | Maxplus linear system simulation.               |
| saturation_graph    |       | Saturation graph.                               |
| show_cr_graph       |       | Flowshop critical graph.                        |
| strong_connex_graph |       | Strong connected component extracting.          |

## Incompatibility

| Function                      | Scilab | NSP | Julia   |
|-------------------------------|--------|-----|---------|
| length(mpzeros(Float64, 2,2)) | 2      | 4   | 4       |


