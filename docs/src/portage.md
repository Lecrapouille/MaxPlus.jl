# Rosetta Stone: Mapping Scilab Function Names to Julia

## (min,+) Core Functions

Scilab/ScicosLab/NSP do not natively support (min,+) algebra.

## (max,+) Core Functions

In Julia, you can type `?` followed by a function's name in the REPL to access its documentation and usage examples. The following table lists functions in the alphabetical order of their ScicosLab names.

| ScicosLab            | Julia                         | Description                                                                                                  |
|----------------------|-------------------------------|--------------------------------------------------------------------------------------------------------------|
| astarb               | [`astarb`](@ref)              | Solves for x in the equation x = A x + b                                                                     |
| / (division)         | Base.:(/)                     | (max,+) division operator for scalars, dense, or sparse matrices                                             |
| full                 | [`full`](@ref), [`dense`](@ref), Array | Converts a (max,+) sparse matrix to a (max,+) dense matrix                                |
| howard               | [`howard`](@ref), [`mpeigen`](@ref)   | Computes eigenvalues and eigenvectors from a sparse matrix using the Howard algorithm                        |
| inv                  | [`Base.inv`](@ref), Base.:(^,-1)      | Computes the inverse of a (max,+) matrix                                                                     |
| karp                 | [`howard`](@ref), [`mpeigen`](@ref)   | Karp algorithm not implemented; Howard algorithm preferred for speed                                         |
|                      | [`LaTeX`](@ref), [`show`](@ref) (MIME"text/latex") | Outputs a (max,+) matrix as LaTeX code                                              |
| # (maxplus)          | [`MP`](@ref)                  | Creates a (max,+) number, dense array, or sparse array                                                       |
| & (min)              | Base.min                      | Returns the minimum of two (max,+) scalars or matrices                                                       |
| - (minus)            | Base.:(-)                     | (max,+) minus (unary and binary) for scalar, dense or sparse matrix                                          |
| mnorm                | [`norm`](@ref)                | Computes the norm of a (max,+) dense or sparse matrix                                                        |
| mptrace              | [`tr`](@ref)                  | Computes the trace of a (max,+) dense or sparse matrix                                                       |
| percent0, %0         | [`Base.zero`](@ref), [`zero`](@ref), [`mp0`](@ref)    | Returns the (max,+) zero element                                        |
| percent1, %1         | [`Base.one`](@ref), [`one`](@ref), [`mp1`](@ref), [`mpe`](@ref) | Returns the (max,+) one element                                 |
| percenteye, %eye     | [`eye`](@ref)                 | Returns the (max,+) identity matrix                                                                          |
|                      | [`mpI`](@ref)                 | Used for creating (max,+) identity matrices (diagonal elements)                                              |
| percentones, %ones   | [`ones`](@ref)                | Returns the (max,+) matrix or column vector of ones                                                          |
| percenttop, %top     | [`mptop`](@ref)               | Returns the min-plus absorbing element (interpreted as infinity in max-plus context)                         |
| percentzeros, %zeros | [`spzeros`](@ref)             | Returns a zero (max,+) sparse matrix or column vector (built-in in Julia)                                    |
|                      | [`zeros`](@ref)               | Returns a zero (max,+) dense matrix or column vector                                                         |
| + (plus)             | [`Base.:(+)`](@ref)           | (max,+) addition for scalars, dense or sparse matrices                                                       |
| plus                 | [`plus`](@ref)                | Computes the (max,+) matrix pseudo-inverse A^+                                                               |
| plustimes            | [`plustimes`](@ref)           | Converts (max,+) numbers or matrices to classic (standard arithmetic) numbers or matrices                    |
| ^ (power)            | [`Base.:(^)`](@ref)           | Computes (max,+) powers of numbers or matrices                                                               |
|                      | [`Base.show`](@ref)           | Displays a (max,+) number or a sparse/dense matrix                                                           |
|                      | [`SparseArrays.sparse`](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#SparseArrays.sparse) | Creates a (max,+) sparse matrix or converts a dense matrix                  |
| \ (residu)           | [`Base.:(\)`](@ref)           | Solves for x in the equation A x = b                                                                         |
| semihoward           | `semihoward`                  | Semi-Markov Howard (`semihoward(S, Tau)`), voir l’API `(max,+)`.                                            |
| star                 | [`star`](@ref)                | Computes the (max,+) matrix closure A^*                                                                      |
| * (times)            | [`Base.:(*)`](@ref)           | (max,+) multiplication for scalars, dense or sparse matrices                                                 |
| typeof               | [`Base.typeof`](https://docs.julialang.org/en/v1/base/base/#Core.typeof)         | Returns the type of a (max,+) number                                    |

## Dynamic Linear (max,+) Systems

| ScicosLab      | Julia                        | Description                                                                          |
|----------------|------------------------------|--------------------------------------------------------------------------------------|
| mpsyslin       | [`MPSysLin`](@ref)           | Structure holding the state matrices of the system (function in Scilab, type in Julia)|
|                | [`LaTeX`](@ref)              | Outputs a (max,+) system in LaTeX format                                             |
| simul          | [`simul`](@ref)            | Simulates a system (implicit form; uses `star(D)` per call)                          |
| explicit       | [`explicit`](@ref)         | Converts to explicit state-space by removing implicit part                            |
| implicit       | [`implicit`](@ref), `implicit` | Canonical implicit form `D = I` with `A ← star(D)*A`, `B ← star(D)*B` (same I/O as `simul`). |
| full(s) mpls   | `mpfull`                     | Converts system matrices to dense representation (Julia deprecated `Base.full`)       |
| sparse(s) mpls | `mpsparse`                   | Converts system matrices to sparse representation                                     |
|                | Base.show                    | Displays the (max,+) system                                                          |
| plus, +        | `Base.:(+)`                  | (max,+) plus operator or combination with scalar/matrix                               |
| times, *       | `Base.:(*)`                  | Series (sequential) composition                                                      |
| \|             | `Base.:(|)`                  | Diagonal (direct sum) composition                                                    |
| /.             | `Base.:(/)`                  | Feedback composition                                                                 |
| [;]            | Base.:vcat, [;]              | Concatenates inputs (inputs in common)                                               |
| [,]            | Base.:hcat, [,]              | Concatenates outputs (outputs in addition)                                           |
| typeof         | typeof                       | Returns the type of the (max,+) system                                               |

## Flowshop

| ScicosLab           | Julia           | Description                                                           |
|---------------------|-----------------|-----------------------------------------------------------------------|
| flowshop            | `flowshop`      | E matrix (machines × jobs in (max,+)); use `mp0` for missing tasks    |
| shift               | `mpshift`       | n-shift system / shift by t (see `shift.sci`)                         |
| flowshop_graph      | `flowshop_graph` | Only `(T, N)` sparse `(max,+)` matrices (no graph / `show_graph`).   |
| flowshop_simu       | `flowshop_simu` | Flowshop system simulation: s = `MPSysLin`; u: rows=inputs, columns=time|
| saturation_graph    | TODO            | Saturation graph (not yet implemented)                                |
| show_cr_graph       | TODO            | Flowshop critical graph (not yet implemented)                         |
| strong_connex_graph | TODO            | Extracts strongly connected components (not yet implemented)          |

## Incompatibilities

| Function             | ScicosLab       | NSP             | Julia       | Comment                 |
|----------------------|-----------------|-----------------|-------------|-------------------------|
| length(mpzeros(2,2)) | 2               | 4               | 4           | 4 is the correct result |
| norm, trace, min     | not implemented | not implemented | implemented |                         |
| MP^(-x)              | not implemented | not implemented | implemented | Negative powers (x > 0) |
