![logo](https://lecrapouille.github.io/icons/juliamaxplus.png)

# Julia (max,+) and (min,+) Algebra Toolbox

[![CI](https://github.com/Lecrapouille/MaxPlus.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Lecrapouille/MaxPlus.jl/actions/workflows/CI.yml) [![codecov](https://codecov.io/gh/Lecrapouille/MaxPlus.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Lecrapouille/MaxPlus.jl) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lecrapouille.github.io/MaxPlus.jl)

This [MaxPlus.jl](https://github.com/Lecrapouille/MaxPlus.jl) package is a Julia port of the [ScicosLab](http://www.scicoslab.org/) (max,+) toolbox (no longer maintained), which provides functions for numerical computations in (max,+) algebra. This Julia toolbox extends the original toolbox by adding (min,+) algebra support.

The (max,+) algebra (max-plus algebra) redefines the operators `plus` and `times` from classical algebra as operators `maximum` (symbolized by ⨁) and `plus` (symbolized by ⨂) in the domain of real numbers ℝ augmented by `-∞`.

The (min,+) algebra (min-plus algebra) redefines operators `plus` and `times` from classical algebra as operators `minimum` (symbolized by ⨁) and `plus` (symbolized by ⨂) in the domain of real numbers ℝ augmented by `+∞`.

Matrix computation in this algebra has been taught since 1960 by J. Kuntzman in his theory of networks. It is used in numerous domains such as operations research (network theory), physics (quantization), probability theory (Cramér's transform), control theory (discrete event systems), computer science (automata theory, Petri nets), and mathematics (algebraic geometry). This algebra is also known as `tropical algebra`.

## Requirements

- **Julia 1.10 or later** (LTS recommended). Older Julia version may suffer of invalid results (sparse matrices)!
- No mandatory external Julia dependencies.
- All Julia standard libraries: `LinearAlgebra`, `SparseArrays`, `Printf`

You may also be interested in this [graphical Petri editor](https://github.com/Lecrapouille/TimedPetriNetEditor), which can generate (max,+) matrices from timed event graphs.

Scalars are `Tropical{Sense,T} <: Real` with `T <: AbstractFloat` (default `Float64`), so you can mix them with much of the Julia `LinearAlgebra` ecosystem; `ε` (`TropicalZero`) is a neutral placeholder compatible with both (max,+) and (min,+).Scalars are `Tropical{Sense,T} <: Real` with `T <: AbstractFloat` (default `Float64`), so you can mix them with much of the Julia `LinearAlgebra` ecosystem; `ε` (`TropicalZero`) is a neutral placeholder compatible with both (max,+) and (min,+).

## Installation

### From Julia Package Registry (stable version)

```julia
julia> ]
pkg> add MaxPlus
```

### From Source (latest development version)

```bash
git clone https://github.com/Lecrapouille/MaxPlus.jl.git --depth=1
cd MaxPlus.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Then from the Julia REPL:

```julia
julia> ]
pkg> dev .
```

## Quick Start

Full documentation about this algebra and its implementation in Julia is available at: [https://lecrapouille.github.io/MaxPlus.jl](https://lecrapouille.github.io/MaxPlus.jl/index.html).

### (max,+) Algebra

In (max,+) algebra, `+` computes the maximum and `*` computes the sum.

```julia
julia> using MaxPlus

julia> MP(5)
(max,+) 1

julia> MP(1) + MP(5)
(max,+) 5

julia> MP(1) * MP(5)
(max,+) 6
```

Here an example on dense matrice:

```julia
julia> MP([1 2; 3 8]) .+ 5
2×2 (max,+) dense matrix:
  5   5
  5   8
```

The `.+` is the Julia element-wise addition operator (broadcasted maximum in (max,+) algebra). The result was computed as:

```julia
[max(1, 5)  max(2, 5)
 max(3, 5)  max(8, 5)]
```

This package also works with sparse matrices where `-∞` are not stored:

```julia
julia> using SparseArrays

julia> S = sparse(MP([1 2; -Inf 4]))
2×2 (max,+) sparse matrix with 3 stored entries:
  [1, 1]  =  1
  [1, 2]  =  2
  [2, 2]  =  4
```

To use MaxPlus matrices with `Graphs.jl`, convert to classical algebra first:

```julia
using Graphs, SimpleWeightedGraphs
A = MP([3.0 7; 2 4])
g = SimpleWeightedDiGraph(A)
```

### (min,+) Algebra

```julia
julia> MI([1 2; 3 8])
2×2 (min,+) dense matrix:
  1   2
  3   8
```

### Constants

| (max,+) | (min,+) | Description |
|---------|---------|-------------|
| `ε` | `ε` | Zero (absorbing element): -∞ for (max,+), +∞ for (min,+) |
| `mp0` | `mi0` | Zero (absorbing element): -∞ for (max,+), +∞ for (min,+) |
| `mp1` / `mpe` | `mi1` / `mie` | One (neutral element): 0 for both algebras |
| `mptop` | `mitop` | Top element: +∞ for (max,+), -∞ for (min,+) |
| `mpI` | `miI` | Identity scaling operator |

### Spectral Analysis (Howard Algorithm)

```julia
julia> using SparseArrays

julia> A = MP([3.0 7; 2 4])
julia> λ, v = mpeigen(sparse(A))
(MP[4.5, 4.5], MP[6.5, 4.0])

julia> A * v == λ[1] * v
true
```

The package also provides `semihoward` for semi-Markov Howard algorithm with time weights.

### Linear Systems

```julia
julia> S = MPSysLin(MP([1 2; 3 4]), MP([1; 1]), MP([1 1]))
julia> [S; S]  # Vertical composition (vcat)
julia> [S S]   # Horizontal composition (hcat)
```

## Documentation

Full documentation is available at:
[https://lecrapouille.github.io/MaxPlus.jl](https://lecrapouille.github.io/MaxPlus.jl/index.html)

Contents:

- [Introduction and tutorials](tutorial) (French and English)
- [API: (max,+) functions](docs/src/maxplus.md)
- [API: (min,+) functions](docs/src/minplus.md)
- [API: Linear systems](docs/src/syslin.md)
- [ScicosLab to MaxPlus.jl translation](docs/src/portage.md)
- [Running tests](docs/src/tests.md)
- [Bibliography](docs/src/bibliography.md)

## Related Projects

- [TimedPetriNetEditor](https://github.com/Lecrapouille/TimedPetriNetEditor): A graphical editor for Timed Petri Nets and Event Graphs with (max,+) algebra integration.

## Contributing

Contributions are welcome! Feel free to:

- Report issues
- Submit pull requests
- Improve documentation
- Add tests

## License

This project is licensed under the GPL-3.0 License.
