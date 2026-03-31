![logo](https://lecrapouille.github.io/icons/juliamaxplus.png)

# Julia (max,+) and (min,+) Algebra Toolbox

[![CI](https://github.com/Lecrapouille/MaxPlus.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Lecrapouille/MaxPlus.jl/actions/workflows/CI.yml) [![codecov](https://codecov.io/gh/Lecrapouille/MaxPlus.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Lecrapouille/MaxPlus.jl) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lecrapouille.github.io/MaxPlus.jl)

The (max,+) algebra (also known as max-plus algebra) redefines the operators `plus` and `times` from classical algebra as operators `maximum` (symbolized by ⨁) and `plus` (symbolized by ⨂) in the domain of real numbers ℝ augmented by minus infinity -∞. The (min,+) algebra (min-plus algebra) redefines operators plus and times from classical algebra as operators minimum (symbolized by ⨁) and plus (symbolized by ⨂) in the domain of real numbers ℝ augmented by plus infinity +∞.

Matrix computation in this algebra has been taught since 1960 by J. Kuntzman in his theory of networks. It is used in numerous domains such as operations research (network theory), physics (quantization), probability theory (Cramér's transform), control theory (discrete event systems), computer science (automata theory, Petri nets), and mathematics (algebraic geometry). This algebra is also known as `tropical algebra`.

This [MaxPlus.jl](https://github.com/Lecrapouille/MaxPlus.jl) package is a Julia port of the [ScicosLab](http://www.scicoslab.org/) (max,+) toolbox (no longer maintained), which provides functions for numerical computations in (max,+) algebra. This Julia toolbox extends the original toolbox by adding (min,+) algebra support.

You may also be interested in this Timed Petri Net and Timed Event Graph graphical [editor](https://github.com/Lecrapouille/TimedPetriNetEditor), which can generate (max,+) matrices from timed event graphs.

## Requirements

- **Julia 1.10 or later** (LTS recommended). Older Julia version may suffer of invalid results!
- All Julia standard libraries: `LinearAlgebra`, `SparseArrays`, `Printf` (no external dependencies).

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

### (max,+) Algebra

```julia
julia> using MaxPlus

julia> MP([1 2; 3 8]) .+ 5
2×2 (max,+) dense matrix:
  5   5
  5   8
```

The result is computed as:

```julia
[max(1, 5)  max(2, 5)
 max(3, 5)  max(8, 5)]
```

Note: In (max,+) algebra, `+` computes the maximum and `*` computes the sum. The `.+` is the Julia element-wise addition operator (broadcasted maximum in (max,+) algebra).

### (min,+) Algebra

```julia
julia> MI([1 2; 3 8])
2×2 (min,+) dense matrix:
  1   2
  3   8
```

### Type Aliases

- `MP` is an alias for `MaxPlus` (Tropical{Max})
- `MI` is an alias for `MinPlus` (Tropical{Min})

### Constants

| (max,+) | (min,+) | Description |
|---------|---------|-------------|
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

## Interoperability with Graphs.jl

To use MaxPlus matrices with `Graphs.jl`, convert to classical algebra first:

```julia
using Graphs, SimpleWeightedGraphs
A = MP([3.0 7; 2 4])
g = SimpleWeightedDiGraph(plustimes(A))  # Convert to Float64 matrix
```

## Contributing

Contributions are welcome! Feel free to:

- Report issues
- Submit pull requests
- Improve documentation
- Add tests

## License

This project is licensed under the GPL-3.0 License.
