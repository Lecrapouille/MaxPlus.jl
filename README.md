![logo](https://lecrapouille.github.io/icons/juliamaxplus.png)

# Julia's (max,+) and (min,+) Algebra Toolbox

[![](https://travis-ci.org/Lecrapouille/MaxPlus.jl.svg?branch=master)](https://travis-ci.org/Lecrapouille/MaxPlus.jl)
[![](https://coveralls.io/repos/github/Lecrapouille/MaxPlus.jl/badge.svg?branch=master)](https://coveralls.io/github/Lecrapouille/MaxPlus.jl?branch=master)
[![](https://codecov.io/gh/Lecrapouille/MaxPlus.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Lecrapouille/MaxPlus.jl)

The (max,+) algebra (or simply max-plus) redefines operators plus and times from
the classical algebra to respective operators maximum (symbolized by ⨁) and
times (symbolized by ⨂) in the domain of real numbers ℝ augmented by the number
minus infinity -∞. The (min,+) algebra (min-plus) redefines operators plus and
times from the classical algebra to respective operators minimum (symbolized by
⨁) and times (symbolized by ⨂) in the domain of real numbers ℝ augmented by the
number minus infinity +∞.

The interest in matrix computation in this algebra is taught since 1960 by
J. Kuntzman in his theory of networks. It is used in numerous domains such as
operational research (network theory), physics (quantization), probabilities
(Cramer's transform), law control (discrete events systems), computer science
(automata theory, Petri nets), and mathematics (algebraic geometry). This
algebra is also known as "tropical algebra".

This current Julia [MaxPlus.jl](https://github.com/Lecrapouille/MaxPlus.jl)
package is a portage in Julia language of the
[ScicosLab](http://www.scicoslab.org/)'s (max,+) toolbox which gave you
functions for doing numerical computations in (max, +) algebra. Due to the young
age of this toolbox, in case of doubt about obtained results, please compare them
with ScicosLab results, and if they are not matching, report an issue.

This Julia toolbox extends the original toolbox by adding the (min, +)
algebra. You may also be interested in this Timed Petri Net and Timed Event
Graphs graphical [editor](https://github.com/Lecrapouille/TimedPetriNetEditor)
which is also bound with (max, +) algebra. This editor can help you generate
(max, +) matrices from timed event graphs.

## Prerequisite

This MaxPlus.jl toolbox depends on the following Julia packages: `Printf,
PrettyTables, LinearAlgebra, SparseArrays`. They are installed automatically by
Julia's packager. The toolbox is supposed to work with any version of Julia >=
0.6.4 but a version >= 1.0 is the most recommended since older Julia versions
are no longer maintained. Depending on the version of your Julia some importants
issue in the core of Julia impacts this toolbox I had to add some fallbacks but
they may interfere with other packages you are using with MaxPlus.jl.

## Installation of the Julia (max,+) package MaxPlus.jl

There are different ways to install the package of this toolbox:

- Get the stable `MaxPlus.jl` version of the package from the official Julia
  packages. Type `]` then type `add MaxPlus`.

- Get the latest code source locally. From your Linux terminal type:
```
git clone https://github.com/Lecrapouille/MaxPlus.jl.git --depth=1
cd MaxPlus.jl
julia
```

Be sure to be inside the root of the git repository. Then, from the Julia REPL
type: `]` then type `add .` The API may be in gestation and not be stable and
changed from the available one in official Julia packages.

## Your first (max,+) lines of code in the REPL

Once the package has been installed, you have to activate the (max,+) package.
From the Julia REPL, type:

```julia
julia> using MaxPlus
```

Now, you can type your first lines of (max, +) code inside the Julia REPL:

```julia
julia> MP([1 2; 3 8]) .+ 5
```

Julia will reply to you:

```julia
2×2 (max,+) dense matrix:
  5   5
  5   8
```

Let's understand how Julia has made its computation:
- First, Julia creates a dense matrix of (max,+) numbers (typed `MP`): `[MP(1)
  MP(2); MP(3) MP(8)])`.
- For each element of the matrix (the `.` operator), the `⨁ 5` will be applied.
- Before the `⨁` is applied, the number `5` is converted to a (max,+)
  number `MP(5)`.
- Julia `Int64` and `Float64` numbers are implicitly promoted to a (max,+)
  number (internally encoded as `Float64`).

Symbols `⨁` and `⨂` are not used to avoid complex formulas hard to type and hard
to read, so keep using `+` and `*` symbols.

The equivalent of `MP([1 2; 3 8]) .+ 5` in classical algebra is:
```julia
  [max(1, 5)  max(2, 5)
   max(3, 5)  max(8, 5)]
```

and shall not be confused with this formula `[1 2; 3 8] .+ 5` in classical algebra
computing:
```julia
  [(1 + 5)  (2 + 5)
   (3 + 5)  (8 + 5)]
```

## Your first (min,+) lines of code in the REPL

This toolbox initially focused on the (max,+) algebra but one thing leading to another, functions for algebra (min,+) have been introduced but the name of the package MaxPlus.jl hasn't changed.

For (min,+) numbers:

```julia
julia> MI([1 2; 3 8])
2×2 (min,+) dense matrix:
  1   2
  3   8
```

## Documentation

Do you want to dive more about programming in (max,+) with this toolbox ? The following
documents are compiled into a single online documentation:
[https://lecrapouille.github.io/MaxPlus.jl](https://lecrapouille.github.io/MaxPlus.jl/index.html)
else within GitHub you can see them as Markdown:

* [Introduction and tutorials](tutorial) of this toolbox and tropical algebra in French and in English.
* The index of [(max,+) functions](docs/src/maxplus.md).
* The index of [(min,+) functions](docs/src/minplus.md).
* The index of [(max,+) linear system](docs/src/syslin.md).
* The [Rosetta Stone](docs/src/portage.md) to translate SicosLab to MaxPlus.jl functions.
* For developpers, you can run [non regression tests](docs/src/tests.md).
* A Timed Petri net and graph event [editor](https://github.com/Lecrapouille/TimedPetriNetEditor), a separate
  project of mine in relation to (max,+) algebra.
* The [bibliography](docs/src/bibliography.md) links to other documentation and ressources.

## Contributing

Feel free to contribute and particularly, since I'm more of a C++ guy, to help
reworking the code using more Julia formalism.
