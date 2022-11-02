![logo](https://lecrapouille.github.io/icons/juliamaxplus.png)

# Julia's (max,+) and (min,+) Algebra Toolbox

[![](https://travis-ci.org/Lecrapouille/MaxPlus.jl.svg?branch=master)](https://travis-ci.org/Lecrapouille/MaxPlus.jl)
[![](https://coveralls.io/repos/github/Lecrapouille/MaxPlus.jl/badge.svg?branch=master)](https://coveralls.io/github/Lecrapouille/MaxPlus.jl?branch=master)
[![](https://codecov.io/gh/Lecrapouille/MaxPlus.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Lecrapouille/MaxPlus.jl)

The (max,+) algebra (or simply max-plus) redefines operators plus and times from
the classical algebra to respective operators maximum (symbolized by ⨁) and
times (symbolized by ⨂) in the domain of real numbers ℝ augmented by the number
minus infinity -∞. The (min,+) algebra (min-plus) redefines operators plus and
times from the classical algebra to respective operators minimum (symbolized
by ⨁) and times (symbolized by ⨂) in the domain of real numbers ℝ augmented by
the number minus infinity +∞.

The interest in matrix computation in this algebra is taught since 1960 by
J. Kuntzman in his theory of networks. It is used in numerous domains such as
operational research (network theory), physics (quantization), probabilities
(Cramer's transform), law control (discrete events systems), computer science
(automata theory, Petri nets), and mathematics (algebraic geometry). This algebra is
also known as "tropical algebra".

This [current Julia package](https://github.com/Lecrapouille/MaxPlus.jl) is a
portage in Julia language of the [ScicosLab](http://www.scicoslab.org/)'s (max,+)
toolbox which gave you functions for doing numerical computations in (max, +)
algebra. This Julia toolbox extends the original toolbox by adding the (min, +)
algebra. You may also be interested by this [Timed Petri Net
Editor](https://github.com/Lecrapouille/TimedPetriNetEditor) which is also
coupled by (max, +) algebra; this editor helps you generating (max, +) matrices
from timed event graphs.

Note: due to the young age of this toolbox, in case of doubt with obtained
results, please compare them with ScicosLab results, and if they are not
matching, report an issue.

## Prerequisite

This toolbox depends on the following Julia packages: `Printf, PrettyTables,
LinearAlgebra, SparseArrays`. They are installed automatically by Julia's
packager. (max,+) toolbox is supposed to work with any version of Julia >=
0.6.4 but a version >= 1.0 is the most recommended since older Julia versions
are no longer maintained.

## Installation of the Julia (max,+) package

Different ways to install the package of this toolbox:

- Get the stable `MaxPlus.jl` version of the package from the official Julia
  packages. Type `]` then type `add MaxPlus`. **Warning:** his API is deprecated
  and you have to follow the documentation of the [master
  branch](https://github.com/Lecrapouille/MaxPlus.jl/tree/master) of this
  repository. Waiting for an update, the next step is for the moment the better way:

- Get the latest code source locally. From your Linux terminal type:
```
git clone https://github.com/Lecrapouille/MaxPlus.jl.git -b dev --depth=1
cd MaxPlus.jl
julia
```

From the Julia REPL type: `]` then type `add /path/to/repository/MaxPlus.jl`
This new API is in gestation and not yet available in official Julia packages.

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

For (min,+) numbers;

```julia
julia> MI([1 2; 3 8])
2×2 (min,+) dense matrix:
  1   2
  3   8
```

## Your first (max,+) or (min,+) lines of code in Jupyter notebook

Julia REPL is fine to prototype but Jupyter notebook offers pretty prints.
This repository offers more detailed [tutorials](tutorial) using Jupyter notebook.
So let's use it.

Inside the Julia REPL:
```julia
using IJulia
notebook()
```

Inside the Jupyter notebook, type:
```julia
push!(LOAD_PATH, pwd())
using MaxPlus
```

Currently, to fix some conflict with Jupyter layout, you have to type these lines before making some computations:
```julia
Base.show(io::IO, ::MIME"text/latex", x::MP) = show(io, MIME"text/plain", x)
Base.show(io::IO, ::MIME"text/latex", A::MPAbstractVecOrMat) = show(io, MIME"text/plain", A)
```

You can type: `MP([1 2; 3 8]) .+ 5` and when pressing enter the answer will be printed.

## Unit tests

If Julia did not complain when installing this package and if you obtained a good result when typing `MP(5)`,
then your (max,+) toolbox has been correctly installed and seems to work correctly. To be totally sure, you can
run unit tests with the following Julia command:

```julia
] activate .
test
```

Hope, you will see:
```
Testing Running tests...
Testing MaxPlus tests passed
```

## Deeper dive with Julia's (max,+) toolbox

Do you want to know more about programming in (max,+)? The following documents are compiled into a single online documentation: [https://lecrapouille.github.io/MaxPlus.jl](https://lecrapouille.github.io/MaxPlus.jl/index.html)
* Introduction and tutorials are given in the [tutorials](tutorial) folder. There is currently work in progress in French.
* The index of (max,+) functions is available at [docs/src/functions.md](docs/src/functions.md)
* A Timed Petri net and graph event [editor](https://github.com/Lecrapouille/TimedPetriNetEditor), a separate
  project of mine in relation to (max,+) algebra.
* Links to other documentation [bibliography](docs/src/bibliography.md).

## Contributing

Feel free to contribute and particularly, since I'm more of a C++ guy, to help rework the code using more Julia formalism.

