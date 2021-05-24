![logo](https://lecrapouille.github.io/icons/juliamaxplus.png)

# Julia's Max-Plus Algebra Toolbox

[![](https://travis-ci.org/Lecrapouille/MaxPlus.jl.svg?branch=master)](https://travis-ci.org/Lecrapouille/MaxPlus.jl)
[![](https://coveralls.io/repos/github/Lecrapouille/MaxPlus.jl/badge.svg?branch=master)](https://coveralls.io/github/Lecrapouille/MaxPlus.jl?branch=master)
[![](https://codecov.io/gh/Lecrapouille/MaxPlus.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Lecrapouille/MaxPlus.jl)

The Max-Plus algebra redefines operators plus and times from the classical algebra to respectively operators maximum (symbolized by ⨁) and times (symbolized by ⨂) in the domain of real numbers ℝ augmented by the number minus infinity -∞.

The interest of matricial computation in this algebra is taugh since 1960 by J. Kuntzman in his theory of networks. It is used in numerous domains such as operational research (network theory), physics (quantization), probabilities (Cramer's transform), law control (discret events systems), computer science (automata theory, Petri nets), mathematics (algebraic geometry).

This [repo](https://github.com/Lecrapouille/MaxPlus.jl) is a portage of the [ScicosLab](http://www.scicoslab.org/) Max-Plus toolbox for Julia and gives you helper functions for doing numerical computations in (max, +) algebra. Note: due to the young age of this toolbox, in case of doubt with obtained results, please compare them with ScicosLab reults, and if they are not matching, report an issue.

## Prerequisite

This toolbox depends on the following Julia packages: `Printf, PrettyTables, LinearAlgebra,  SparseArrays`. They are installed from the Julia's packager. Max-Plus toolbox is supposed to work with any version of Julia >= 0.6.4 but a version >= 1.0 is the most recommended since older Julia versions are no longer maintained.

## Installation and your first Max-Plus code

Different ways:
- Get the stable version of the package from the Julia package manager: type `]` then type `add MaxPlus`.
- Get the latest code source localy: from your Linux terminal type `git clone https://github.com/Lecrapouille/MaxPlus.jl.git`, go inside the root of the project `cd MaxPlus.jl` then call `julia`. From the Julia REPL type: `push!(LOAD_PATH, pwd()); using MaxPlus`.

From the Julia REPL, you can type you first (max, +) code:
```julia
MP([1 2; 3 8]) .+ 5
```

This will, firstly, create a dense matrix of Max-Plus numbers (MP), and secondely, for each elements of it, the ⨁ 5 will be applied and therefore the result will the maximum between this element and the the Max-Plus number 5. Note that Julia Int64 and Float64 numbers are implictely promoted to a Max-Plus number. So the expected the result is :
```julia
2×2 Max-Plus dense matrix:
  5   5
  5   8
```

If this is the case, then your Max-Plus toolbox has been correctly installed and seem to work correctly. Optionnaly, to be totally sure, you can run unit-tests with the following Julia command:
```julia
] activate .
test
```

## Deeper dive with Julia's Max-Plus toolbox

* The index of Max-Plus functions available: [docs/src/functions.md](docs/src/functions.md)
* Introduction to this Max-Plus toolbox are given in the [tutorials](tutorials) folder. This is currently work in progress.

## Contribution

Feel free to contribute and particurlarly, since I'm more a C++ guy, to help reworking the code using more Julia formalism.
