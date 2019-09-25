<div align="center"><img src="docs/logo/julia-max-plus.png" alt="Julia MaxPlus logo" width="150"/></div>

# Max-Plus Algebra

[![](https://travis-ci.org/Lecrapouille/MaxPlus.jl.svg?branch=master)](https://travis-ci.org/Lecrapouille/MaxPlus.jl)
[![](https://coveralls.io/repos/github/Lecrapouille/MaxPlus.jl/badge.svg?branch=master)](https://coveralls.io/github/Lecrapouille/MaxPlus.jl?branch=master)
[![](https://codecov.io/gh/Lecrapouille/MaxPlus.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Lecrapouille/MaxPlus.jl)

This repo contains a max-plus toolbox for Julia >= 1.0.3 for doing numerical computations with the tropical
semi-ring max,+ (ℝ ∪ {-∞}, ⊕, ⊙) where ⊕ is the usual multiplication and ⊙ is the usual maximum.

This repo is a portage of the [ScicosLab](http://www.scicoslab.org/) max-plus toolbox with the agreement of the original authors.
Note: due to the young age of this Julia toolbox, in case of doubt with the result of the function please compare it with ScicosLab. If
the result is not equal please report an issue on https://github.com/Lecrapouille/MaxPlus.jl

* [Functions:](docs/src/functions.md)
  description of Max-Plus functions.
* [Tutorial:](tutorial)
  introduces Max-Plus algebra via examples.
