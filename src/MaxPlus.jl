# ==============================================================================
# Max-Plus Algebra toolbox for Julia >= 1.0.3
# A portage of the ScicosLab Max-Plus toolbox http://www.scicoslab.org/
# License: public domain
#
# Note: the documentation of functions for the REPL are placed in docstrings.jl
# ==============================================================================

module MaxPlus

include("types.jl")
include("fallbacks.jl")
include("show.jl")
include("core.jl")
include("howard.jl")
include("syslin.jl")
#include("flowshop.jl")
include("docstrings.jl")

end
