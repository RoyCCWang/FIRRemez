
module FIRRemez

using FFTW

using Distributed
using SharedArrays

import Random
import LinearAlgebra
#import PyPlot
import Printf

import FiniteDiff

import DelimitedFiles


include("../src/misc/declarations.jl")
include("./minimax/Chebyshev.jl")
include("./minimax/interpolators.jl")
include("./minimax/eval.jl")
include("./minimax/filter.jl")
include("./minimax/exchange.jl")
include("./minimax/exchange_strategy.jl")
include("./minimax/engine.jl")
include("./minimax/extrema_search.jl")

include("../src/misc/utilities.jl")
include("../src/misc/template_signals.jl")
include("../src/misc/bounds.jl")
include("../src/misc/IO.jl")
include("../src/misc/frontend.jl")

include("./FIR/design.jl")
include("./FIR/type_conversion.jl")

include("../src/misc/visualize.jl")

export frontendoptimization, visualizeoptimzationsolution,
        naiveimpulsersp, getdefaultoptparameters,
        MiniMaxConfigType,
        Barycentric2nditp,
        frontendfilterdesign,
        getp

end
