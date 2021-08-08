
using Profile

using FFTW

using Distributed
using SharedArrays

@everywhere import DelimitedFiles

@everywhere import Random
@everywhere import LinearAlgebra
@everywhere import PyPlot
@everywhere import Printf

@everywhere import Calculus
@everywhere import ForwardDiff

@everywhere include("../src/misc/declarations.jl")
@everywhere include("../src/minimax/Chebyshev.jl")
#@everywhere include("../src/minimax/APF.jl")
@everywhere include("../src/minimax/interpolators.jl")
@everywhere include("../src/minimax/eval.jl")
@everywhere include("../src/minimax/filter.jl")

@everywhere include("../src/misc/utilities.jl")
@everywhere include("../src/misc/template_signals.jl")
@everywhere include("../src/misc/bounds.jl")
@everywhere include("../src/minimax/exchange.jl")
@everywhere include("../src/minimax/exchange_strategy.jl")
@everywhere include("../src/misc/IO.jl")

@everywhere include("../src/minimax/engine.jl")
@everywhere include("../src/minimax/extrema_search.jl")

@everywhere include("../src/FIR/design.jl")
@everywhere include("../src/FIR/type_conversion.jl")

@everywhere include("../src/misc/frontend.jl")
@everywhere include("../test/visualize.jl")
@everywhere include("../test/filter_helpers.jl")


ig_num = 1
PyPlot.close("all")

Random.seed!(25)

fig_num = 1

#L = 800
L = 25

# bandpass.
# rs = 1.05664
# rp = 1.45664
# fp = 1.68496
# fs = 2.08496
# f = xx->cosinetransitionbandpassfunction(xx, BigFloat(rs),
#                                             BigFloat(rp),
#                                             BigFloat(fp),
#                                             BigFloat(fs))

# rs = 0.428319
# rp = 0.828319
# fp = 1.05664
# fs = 1.45664
# f = xx->cosinetransitionbandpassfunction(xx, BigFloat(rs), BigFloat(rp),
#                                              BigFloat(fp), BigFloat(fs))

# lowpass.
fp = 0.428319
fs = 0.828319
f = xx->cosinetransitionlowpassfunction(xx, BigFloat(fp), BigFloat(fs))

filter_type_num = 1 # symmetric, odd length.
max_iters = 50
verbose_flag = true
debug_flag = false # plot reference update at every iteration.

savefig_flag = false
plot_output_flag = true # plot mag rsp.

profile_flag = false

#opt_params_file_name = "0" # no plot between iterations.
#opt_params_file_name = "/home/roy/Documents/repos/FIRRemez/test/config_opt.txt" # plot between iterations.


h_BigFloat = Vector{BigFloat}(undef,0)
h = Vector{Float64}(undef,0)
X = Vector{BigFloat}(undef,0)

if profile_flag
    io = open("./outputs/profile_output.txt", "w")

    Profile.clear()
    Profile.init(n=15_000)
    #Profile.init()

    @profile frontendfilterdesign(filter_type_num, f,
                                L, fig_num, savefig_flag, plot_output_flag,
                                max_iters, verbose_flag, debug_flag)
    Profile.print(io);

    # using ProfileView
    # ProfileView.view()

else

    config_opt = getdefaultoptparameters(L)
    config_opt.max_iters = max_iters
    config_opt.verbose_flag = verbose_flag
    config_opt.plot_flag = debug_flag

    h_BigFloat, X, fig_num = frontendfilterdesignfast(filter_type_num, f,
                                L, fig_num, false, false,
                                config_opt)
    h = convert(Vector{Float64}, h_BigFloat)

    Printf.@printf("Order of Chebyshev is %d, length of filter is %d, FIR type %d\n", L, length(h), filter_type_num)
end
