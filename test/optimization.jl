# mini-max fit to a generic function on the interval [-1,1].
# Ignores the domain outside of this interval.

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


fig_num = 1
PyPlot.close("all")

Random.seed!(25)

fig_num = 1

L = 25
f = xx->sin(xx*π+3*xx)+xx

𝓧_new = frontendoptimization(f, L)

xq = convert(Vector{BigFloat},collect(LinRange(-2,1,10*L)))
f_itp_xq, f_itp_𝓧, 𝑒_xq, 𝑒_𝓧, fig_num = visualizeoptimzationsolution(xq, 𝓧_new, f, fig_num)



@assert 1==232



L = 25 # 55 # set to odd.
# L_interpolant = L
# if isodd(L)
#     L_interpolant += 1
# end
L_interpolant = 4

max_iters = 200
filter_type_num = 2
tol_converged = convert(BigFloat,1e-3)
tol_𝓧_spacing = convert(BigFloat,1e-6/(L+2))
tol_derivative_zero = convert(BigFloat,0.2)
tol_no_update = convert(BigFloat,1e-6)
tol_h = convert(BigFloat,1e-17)

N_samples_interval = 10
N_test_positions = 200

verbose_flag = true
plot_flag = false #true

savefig_flag = true
save_delay = 0.5

config = MiniMaxConfigType( tol_converged,
                            tol_𝓧_spacing,
                            tol_derivative_zero,
                            tol_no_update,
                            tol_h,
                            N_samples_interval,
                            N_test_positions,
                            max_iters,
                            L_interpolant,
                            verbose_flag,
                            plot_flag)


candidate_multiple = 20

## set up functions.
wfunc = getweightfunc(filter_type_num)


zfunc = xx->sin(xx*π+3*xx)+xx
dfunc = xx->zfunc(acos(xx))
#dfunc = xx->sin(xx*π+3*xx)+xx

#dfunc = xx->rectfunc(xx, BigFloat("0.5"))
#dfunc = xx->sinc(xx)

# β = BigFloat("1")
# 𝑇 = BigFloat("2")
# dfunc = xx->raisedcosinefreqrsp(xx,β,𝑇)

𝓧_selection_string = "Chebyshev"

@time 𝓧_new, h_history,
        iters_ran, status_msg, fig_num = designfilter(  config,
                                                        L,
                                                        candidate_multiple,
                                                        wfunc,
                                                        dfunc,
                                                        𝓧_selection_string,
                                                        fig_num)
#

println("Status: ", status_msg)

include("visualize.jl")
