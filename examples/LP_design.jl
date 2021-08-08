
using Distributed
using SharedArrays

@everywhere import FIRRemez
@everywhere import PyPlot
@everywhere import Random
@everywhere import Printf

import FileIO
import Printf

import VisualizationTools

using FFTW
import JLD


fig_num = 1
PyPlot.close("all")

Random.seed!(25)

fig_num = 1

L = 1000 # change this.

# try 1
#𝑟0 = 0.00625*2 # relative central transition band freq, in rads, but /π.

𝑟0 = 0.00625*3 # try 3

#κ = 𝑟0*0.2
κ = 𝑟0*0.5 # try 2, 3
passband_ω = BigFloat(𝑟0*π-κ) # stopband_ω*0.9
stopband_ω = BigFloat(𝑟0*π+κ)
passband_ω_Float64 = convert(Float64, passband_ω)
stopband_ω_Float64 = convert(Float64, stopband_ω)

f = xx->FIRRemez.cosinetransitionlowpassfunction(xx,passband_ω,stopband_ω)
save_name_tag = "lowpass"

println("passband is ", passband_ω)
println("stopband is ", stopband_ω)

filter_type_num = 1 # symmetric, odd length.

savefig_flag = false
plot_output_flag = false

opt_params_file_name = "../default/config_opt.txt" # location of the configuration file.

# omit opt_params_file_name is using default configuration file provided by the package.
h_BigFloat, X_BigFloat, fig_num = FIRRemez.frontendfilterdesign(filter_type_num, f,
                            L, fig_num, savefig_flag, plot_output_flag,
                            opt_params_file_name)

h = convert(Vector{Float64},h_BigFloat)
X = convert(Vector{Float64}, X_BigFloat)

Printf.@printf("Order of Chebyshev is %d, length of filter is %d, FIR type %d\n", L, length(h), filter_type_num)

band_info_string = Printf.@sprintf("passband_%f_stopband_%f", passband_ω_Float64, stopband_ω_Float64)
save_name = Printf.@sprintf("type%d_%s_L%d_%s.jld", filter_type_num,
                                                    save_name_tag,
                                                    L,
                                                    band_info_string)
passband = passband_ω_Float64
stopband = stopband_ω_Float64
FileIO.save(save_name,  "h", h,
                        "X", X,
                        "passband", passband,
                        "stopband", stopband)

println("Plotting")
@time fig_num = VisualizationTools.plotmagnitudersp(h, fig_num, "filter's magnitude response")

# this was the output on threadripper.
# Reference set no longer changing.
# 6738.074528 seconds (1.68 G allocations: 80.051 GiB, 0.43% gc time)
# Completed mini-max algorithm. Status: Reference set no longer changing.
# Start to package filter solution.
# 1561.200300 seconds (32.29 G allocations: 1.643 TiB, 20.47% gc time)
# conversion discrepany between barycenter and Chebyshev is 8.833208284570889182364488789150169650254777059595584072429122220208126304367404e-19
# 1561.382912 seconds (32.29 G allocations: 1.644 TiB, 20.47% gc time)
# Finished packaging filter solution.
# Order of Chebyshev is 1000, length of filter is 2003, FIR type 1
# Plotting
#   1.885682 seconds (80.37 k allocations: 3.159 MiB)
