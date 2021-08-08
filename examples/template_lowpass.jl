
@everywhere import FIRRemez
@everywhere import PyPlot
@everywhere import Random
@everywhere import Printf

import FileIO
import Printf

using FFTW
import JLD


fig_num = 1
PyPlot.close("all")

Random.seed!(25)

fig_num = 1

#L = 25 # change this.
L = 1000

## low-pass.
passband_ω_Float64 = 0.025*π # change this.
passband_ω = convert(BigFloat, passband_ω_Float64)

stopband_ω_Float64 = 0.025*2*π # change this.
stopband_ω = convert(BigFloat, stopband_ω_Float64)

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

fig_num = VisualizationTools.plotmagnitudersp(h, fig_num, "filter's magnitude response")
