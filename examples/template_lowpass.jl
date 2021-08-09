
using Distributed
using SharedArrays

@everywhere import FIRRemez
@everywhere import PyPlot
@everywhere import Random
@everywhere import Printf

import FileIO
import Printf

using FFTW
#import JLD
import XLSX
import BSON

fig_num = 1
PyPlot.close("all")

Random.seed!(25)

fig_num = 1


L = 25
#L = 1000
if length(ARGS) > 1
    L = tryparse(Int, ARGS[1])
end


## low-pass.
passband_ω_Float64 = 0.025*π # change this.
if length(ARGS) > 2
    passband_ω_Float64 = tryparse(Float64, ARGS[2])
end

stopband_ω_Float64 = 0.025*2*π # change this.
if length(ARGS) > 3
    stopband_ω_Float64 = tryparse(Float64, ARGS[3])
end

save_name_tag = "lowpass"
if length(ARGS) > 4
    save_name_tag = ARGS[4]
end

passband_ω = convert(BigFloat, passband_ω_Float64)
stopband_ω = convert(BigFloat, stopband_ω_Float64)

f = xx->FIRRemez.cosinetransitionlowpassfunction(xx,passband_ω,stopband_ω)


println("passband is ", passband_ω)
println("stopband is ", stopband_ω)

filter_type_num = 1 # symmetric, odd length.


savefig_flag = false
plot_output_flag = false

output_folder = "./output"
opt_params_file_name = "../default/config_opt.txt" # location of the configuration file.

# omit opt_params_file_name is using default configuration file provided by the package.
h_BigFloat, X_BigFloat, fig_num = FIRRemez.frontendfilterdesign(filter_type_num, f,
                            L, fig_num, savefig_flag, plot_output_flag,
                            opt_params_file_name)

h = convert(Vector{Float64},h_BigFloat)
X = convert(Vector{Float64}, X_BigFloat)

Printf.@printf("Order of Chebyshev is %d, length of filter is %d, FIR type %d\n", L, length(h), filter_type_num)

band_info_string = Printf.@sprintf("passband_%f_stopband_%f", passband_ω_Float64, stopband_ω_Float64)
save_name = Printf.@sprintf("type%d_%s_L%d_%s.bson",
    filter_type_num, save_name_tag, L, band_info_string)

passband = passband_ω_Float64
stopband = stopband_ω_Float64

BSON.bson(joinpath(output_folder, save_name),
    h = h, X = X, passband = passband, stopband = stopband)
#JLD.save(joinpath(output_folder, save_name), "h" = h, "X" = X, "passband" = passband, "stopband" = stopband)
#
save_name = Printf.@sprintf("type%d_%s_L%d_%s.xlsx",
    filter_type_num, save_name_tag, L, band_info_string)
XLSX.openxlsx(joinpath(output_folder, save_name), mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "filter")

    sheet["A1", dim = 1] = ["passband"; passband]
    sheet["B1", dim = 1] = ["passband"; passband]
    sheet["C1", dim = 1] = ["Chebyshev node positions"; X]
    sheet["D1", dim = 1] = ["coefficients"; h]
end

#fig_num = FIRRemez.plotmagnitudersp(h, fig_num, "filter's magnitude response")

import Plots
Plots.plotly()

ω_set_fft, DFT_evals, ω_set, DTFT_evals = FIRRemez.getfreqrsp(h)

mag_rsp_G = abs.(DTFT_evals)
mag_rsp_fft = abs.(DFT_evals)

plot_obj = Plots.plot( ω_set_fft,
mag_rsp_fft,
title = "Magnitude spectrum",
label = "DFT",
seriestype = :line,
ticks = :native,
hover = ω_set_fft,
linewidth = 4,
size = (1200, 900))

Plots.plot!(plot_obj, ω_set, mag_rsp_G, label = "DTFT",
seriestype = :line,
linestyle = :dash,
linewidth = 4)

save_name = Printf.@sprintf("type%d_%s_L%d_%s.html",
    filter_type_num, save_name_tag, L, band_info_string)
Plots.savefig(plot_obj, joinpath(output_folder, save_name))
