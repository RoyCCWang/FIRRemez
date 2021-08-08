
using Distributed

@everywhere import FIRRemez
@everywhere import PyPlot
@everywhere import Random
@everywhere import Printf

import FileIO
import Printf

using FFTW
import BSON


fig_num = 1
PyPlot.close("all")

Random.seed!(25)

fig_num = 1

L = 50 # change this.

## low-pass.
rp_ω_Float64 = 0.05*π # change this.
rp_ω = convert(BigFloat, rp_ω_Float64)

rs_ω_Float64 = 0.2*π # change this.
rs_ω = convert(BigFloat, rs_ω_Float64)

fp_ω_Float64 = 0.8*π # change this.
fp_ω = convert(BigFloat, fp_ω_Float64)

fs_ω_Float64 = 0.95*π # change this.
fs_ω = convert(BigFloat, fs_ω_Float64)

#f = xx->FIRRemez.cosinetransitionlowpassfunction(xx,passband_ω,stopband_ω)

@everywhere function cosinetransitionbandpassfunction(y::T, rs, rp, fp, fs)::T where T
    abs_y = abs(y)

    # rising edge.
    if abs_y <= rs
        return zero(T)
    elseif rs < abs_y <= rp
        #return cos( (rp-y)/(rp-rs)*π/2 )
        return cos( (rp-abs_y)/(rp-rs)*π/2 )
    end

    # falling edge.
    if rp <= abs_y <= fp
        return one(T)
    elseif fp < abs_y <= fs
        #return cos( (fp-y)/(fs-fp)*π/2 )
        return cos( (fp-abs_y)/(fs-fp)*π/2 )
    end

    return zero(T)
end

f = xx->cosinetransitionbandpassfunction(xx, rp_ω,
                                    rs_ω, fp_ω, fs_ω)
save_name_tag = "bandpass"

N_viz = 500
ω_range = LinRange(0.0, π, N_viz)
f_ω = f.(ω_range)
PyPlot.plot(ω_range, f_ω)


println("rising passband is ", rp_ω_Float64)
println("rising stopband is ", rs_ω_Float64)
println("falling passband is ", fp_ω_Float64)
println("falling stopband is ", fs_ω_Float64)

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

band_info_string = Printf.@sprintf("rp_%f_rs_%f_fp_%f_fs_%f", rp_ω_Float64,
                            rs_ω_Float64, fp_ω_Float64, fs_ω_Float64)

save_name = Printf.@sprintf("bp_type%d_%s_L%d_%s.bson", filter_type_num,
                                                    save_name_tag,
                                                    L,
                                                    band_info_string)
# rp = rp_ω_Float64
# rs = rs_ω_Float64
# passband = rp_ω_Float64
# stopband = rs_ω_Float64

BSON.bson(save_name,  filter_coeffs = h,
                        Chebyshev_nodes = X,
                        rising_passband_angular = rp_ω_Float64,
                        rising_stopband_angular = rs_ω_Float64,
                        falling_passband_angular = fp_ω_Float64,
                        falling_stopband_angular = fs_ω_Float64)

import VisualizationTools
fig_num = FIRRemez.plotmagnitudersp(h, fig_num, "filter's magnitude response")
