### visualize a solution.


function visualizefiltersolution(𝓧, xq, wfunc, dfunc,f)
    𝑤_new = get𝑤scaled(𝓧)
    h_new = geth(𝑤_new, wfunc, dfunc, 𝓧)

    f_itp = xx->getp(xx, 𝑤_new, h_new, wfunc, dfunc, 𝓧)
    f_itp_xq = collect( f_itp(xq[i]) for i = 1:length(xq) )
    f_itp_𝓧 = collect( f_itp(𝓧[i]) for i = 1:length(𝓧) )

    𝑒 = xx->(f(xx)-f_itp(xx))
    𝑒_xq = collect( 𝑒(xq[i]) for i = 1:length(xq) )
    𝑒_𝓧 = collect( 𝑒(𝓧[i]) for i = 1:length(𝓧) )

    return f_itp_xq, f_itp_𝓧, 𝑒_xq, 𝑒_𝓧
end

function getfreqrsp(h::Vector{Float64}, resolution_multiple::Int = 20)
    N_samples = length(h)

    ω_set_fft = collect( LinRange(0,2*π-2*π/N_samples,N_samples))
    DFT_evals = fft(h)

    ω_set = collect( LinRange(0,2*π,N_samples*resolution_multiple) )
    DTFT_evals = collect( computeDTFTviaformula(h,ω_set[i]) for i = 1:length(ω_set) )

    return ω_set_fft, DFT_evals, ω_set, DTFT_evals
end

function plotmagnitudersp(h::Vector{Float64}, fig_num::Int, title_string::String = "Magnitude response")

    ω_set_fft, DFT_evals, ω_set, DTFT_evals = getfreqrsp(h)

    # visualize.
    mag_rsp_G = abs.(DTFT_evals)
    mag_rsp_fft = abs.(DFT_evals)


    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(ω_set_fft, mag_rsp_fft, ".", label = "DFT")
    PyPlot.plot(ω_set, mag_rsp_G, label = "DTFT")

    PyPlot.title(title_string)
    PyPlot.legend()

    return fig_num
end
