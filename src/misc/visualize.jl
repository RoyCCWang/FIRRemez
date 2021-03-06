### visualize a solution.


function visualizefiltersolution(π§, xq, wfunc, dfunc,f)
    π€_new = getπ€scaled(π§)
    h_new = geth(π€_new, wfunc, dfunc, π§)

    f_itp = xx->getp(xx, π€_new, h_new, wfunc, dfunc, π§)
    f_itp_xq = collect( f_itp(xq[i]) for i = 1:length(xq) )
    f_itp_π§ = collect( f_itp(π§[i]) for i = 1:length(π§) )

    π = xx->(f(xx)-f_itp(xx))
    π_xq = collect( π(xq[i]) for i = 1:length(xq) )
    π_π§ = collect( π(π§[i]) for i = 1:length(π§) )

    return f_itp_xq, f_itp_π§, π_xq, π_π§
end

function computeDTFTviaformula(h::AbstractArray{T}, Ο::T)::Complex{T} where T <: Real
    N = length(h)
    return sum( h[n+1]*exp(-im*Ο*n) for n = 0:N-1 )
end

function computeDTFTviaformula(h::AbstractArray{T}, Ο::Vector{T}) where T

    return collect( computeDTFTviaformula(h, Ο[i]) for i = 1:length(Ο))
end

function getfreqrsp(h::Vector{Float64}, resolution_multiple::Int = 20)
    N_samples = length(h)

    Ο_set_fft = collect( LinRange(0,2*Ο-2*Ο/N_samples,N_samples))
    DFT_evals = fft(h)

    Ο_set = collect( LinRange(0,2*Ο,N_samples*resolution_multiple) )
    DTFT_evals = collect( computeDTFTviaformula(h, Ο_set[i]) for i = 1:length(Ο_set) )

    return Ο_set_fft, DFT_evals, Ο_set, DTFT_evals
end

function plotmagnitudersp(h::Vector{Float64}, fig_num::Int, title_string::String = "Magnitude response")

    Ο_set_fft, DFT_evals, Ο_set, DTFT_evals = getfreqrsp(h)

    # visualize.
    mag_rsp_G = abs.(DTFT_evals)
    mag_rsp_fft = abs.(DFT_evals)


    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(Ο_set_fft, mag_rsp_fft, ".", label = "DFT")
    PyPlot.plot(Ο_set, mag_rsp_G, label = "DTFT")

    PyPlot.title(title_string)
    PyPlot.legend()

    return fig_num
end
