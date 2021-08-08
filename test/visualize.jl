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


function packagefiltersolutionfast( 𝓧_new,
                                    h_history,
                                    wfunc,
                                    zfunc,
                                    dfunc,
                                    filter_type_num::Int)

##
    c,ν, discrepancy, xq, yq_bary, yq_cheby = BarycentrictoChebyshevfast(𝓧_new, wfunc, dfunc, filter_type_num)
    println("conversion discrepany between barycenter and Chebyshev is ", discrepancy)

    ######
    L = length(c)

    L_old = L

    # get impulse response.
    𝑏 = c
    if filter_type_num == 2
        𝑏 = converttotype2(c)
    elseif filter_type_num == 3
        𝑏 = converttotype3(c)
    elseif filter_type_num ==4
        𝑏 = converttotype4(c)
    end

    q = 1
    if filter_type_num == 1 || filter_type_num == 2
        q = 0
    end

    gethfunc = gethtype4
    if filter_type_num == 1
        L = length(c)-1

        gethfunc = gethtype1
    elseif filter_type_num == 2
        gethfunc = gethtype2
    elseif filter_type_num ==3
        gethfunc = gethtype3
    end

    h = gethfunc(𝑏)

    return h
end

function packagefiltersolution(𝓧_new,
                                h_history,
                                wfunc,
                                zfunc,
                                dfunc,
                                filter_type_num::Int,
                                fig_num::Int,
                                savefig_flag::Bool,
                                plot_flag::Bool)


    @time c,ν, discrepancy, xq, yq_bary, yq_cheby = BarycentrictoChebyshev(𝓧_new, wfunc, dfunc, filter_type_num)
    println("conversion discrepany between barycenter and Chebyshev is ", discrepancy)

    ######
    L = length(c)

    L_old = L

    # get impulse response.
    𝑏 = c
    if filter_type_num == 2
        𝑏 = converttotype2(c)
    elseif filter_type_num == 3
        𝑏 = converttotype3(c)
    elseif filter_type_num ==4
        𝑏 = converttotype4(c)
    end

    q = 1
    if filter_type_num == 1 || filter_type_num == 2
        q = 0
    end

    gethfunc = gethtype4
    if filter_type_num == 1
        L = length(c)-1

        gethfunc = gethtype1
    elseif filter_type_num == 2
        gethfunc = gethtype2
    elseif filter_type_num ==3
        gethfunc = gethtype3
    end


    h = gethfunc(𝑏)

    if !plot_flag
        return h, fig_num
    end

    ######### plot to interpret results.

    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(xq, yq_bary, "--", label = "yq_bary")
    PyPlot.plot(xq, yq_cheby, label = "yq_cheby")
    PyPlot.title("Barycenter vs. Chebyshev")
    PyPlot.legend()

    if savefig_flag
        save_name = Printf.@sprintf("./outputs/type_%s_Barycenter_to_Chebyshev.png",
                                        filter_type_num)
        PyPlot.savefig(save_name)
        sleep(save_delay)
        PyPlot.close(fig_num)
    end


    ## plot level error history.

    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(collect(1:length(h_history)), abs.(h_history), ".-")
    PyPlot.title("level error history")


    if savefig_flag
        save_name = Printf.@sprintf("./outputs/type_%s_level_error_history.png",
                                        filter_type_num)
        PyPlot.savefig(save_name)
        sleep(save_delay)
        PyPlot.close(fig_num)
    end

    ## plot resultant filter's response.
    f(xx) = wfunc(xx)*dfunc(xx)
    L = length(𝓧_new)-2
    xq = convert(Vector{BigFloat},collect(LinRange(-1,1,5*L)))

    f_itp_xq, f_itp_𝓧_new, 𝑒_xq, 𝑒_𝓧_new = visualizefiltersolution(𝓧_new,xq,wfunc,dfunc,f)


    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(xq, f.(xq), "--", label = "f = w*d")
    PyPlot.plot(xq, f_itp_xq, label = "solution")
    PyPlot.plot(xq, dfunc.(xq), label = "d")
    PyPlot.plot(𝓧_new, f_itp_𝓧_new,".", label = "reference positions")
    PyPlot.title("solution filter's response")
    PyPlot.legend()

    if savefig_flag
        save_name = Printf.@sprintf("./outputs/type_%s_solution.png",
                                        filter_type_num)
        PyPlot.savefig(save_name)
        sleep(save_delay)
        PyPlot.close(fig_num)
    end




    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(xq, 𝑒_xq, label = "error")
    PyPlot.plot(𝓧_new, 𝑒_𝓧_new,".", label = "reference positions")
    PyPlot.title("error")
    PyPlot.legend()

    if savefig_flag
        save_name = Printf.@sprintf("./outputs/type_%s_solution_error.png",
                                        filter_type_num)
        PyPlot.savefig(save_name)
        sleep(save_delay)
        PyPlot.close(fig_num)
    end


    ###### analyze

    # prepare.
    N = length(h)

    a = 0.0 #-π
    b = 2*π # π
    ω_set = convert(Vector{BigFloat}, collect( LinRange(a,b,N*20) ) )
    ω_set_fft = collect( LinRange(a,b-2*π/N,N))

    # this is the verified one.
    y = Vector{Complex{BigFloat}}(undef,0)
    y2 = Vector{BigFloat}(undef,0)

    if filter_type_num == 1
        y = collect( evalunivariateChebyshevpolynomial(c,cos(ω_set[i]))*exp( im*(q*π - (N-1)*ω_set[i])/2 ) for i = 1:length(ω_set) )
        y2 = collect( sum( c[k+1]*cos(k*ω_set[i]) for k = 0:L ) for i = 1:length(ω_set) )

    elseif filter_type_num == 2
        y = collect( cos(ω_set[i]/2)*evalunivariateChebyshevpolynomial(c,cos(ω_set[i]))*exp( im*(q*π - (N-1)*ω_set[i])/2 ) for i = 1:length(ω_set) )
        y2 = collect( sum( 𝑏[k]*cos((k-0.5)*ω_set[i]) for k = 1:L ) for i = 1:length(ω_set) )

    elseif filter_type_num == 3
       y = collect( sin(ω_set[i])*evalunivariateChebyshevpolynomial(c,cos(ω_set[i]))*exp( im*(q*π - (N-1)*ω_set[i])/2 ) for i = 1:length(ω_set) )
       y2 = collect( sum( 𝑏[k]*sin(k*ω_set[i]) for k = 1:L ) for i = 1:length(ω_set) )

    elseif filter_type_num == 4
       y = collect( sin(ω_set[i]/2)*evalunivariateChebyshevpolynomial(c,cos(ω_set[i]))*exp( im*(q*π - (N-1)*ω_set[i])/2 ) for i = 1:length(ω_set) )
       y2 = collect( sum( 𝑏[k]*sin((k-0.5)*ω_set[i]) for k = 1:L ) for i = 1:length(ω_set) )
    end


    y = abs.(y)
    y2 = abs.(y2)
    discrepancy = LinearAlgebra.norm(y-y2)
    println("conversion discrepany between y and y2 is ", discrepancy)


    ##### Visualize impulse response.
    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(h, label = "h")
    PyPlot.plot(h, ".")
    title_string = Printf.@sprintf("check if the recovered h has the characteristics of type %d", filter_type_num)
    PyPlot.title(title_string)
    PyPlot.legend()

    if savefig_flag
        save_name = Printf.@sprintf("./outputs/type_%s_solution_impulse_rsp.png",
                                        filter_type_num)
        PyPlot.savefig(save_name)
        sleep(save_delay)
        PyPlot.close(fig_num)
    end

    ##### prepare fft and DFT from impulse response.

    impulse_rsp = fft(convert(Vector{Float64},h))
    mag_rsp_fft = abs.(impulse_rsp)
    phase_rsp_fft = angle.(impulse_rsp)

    impulse_rsp_G = collect( naiveimpulsersp(h,ω_set[i]) for i = 1:length(ω_set) )
    mag_rsp_G = abs.(impulse_rsp_G)
    phase_rsp_G = angle.(impulse_rsp_G)

    println("discrepancy between y and G is ", LinearAlgebra.norm(y-mag_rsp_G))

    ##### visualize.

    ## plot the polynomial specified by b, see if it is in agreement with
    #   the fft/DFT response.
    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(ω_set_fft, mag_rsp_fft, label = "fft: solution discrete")
    PyPlot.plot(ω_set_fft, mag_rsp_fft, ".")

    PyPlot.plot(ω_set, mag_rsp_G, label = "G: solution continuous")

    y_d = abs.(wfunc.(cos.(ω_set)).*dfunc.(cos.(ω_set)))
    PyPlot.plot(ω_set, y_d, label = "y_d: target")

    # # same as y_d but without the abs().
    # zq = convert(Vector{BigFloat}, collect( LinRange(0,π,5*L) ) )
    # z = wfunc.(cos.(zq)).*zfunc.(zq)
    # PyPlot.plot(zq, z, "--", label = "(w∘cos).*targetfunc")

    PyPlot.title("Magnitude response")
    PyPlot.legend()

    if savefig_flag
        save_name = Printf.@sprintf("./outputs/type_%s_solution_magnitude_rsp.png",
                                        filter_type_num)
        PyPlot.savefig(save_name)
        sleep(save_delay)
        PyPlot.close(fig_num)
    end

    return h, fig_num
end
