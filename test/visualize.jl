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


function packagefiltersolutionfast( π§_new,
                                    h_history,
                                    wfunc,
                                    zfunc,
                                    dfunc,
                                    filter_type_num::Int)

##
    c,Ξ½, discrepancy, xq, yq_bary, yq_cheby = BarycentrictoChebyshevfast(π§_new, wfunc, dfunc, filter_type_num)
    println("conversion discrepany between barycenter and Chebyshev is ", discrepancy)

    ######
    L = length(c)

    L_old = L

    # get impulse response.
    π = c
    if filter_type_num == 2
        π = converttotype2(c)
    elseif filter_type_num == 3
        π = converttotype3(c)
    elseif filter_type_num ==4
        π = converttotype4(c)
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

    h = gethfunc(π)

    return h
end

function packagefiltersolution(π§_new,
                                h_history,
                                wfunc,
                                zfunc,
                                dfunc,
                                filter_type_num::Int,
                                fig_num::Int,
                                savefig_flag::Bool,
                                plot_flag::Bool)


    @time c,Ξ½, discrepancy, xq, yq_bary, yq_cheby = BarycentrictoChebyshev(π§_new, wfunc, dfunc, filter_type_num)
    println("conversion discrepany between barycenter and Chebyshev is ", discrepancy)

    ######
    L = length(c)

    L_old = L

    # get impulse response.
    π = c
    if filter_type_num == 2
        π = converttotype2(c)
    elseif filter_type_num == 3
        π = converttotype3(c)
    elseif filter_type_num ==4
        π = converttotype4(c)
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


    h = gethfunc(π)

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
    L = length(π§_new)-2
    xq = convert(Vector{BigFloat},collect(LinRange(-1,1,5*L)))

    f_itp_xq, f_itp_π§_new, π_xq, π_π§_new = visualizefiltersolution(π§_new,xq,wfunc,dfunc,f)


    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(xq, f.(xq), "--", label = "f = w*d")
    PyPlot.plot(xq, f_itp_xq, label = "solution")
    PyPlot.plot(xq, dfunc.(xq), label = "d")
    PyPlot.plot(π§_new, f_itp_π§_new,".", label = "reference positions")
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
    PyPlot.plot(xq, π_xq, label = "error")
    PyPlot.plot(π§_new, π_π§_new,".", label = "reference positions")
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

    a = 0.0 #-Ο
    b = 2*Ο # Ο
    Ο_set = convert(Vector{BigFloat}, collect( LinRange(a,b,N*20) ) )
    Ο_set_fft = collect( LinRange(a,b-2*Ο/N,N))

    # this is the verified one.
    y = Vector{Complex{BigFloat}}(undef,0)
    y2 = Vector{BigFloat}(undef,0)

    if filter_type_num == 1
        y = collect( evalunivariateChebyshevpolynomial(c,cos(Ο_set[i]))*exp( im*(q*Ο - (N-1)*Ο_set[i])/2 ) for i = 1:length(Ο_set) )
        y2 = collect( sum( c[k+1]*cos(k*Ο_set[i]) for k = 0:L ) for i = 1:length(Ο_set) )

    elseif filter_type_num == 2
        y = collect( cos(Ο_set[i]/2)*evalunivariateChebyshevpolynomial(c,cos(Ο_set[i]))*exp( im*(q*Ο - (N-1)*Ο_set[i])/2 ) for i = 1:length(Ο_set) )
        y2 = collect( sum( π[k]*cos((k-0.5)*Ο_set[i]) for k = 1:L ) for i = 1:length(Ο_set) )

    elseif filter_type_num == 3
       y = collect( sin(Ο_set[i])*evalunivariateChebyshevpolynomial(c,cos(Ο_set[i]))*exp( im*(q*Ο - (N-1)*Ο_set[i])/2 ) for i = 1:length(Ο_set) )
       y2 = collect( sum( π[k]*sin(k*Ο_set[i]) for k = 1:L ) for i = 1:length(Ο_set) )

    elseif filter_type_num == 4
       y = collect( sin(Ο_set[i]/2)*evalunivariateChebyshevpolynomial(c,cos(Ο_set[i]))*exp( im*(q*Ο - (N-1)*Ο_set[i])/2 ) for i = 1:length(Ο_set) )
       y2 = collect( sum( π[k]*sin((k-0.5)*Ο_set[i]) for k = 1:L ) for i = 1:length(Ο_set) )
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

    impulse_rsp_G = collect( naiveimpulsersp(h,Ο_set[i]) for i = 1:length(Ο_set) )
    mag_rsp_G = abs.(impulse_rsp_G)
    phase_rsp_G = angle.(impulse_rsp_G)

    println("discrepancy between y and G is ", LinearAlgebra.norm(y-mag_rsp_G))

    ##### visualize.

    ## plot the polynomial specified by b, see if it is in agreement with
    #   the fft/DFT response.
    PyPlot.figure(fig_num)
    fig_num += 1
    PyPlot.plot(Ο_set_fft, mag_rsp_fft, label = "fft: solution discrete")
    PyPlot.plot(Ο_set_fft, mag_rsp_fft, ".")

    PyPlot.plot(Ο_set, mag_rsp_G, label = "G: solution continuous")

    y_d = abs.(wfunc.(cos.(Ο_set)).*dfunc.(cos.(Ο_set)))
    PyPlot.plot(Ο_set, y_d, label = "y_d: target")

    # # same as y_d but without the abs().
    # zq = convert(Vector{BigFloat}, collect( LinRange(0,Ο,5*L) ) )
    # z = wfunc.(cos.(zq)).*zfunc.(zq)
    # PyPlot.plot(zq, z, "--", label = "(wβcos).*targetfunc")

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
