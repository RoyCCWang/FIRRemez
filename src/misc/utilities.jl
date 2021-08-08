
# first sorts x, then forces the entries in x such that f(x) has alternating sign.
function forcealternating(x_in::Vector{T}, f::Function) where T
    @assert length(x_in) > 1

    x = sort(x_in)

    i = 2
    while i <= length(x)

        fx_im1 = f(x[i-1])
        u::Bool = sign(fx_im1) > 0

        fx_i = f(x[i])
        v = sign(fx_i) > 0

        if !xor(u,v)
            # same sign for error evaluations of x[i-1] and x[i].
            # Remove the position with smaller absolute error.
            if abs(fx_i) > abs(fx_im1)
                deleteat!(x,i-1)

            else
                deleteat!(x,i)
            end
        else
            i += 1
        end
    end

    # the above loop sometimes miss the end points of the position.
    # check if endpoints have different error signs. If so, add them back.
    if xor( sign(f(x_in[1])) > 0, sign(f(x[1])) > 0 )
        insert!(x, 1, x_in[1])
    end

    if xor( sign(f(x_in[end])) > 0, sign(f(x[end])) > 0 )
        push!(x, x_in[end])
    end

    return x
end




# try ForwardDiff.derivative(f,x) someday.
function numericalchkextrema(f::Function, x::T, tol)::Bool where T
    return numericalchkextrema(f,x,convert(T,tol))
end

# TODO try passing df instead of f.
function numericalchkextrema(f::Function, x::T, tol::T)::Bool where T

    if isnumericallyclose(FiniteDiff.finite_difference_derivative(f,x), zero(T), tol)
        #δ = 1e-3
        df_δ = FiniteDiff.finite_difference_derivative(f,x+tol)
        df_mδ = FiniteDiff.finite_difference_derivative(f,x-tol)
        #return true
        return xor(df_δ > zero(T),df_mδ > zero(T))
    end

    return false
end



function returnNclosestzeroderivatives(a::Vector{T}, 𝑒::Function, N::Int, tol::T) where T
    @assert length(a) >= N

    abs_𝑑out_𝑑x = collect( abs(FiniteDiff.finite_difference_derivative(𝑒, a[i])) for i = 1:length(a) )

    sort_ind = sortperm(abs_𝑑out_𝑑x)

    out = a[sort_ind[1:N]]

    if out[end] > tol

        # if the smallest zero is larger than the input tolerance, then just return the original array unmodified.
        return a, false
    end

    return out, true
end



# array version of samplepoints(f::Function,...). F is f(k,L), k = 0,1,...,L
function samplepoints(F::Vector{T}, a::T, b::T)::Vector{T} where T
    @assert a < b

    two = one(T) + one(T)
    Δ = b-a
    return collect( (F[i]+one(T))/two*Δ+a for i = 1:length(F) )
end

# f generates N points between the interval [-1,1].
function samplepoints(f::Function, a::T, b::T, N::Int)::Vector{T} where T
    @assert a < b
    L = N-1

    two = one(T) + one(T)
    Δ = b-a
    return collect( (f(k,L)+one(T))/two*Δ+a for k = 0:L )
end


function packagefiltersolutionfast( 𝓧_new,
                                    h_history,
                                    wfunc,
                                    zfunc,
                                    dfunc,
                                    filter_type_num::Int)


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

    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(xq, yq_bary, "--", label = "yq_bary")
    # PyPlot.plot(xq, yq_cheby, label = "yq_cheby")
    # PyPlot.title("Barycenter vs. Chebyshev")
    # PyPlot.legend()
    #
    # if savefig_flag
    #     save_name = Printf.@sprintf("./outputs/type_%s_Barycenter_to_Chebyshev.png",
    #                                     filter_type_num)
    #     PyPlot.savefig(save_name)
    #     sleep(save_delay)
    #     PyPlot.close(fig_num)
    # end


    ## plot level error history.

    # PyPlot.figure(fig_num)
    # fig_num += 1
    #
    # PyPlot.plot(collect(1:length(h_history)), abs.(h_history), ".-")
    # PyPlot.title("level error history")
    #
    #
    # if savefig_flag
    #     save_name = Printf.@sprintf("./outputs/type_%s_level_error_history.png",
    #                                     filter_type_num)
    #     PyPlot.savefig(save_name)
    #     sleep(save_delay)
    #     PyPlot.close(fig_num)
    # end

    ## plot resultant filter's response.
    f(xx) = wfunc(xx)*dfunc(xx)
    L = length(𝓧_new)-2
    xq = convert(Vector{BigFloat},collect(LinRange(-1,1,5*L)))

    f_itp_xq, f_itp_𝓧_new, 𝑒_xq, 𝑒_𝓧_new = visualizefiltersolution(𝓧_new,xq,wfunc,dfunc,f)


    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(xq, f.(xq), "--", label = "f = w*d")
    # PyPlot.plot(xq, f_itp_xq, label = "solution")
    # PyPlot.plot(xq, dfunc.(xq), label = "d")
    # PyPlot.plot(𝓧_new, f_itp_𝓧_new,".", label = "reference positions")
    # PyPlot.title("solution filter's response")
    # PyPlot.legend()
    #
    # if savefig_flag
    #     save_name = Printf.@sprintf("./outputs/type_%s_solution.png",
    #                                     filter_type_num)
    #     PyPlot.savefig(save_name)
    #     sleep(save_delay)
    #     PyPlot.close(fig_num)
    # end




    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(xq, 𝑒_xq, label = "error")
    # PyPlot.plot(𝓧_new, 𝑒_𝓧_new,".", label = "reference positions")
    # PyPlot.title("error")
    # PyPlot.legend()
    #
    # if savefig_flag
    #     save_name = Printf.@sprintf("./outputs/type_%s_solution_error.png",
    #                                     filter_type_num)
    #     PyPlot.savefig(save_name)
    #     sleep(save_delay)
    #     PyPlot.close(fig_num)
    # end


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
    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(h, label = "h")
    # PyPlot.plot(h, ".")
    # title_string = Printf.@sprintf("check if the recovered h has the characteristics of type %d", filter_type_num)
    # PyPlot.title(title_string)
    # PyPlot.legend()
    #
    # if savefig_flag
    #     save_name = Printf.@sprintf("./outputs/type_%s_solution_impulse_rsp.png",
    #                                     filter_type_num)
    #     PyPlot.savefig(save_name)
    #     sleep(save_delay)
    #     PyPlot.close(fig_num)
    # end

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
    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(ω_set_fft, mag_rsp_fft, label = "fft: solution discrete")
    # PyPlot.plot(ω_set_fft, mag_rsp_fft, ".")
    #
    # PyPlot.plot(ω_set, mag_rsp_G, label = "G: solution continuous")
    #
    # y_d = abs.(wfunc.(cos.(ω_set)).*dfunc.(cos.(ω_set)))
    # PyPlot.plot(ω_set, y_d, label = "y_d: target")
    #
    # PyPlot.title("Magnitude response")
    # PyPlot.legend()
    #
    # if savefig_flag
    #     save_name = Printf.@sprintf("./outputs/type_%s_solution_magnitude_rsp.png",
    #                                     filter_type_num)
    #     PyPlot.savefig(save_name)
    #     sleep(save_delay)
    #     PyPlot.close(fig_num)
    # end

    return h, fig_num
end
