
function frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                savefig_flag, plot_output_flag::Bool,
                                max_iters::Int, verbose_flag::Bool,
                                debug_flag::Bool)

    config = getdefaultoptparameters(L)
    config.max_iters = max_iters
    config.verbose_flag = verbose_flag
    config.plot_flag = debug_flag

    return frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                    savefig_flag, plot_output_flag,
                                    config)

end

function frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                savefig_flag, plot_output_flag::Bool)

    config = getdefaultoptparameters(L)

    return frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                    savefig_flag, plot_output_flag,
                                    config)

end

function frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                savefig_flag, plot_output_flag::Bool,
                                opt_parameters_file_name::String)

    config = getoptparameters(L,opt_parameters_file_name)

    return frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                    savefig_flag, plot_output_flag,
                                    config)
end

function frontendfilterdesign(filter_type_num, zfunc, L, fig_num,
                                savefig_flag, plot_output_flag,
                                config::MiniMaxConfigType)

    𝓧 = collect( Chebyshev2ndnode(k,L+1) for k = 0:L+1 )

    ##
    wfunc = getweightfunc(filter_type_num)
    dfunc = xx->zfunc(acos(xx))

    𝓧_selection_string = "Chebyshev"
    candidate_multiple = 1

    @time 𝓧_new, h_history,
            iters_ran, status_msg, fig_num = designfilter(  config,
                                                            L,
                                                            candidate_multiple,
                                                            wfunc,
                                                            dfunc,
                                                            𝓧_selection_string,
                                                            fig_num)

    #
    println("Completed mini-max algorithm. Status: ", status_msg)

    println("Start to package filter solution.")
    @time h, fig_num = packagefiltersolution(𝓧_new, h_history, wfunc, zfunc, dfunc,
                                        filter_type_num, fig_num, savefig_flag, plot_output_flag)
    println("Finished packaging filter solution.")

    return h, 𝓧_new, fig_num
end







function frontendfilterdesignfast(filter_type_num, zfunc, L, fig_num,
                                savefig_flag, plot_output_flag,
                                config::MiniMaxConfigType)

    𝓧_tmp = SharedArray{Float64}(L+2)
    @sync @distributed for k = 0:L+1
        𝓧_tmp[k+1] = convert(Float64, Chebyshev2ndnode(k,L+1))
    end
    𝓧 = convert(Vector{BigFloat}, 𝓧_tmp)

    ##
    wfunc = getweightfunc(filter_type_num)
    dfunc = xx->zfunc(acos(xx))

    𝓧_selection_string = "Chebyshev"
    candidate_multiple = 1

    println("Begin mini-max algorithm.")
    @time 𝓧_new, h_history,
            iters_ran, status_msg, fig_num = designfilter(  config,
                                                            L,
                                                            candidate_multiple,
                                                            wfunc,
                                                            dfunc,
                                                            𝓧_selection_string,
                                                            fig_num)

    #

    println("Completed mini-max algorithm. Status: ", status_msg)

    println("Start to package filter solution.")
    @time h = packagefiltersolutionfast(𝓧_new, h_history, wfunc, zfunc, dfunc,
                                        filter_type_num)
    println("Finished packaging filter solution.")

    return h, 𝓧_new, fig_num
end

function visualizeoptimzationsolution(xq, 𝓧, f, fig_num)

    𝑤_new = get𝑤scaled(𝓧)
    h_new = geth(𝑤_new, constantfunc, f, 𝓧)

    f_itp = xx->getp(xx, 𝑤_new, h_new, constantfunc, f, 𝓧)
    f_itp_xq = collect( f_itp(xq[i]) for i = 1:length(xq) )
    f_itp_𝓧 = collect( f_itp(𝓧[i]) for i = 1:length(𝓧) )

    𝑒 = xx->(f(xx)-f_itp(xx))
    𝑒_xq = collect( 𝑒(xq[i]) for i = 1:length(xq) )
    𝑒_𝓧 = collect( 𝑒(𝓧[i]) for i = 1:length(𝓧) )

    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(xq, f.(xq), "--", label = "target")
    # PyPlot.plot(xq, f_itp_xq, label = "solution")
    # PyPlot.plot(𝓧, f_itp_𝓧, ".")
    # PyPlot.title("f vs. solution")
    # PyPlot.legend()
    #
    # PyPlot.figure(fig_num)
    # fig_num += 1
    # PyPlot.plot(xq, 𝑒_xq)
    # PyPlot.plot(𝓧, 𝑒_𝓧, ".")
    # PyPlot.title("error")
    # PyPlot.legend()

    return f_itp_xq, f_itp_𝓧, 𝑒_xq, 𝑒_𝓧, fig_num
end

function frontendoptimization(targetfunc, L,
            opt_parameters_file_name::String = "0")

    𝓧 = collect( Chebyshev2ndnode(k,L+1) for k = 0:L+1 )

    config = getdefaultoptparameters(L)
    if opt_parameters_file_name != "0"
        config = getoptparameters(L,opt_parameters_file_name)
    end

    𝓧_new, h_history, iters_ran, status_msg, fig_num_unused = fitminimax(  constantfunc,
                                                                targetfunc,
                                                                𝓧,
                                                                1,
                                                                config)
    return 𝓧_new
end
