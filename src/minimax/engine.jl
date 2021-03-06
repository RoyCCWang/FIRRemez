# Nq is the number of positions used to test convergence.
function fitminimax( wfunc::Function,
                    dfunc::Function,
                    ๐ง0::Vector{T},
                    fig_num,
                    config) where T

    ##### parse.
    tol_converged = config.tol_converged
    tol_๐ง_spacing = config.tol_๐ง_spacing
    tol_derivative_zero = config.tol_derivative_zero
    tol_no_update = config.tol_no_update
    N_samples_interval = config.N_samples_interval
    N_test_positions = config.N_test_positions
    max_iters = config.max_iters
    L_interval = config.L_interval

    ##### set up.
    ๐ง = copy(๐ง0)
    L = length(๐ง) - 2

    # diagnostics.
    level_error_history = Vector{T}(undef,max_iters)
    xq = convert(Vector{BigFloat},collect(LinRange(-one(T),one(T),N_test_positions)))

    # pre-allocate.
    ๐ค = get๐คscaled(๐ง)
    h = geth(๐ค, wfunc, dfunc, ๐ง)
    ๐ = collect( dfunc(๐ง[i]) - (-1)^(i-1)*h/wfunc(๐ง[i]) for i = 1:length(๐ง))
    int_storage_Lโ = SharedArray{Float64}(length(xq))

    f = xx->wfunc(xx)*dfunc(xx)
    f_itp = xx->getp(xx, ๐, ๐ค, h, wfunc, dfunc, ๐ง)
    ๐ = xx->(f(xx)-f_itp(xx))

    sample_position_template = collect( Chebyshev2ndnode(k,config.N_samples_interval) for k = 0:config.N_samples_interval )

    ##### run algorithm 2.
    iter = 1

    # repeat until convergence.
    while iter <= max_iters

        ##### update interpolant.
        get๐คscaled!(๐ค, ๐ง)
        h = geth(๐ค, wfunc, dfunc, ๐ง)

        for i = 1:length(๐ง)
            ๐[i] = dfunc(๐ง[i]) - (-1)^(i-1)*h/wfunc(๐ง[i])
        end

        ##### diagnostics
        level_error_history[iter] = h

        ### brute-force uniform sampling method to find the position in [-1,1]
        #       with the largest absolute error.
        Lโ_value, Lโ_ind, ๐_at_Lโ_position = discretizeLโnormparallel!(int_storage_Lโ, f,f_itp,xq)

        sign_๐_at_Lโ_position = sign(๐_at_Lโ_position) > 0

        Lโ_position = xq[Lโ_ind]
        score = (Lโ_value - abs(h))/Lโ_value

        ### outputs.
        if config.verbose_flag
            println("iter: ", iter)
            println("h: ", h)
            println("score: ", score)
        end

        if score <= tol_converged
            resize!(level_error_history, iter)

            status_msg = "Converged"
            return ๐ง, level_error_history, iter, status_msg, fig_num
        end

        ### find next set of reference positions.
        fig_num, proceed_flag, status_msg = innerstep(๐ง, h,
                                                        f, f_itp, ๐,
                                                        wfunc, dfunc, L_interval,
                                                        L,
                                                        xq,
                                                        Lโ_value,
                                                        Lโ_position,
                                                        sign_๐_at_Lโ_position,
                                                        sample_position_template,
                                                        iter,
                                                        fig_num,
                                                        config)

        if config.verbose_flag
            println()
        end


        if !proceed_flag
            resize!(level_error_history, iter)
            println(status_msg)
            return ๐ง, level_error_history, iter, status_msg, fig_num
        end

        iter += 1
    end

    status_msg = "Max. iterations reached."
    println(status_msg)
    return ๐ง, level_error_history, iter, status_msg, fig_num
end

# interval extrema search.
function innerstep(๐ง::Vector{T},
                    h::T,
                    f::Function,
                    f_itp::Function,
                    ๐::Function,
                    wfunc::Function,
                    dfunc::Function,
                    L_interval::Int,
                    L::Int,
                    xq,
                    Lโ_value::Float64,
                    Lโ_position::T,
                    sign_๐_at_Lโ_position::Bool,
                    sample_position_template,
                    iter,
                    fig_num,
                    config) where T

    tol_๐ง_spacing = config.tol_๐ง_spacing
    tol_derivative_zero = config.tol_derivative_zero
    tol_no_update = config.tol_no_update

    if config.plot_flag
        # PyPlot.figure(fig_num)
        # fig_num += 1
    end

    if config.plot_flag
        min_dist, j = findmin(abs.( ๐ง .- Lโ_position ))
        x_j = ๐ง[j]
    end

    ##### find error extrema positions.
    ๐, ๐_๐, d_set, X_set = getcandidatesfromextremas( ๐,
                    L_interval, ๐ง,
                    sample_position_template,
                    tol_๐ง_spacing, tol_derivative_zero, h)

    if config.plot_flag
        ๐_before_processing = copy(๐)
    end

    ๐ง_new, status_flag = processreferenceset(  ๐ง,
                                                h,
                                                ๐,
                                                ๐,
                                                ๐_๐,
                                                Lโ_value,
                                                Lโ_position,
                                                sign_๐_at_Lโ_position,
                                                tol_๐ง_spacing,
                                                #config.tol_h,
                                                iter,
                                                config.verbose_flag)


    if config.plot_flag
        # xq = convert(Vector{BigFloat},collect(LinRange(-1,1,20*length(๐ง))))
        # boundary_positions = collect( X_set[i][1] for i = 1:length(X_set))
        # push!(boundary_positions, X_set[end][end])
        #
        # d, ฮฝ = getChebyshevinterpolator(๐, L)
        # ๐_Chebyitp_interval = xx->evalChebyshevpolynomialinterval(xx, d_set, X_set)
        #
        # PyPlot.plot(๐ง, ๐.(๐ง), "+", label = "e(X)")
        #
        # PyPlot.plot(xq, ๐.(xq), label = "e()")
        #
        # PyPlot.plot(๐ง_new, ๐.(๐ง_new), "x", label = "X_new")
        #
        # println("x_j is ", x_j)
        # println("Lโ_position is ", Lโ_position)
        #
        #
        # title_string = Printf.@sprintf("๐_Chebyitp_interval: boundary and extrema, iter = %d", iter)
        # PyPlot.title(title_string)
        # PyPlot.legend()
    end

    ### update.

    # check if we're stuck.
    if sum(abs.(๐ง-๐ง_new)) < tol_no_update
        println("Reference set no longer changing.")
        return fig_num, false, "Reference set no longer changing."
    end

    ๐ง[:] = ๐ง_new

    return fig_num, true, ""
end
