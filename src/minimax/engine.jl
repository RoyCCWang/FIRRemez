# Nq is the number of positions used to test convergence.
function fitminimax( wfunc::Function,
                    dfunc::Function,
                    𝓧0::Vector{T},
                    fig_num,
                    config) where T

    ##### parse.
    tol_converged = config.tol_converged
    tol_𝓧_spacing = config.tol_𝓧_spacing
    tol_derivative_zero = config.tol_derivative_zero
    tol_no_update = config.tol_no_update
    N_samples_interval = config.N_samples_interval
    N_test_positions = config.N_test_positions
    max_iters = config.max_iters
    L_interval = config.L_interval

    ##### set up.
    𝓧 = copy(𝓧0)
    L = length(𝓧) - 2

    # diagnostics.
    level_error_history = Vector{T}(undef,max_iters)
    xq = convert(Vector{BigFloat},collect(LinRange(-one(T),one(T),N_test_positions)))

    # pre-allocate.
    𝑤 = get𝑤scaled(𝓧)
    h = geth(𝑤, wfunc, dfunc, 𝓧)
    𝑓 = collect( dfunc(𝓧[i]) - (-1)^(i-1)*h/wfunc(𝓧[i]) for i = 1:length(𝓧))
    int_storage_L∞ = SharedArray{Float64}(length(xq))

    f = xx->wfunc(xx)*dfunc(xx)
    f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝓧)
    𝑒 = xx->(f(xx)-f_itp(xx))

    sample_position_template = collect( Chebyshev2ndnode(k,config.N_samples_interval) for k = 0:config.N_samples_interval )

    ##### run algorithm 2.
    iter = 1

    # repeat until convergence.
    while iter <= max_iters

        ##### update interpolant.
        get𝑤scaled!(𝑤, 𝓧)
        h = geth(𝑤, wfunc, dfunc, 𝓧)

        for i = 1:length(𝓧)
            𝑓[i] = dfunc(𝓧[i]) - (-1)^(i-1)*h/wfunc(𝓧[i])
        end

        ##### diagnostics
        level_error_history[iter] = h

        ### brute-force uniform sampling method to find the position in [-1,1]
        #       with the largest absolute error.
        L∞_value, L∞_ind, 𝑒_at_L∞_position = discretizeL∞normparallel!(int_storage_L∞, f,f_itp,xq)

        sign_𝑒_at_L∞_position = sign(𝑒_at_L∞_position) > 0

        L∞_position = xq[L∞_ind]
        score = (L∞_value - abs(h))/L∞_value

        ### outputs.
        if config.verbose_flag
            println("iter: ", iter)
            println("h: ", h)
            println("score: ", score)
        end

        if score <= tol_converged
            resize!(level_error_history, iter)

            status_msg = "Converged"
            return 𝓧, level_error_history, iter, status_msg, fig_num
        end

        ### find next set of reference positions.
        fig_num, proceed_flag, status_msg = innerstep(𝓧, h,
                                                        f, f_itp, 𝑒,
                                                        wfunc, dfunc, L_interval,
                                                        L,
                                                        xq,
                                                        L∞_value,
                                                        L∞_position,
                                                        sign_𝑒_at_L∞_position,
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
            return 𝓧, level_error_history, iter, status_msg, fig_num
        end

        iter += 1
    end

    status_msg = "Max. iterations reached."
    println(status_msg)
    return 𝓧, level_error_history, iter, status_msg, fig_num
end

# interval extrema search.
function innerstep(𝓧::Vector{T},
                    h::T,
                    f::Function,
                    f_itp::Function,
                    𝑒::Function,
                    wfunc::Function,
                    dfunc::Function,
                    L_interval::Int,
                    L::Int,
                    xq,
                    L∞_value::Float64,
                    L∞_position::T,
                    sign_𝑒_at_L∞_position::Bool,
                    sample_position_template,
                    iter,
                    fig_num,
                    config) where T

    tol_𝓧_spacing = config.tol_𝓧_spacing
    tol_derivative_zero = config.tol_derivative_zero
    tol_no_update = config.tol_no_update

    if config.plot_flag
        # PyPlot.figure(fig_num)
        # fig_num += 1
    end

    if config.plot_flag
        min_dist, j = findmin(abs.( 𝓧 .- L∞_position ))
        x_j = 𝓧[j]
    end

    ##### find error extrema positions.
    𝓟, 𝑒_𝓟, d_set, X_set = getcandidatesfromextremas( 𝑒,
                    L_interval, 𝓧,
                    sample_position_template,
                    tol_𝓧_spacing, tol_derivative_zero, h)

    if config.plot_flag
        𝓟_before_processing = copy(𝓟)
    end

    𝓧_new, status_flag = processreferenceset(  𝓧,
                                                h,
                                                𝑒,
                                                𝓟,
                                                𝑒_𝓟,
                                                L∞_value,
                                                L∞_position,
                                                sign_𝑒_at_L∞_position,
                                                tol_𝓧_spacing,
                                                #config.tol_h,
                                                iter,
                                                config.verbose_flag)


    if config.plot_flag
        # xq = convert(Vector{BigFloat},collect(LinRange(-1,1,20*length(𝓧))))
        # boundary_positions = collect( X_set[i][1] for i = 1:length(X_set))
        # push!(boundary_positions, X_set[end][end])
        #
        # d, ν = getChebyshevinterpolator(𝑒, L)
        # 𝑒_Chebyitp_interval = xx->evalChebyshevpolynomialinterval(xx, d_set, X_set)
        #
        # PyPlot.plot(𝓧, 𝑒.(𝓧), "+", label = "e(X)")
        #
        # PyPlot.plot(xq, 𝑒.(xq), label = "e()")
        #
        # PyPlot.plot(𝓧_new, 𝑒.(𝓧_new), "x", label = "X_new")
        #
        # println("x_j is ", x_j)
        # println("L∞_position is ", L∞_position)
        #
        #
        # title_string = Printf.@sprintf("𝑒_Chebyitp_interval: boundary and extrema, iter = %d", iter)
        # PyPlot.title(title_string)
        # PyPlot.legend()
    end

    ### update.

    # check if we're stuck.
    if sum(abs.(𝓧-𝓧_new)) < tol_no_update
        println("Reference set no longer changing.")
        return fig_num, false, "Reference set no longer changing."
    end

    𝓧[:] = 𝓧_new

    return fig_num, true, ""
end
