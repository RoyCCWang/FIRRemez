
function getdefaultoptparameters(L)
    tol_converged = 1e-3
    tol_𝓧_spacing = 1e-6
    tol_derivative_zero = 0.2
    tol_no_update = 1e-6
    #tol_h = 1e-17
    N_samples_interval = 10
    N_test_positions = 200
    max_iters = 200
    L_interpolant = 4

    tol_converged = convert(BigFloat, tol_converged )
    tol_𝓧_spacing = convert(BigFloat, tol_𝓧_spacing/(L+2))
    tol_derivative_zero = convert(BigFloat, tol_derivative_zero )
    tol_no_update = convert(BigFloat, tol_no_update )
    #tol_h = convert(BigFloat, tol_h )
    N_samples_interval = convert(Int,N_samples_interval)
    N_test_positions = convert(Int,N_test_positions)
    max_iters = convert(Int,max_iters)
    L_interpolant = convert(Int,L_interpolant)
    verbose_flag = true
    plot_flag = false

    config = MiniMaxConfigType( tol_converged,
                                tol_𝓧_spacing,
                                tol_derivative_zero,
                                tol_no_update,
                                #tol_h,
                                N_samples_interval,
                                N_test_positions,
                                max_iters,
                                L_interpolant,
                                verbose_flag,
                                plot_flag)

    return config
end

function getoptparameters(L,file_name::String)

    tol_converged,
    tol_𝓧_spacing,
    tol_derivative_zero,
    tol_no_update,
    #tol_h,
    N_samples_interval,
    N_test_positions,
    max_iters,
    L_interpolant,
    verbose_flag_input,
    plot_flag_input = DelimitedFiles.readdlm(file_name, ',', Float64)

    tol_converged = convert(BigFloat, tol_converged )
    tol_𝓧_spacing = convert(BigFloat, tol_𝓧_spacing/(L+2))
    tol_derivative_zero = convert(BigFloat, tol_derivative_zero )
    tol_no_update = convert(BigFloat, tol_no_update )
    #tol_h = convert(BigFloat, tol_h )
    N_samples_interval = convert(Int,N_samples_interval)
    N_test_positions = convert(Int,N_test_positions)
    max_iters = convert(Int,max_iters)
    L_interpolant = convert(Int,L_interpolant)
    verbose_flag = verbose_flag_input==1 ? true : false
    plot_flag = plot_flag_input==1 ? true : false

    config = MiniMaxConfigType( tol_converged,
                                tol_𝓧_spacing,
                                tol_derivative_zero,
                                tol_no_update,
                                #tol_h,
                                N_samples_interval,
                                N_test_positions,
                                max_iters,
                                L_interpolant,
                                verbose_flag,
                                plot_flag)
    return config
end
