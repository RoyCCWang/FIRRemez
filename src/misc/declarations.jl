

mutable struct MiniMaxConfigType{T}
      tol_converged::T # used to compare the score; this is the convergence criterion in step 4 of algorithm 2 from Filip's thesis.
      tol_𝓧_spacing::T # minimum distance between reference positions.
      tol_derivative_zero::T
      tol_no_update::T
      #tol_h::T # new reference points are allowed to be this much
               #  smaller than the absolute level error.

      N_samples_interval::Int # number of samples to take for each error function's interval.
      N_test_positions::Int # number of samples used to evaluate the convergence score.
      max_iters::Int
      L_interval::Int

      verbose_flag::Bool
      plot_flag::Bool # if true, plots error function and new reference positions every iteration.
 end
