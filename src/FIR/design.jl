
##### filter design routines.


function designfilter(  config::MiniMaxConfigType,
                        L::Int,
                        candidate_multiple::Int,
                        wfunc::Function,
                        dfunc::Function,
                        𝓧_selection_string::String,
                        fig_num::Int)


#
    ### prepare initial reference positions.
    𝓧 = Vector{BigFloat}(undef,0)
    if 𝓧_selection_string == "Chebyshev"
        𝓧 = collect( Chebyshev2ndnode(k,L+1) for k = 0:L+1 )

    elseif 𝓧_selection_string == "maxvol via QR"
        𝓧, 𝑖_unused = getinitialpositions(L, candidate_multiple)

    else
        println("Unknown 𝓧 selection string. Using Chebyshev.")
        𝓧 = collect( Chebyshev2ndnode(k,L+1) for k = 0:L+1 )
    end
    @assert length(unique(𝓧)) == length(𝓧) # sanity check

    ### run mini-max.
    𝓧_new, h_history, iters_ran, status_msg, fig_num = fitminimax(  wfunc,
                                                                dfunc,
                                                                𝓧,
                                                                fig_num,
                                                                config)


    ### get filter coefficients.

    # pose the solution as a Lagrange interpolated function.




    return 𝓧_new, h_history, iters_ran, status_msg, fig_num
end
