## this script contain functions related to the reference set
#   exchange strategy. It is a multi-point exchange strategy, that
#   encourages alternation in error signs.
# given:
#   - an array, a, of pairs of real numbers (x,y)
#   - a positive number h
#   - an integer N such that N+2 <= length of a.
# find the N+2 entries in that encourages:
#   - the y values alternate in sign, when ordered by x.
#   - the y values are among the largest.



function getcandidatesfromextremas( 𝑒::Function,
                                    L::Int,
                                    𝓧::Vector{T},
                                    sample_position_template::Vector{T},
                                    tol_𝓧_spacing, # was defaulted to convert(T,1e-6/length(𝓟_in)
                                    tol_derivative_zero,
                                    h) where T

    ### get roots.
    extremas, d_set, X_set = getextremafromrootsparallel( 𝑒, L, 𝓧,
                                            sample_position_template,
                                            tol_𝓧_spacing, tol_derivative_zero)
    #
    # if not enough extrema, there might be a numerical issue,
    #   or that the interpolant degree is not high enough, or that the extrema is
    #   at one of the boundary points in 𝓟 (i.e. in the reference set 𝓧).
    # add the reference set.

    # enforce the extremas positions must have an absolute error that is larger than h.
    N_samples = 5
    𝓟, 𝑒_𝓟 = refinecandidates!( extremas,
                                𝑒,
                                h,
                                N_samples)

    return 𝓟, 𝑒_𝓟, d_set, X_set
end


# Adds uniformly spaced positions between the closest candidate to -1. Similarly for 1.
# Removes candidates that have an error that is less than the level error h.
# These are incorrect results, since the local extremas
#   should be greater or equal to the level error, which is non-zero
#   unless the target function is in the space of Chebyshev approximations.
#
# ϵ is used for comparing two floating point numbers.
# Keep only the entries in a that are larger than the level error h.
# add new entries that are between the kept entires.
# modifies candidate_array.
# f here is the error function.
# The output is sorted in terms of the level error, not the candidate position.
function refinecandidates!( candidate_array::Vector{T},
                            f::Function,
                            h,
                            N_samples::Int,
                            ϵ::T = eps(T)*2) where T
    #
    unique!(candidate_array)

    # add endpoints and equally spaced points between each end point to the
    #   closest candidate, to the candidate list.
    if !isnumericallyclose(-one(T), candidate_array[1])
        tmp = collect( LinRange(-one(T), candidate_array[1], N_samples) )
        push!(candidate_array, tmp[1:end-1]...)
    end

    if !isnumericallyclose(one(T), candidate_array[end])

        tmp = collect( LinRange(candidate_array[end], one(T), N_samples) )
        push!(candidate_array, tmp[2:end]...)
    end

    # Evaluate the error for each candidate.
    error_array = SharedArray{Float64}(length(candidate_array))
    @sync @distributed for i = 1:length(candidate_array)
        error_array[i] = convert(Float64, f(candidate_array[i]))
    end


    ## discard the candidates that have an absolute error less than the level error.
    sort_ind = sortperm(abs.(error_array), rev=true) # descending order.
    candidate_𝓧_array = candidate_array[sort_ind]
    candidate_𝑒_array = error_array[sort_ind]

    k = 1
    while k <= length(candidate_𝑒_array) && abs(candidate_𝑒_array[k]) > h
        k += 1
    end
    k -= 1

    if k > 0
        # among the candidates, at least one is as large as the level error.
        # Keep the ones that are larger or equal to the level error.
        candidate_𝑒_array = convert(Vector{T}, candidate_𝑒_array[1:k])
        candidate_𝓧_array = candidate_𝓧_array[1:k]

        return candidate_𝓧_array, candidate_𝑒_array
    end

    # case: no candidates that have absolute error that is > the level error.
    return Vector{T}(undef,0), Vector{T}(undef,0)
end




# Use the candidates and the current reference set to update the reference set.
function processreferenceset!(  𝓧::Vector{T},
                                h,
                                𝑒::Function,
                                𝓟::Vector{T},
                                𝑒_𝓟::Vector{T},
                                L∞_value::Float64,
                                L∞_position::T,
                                sign_𝑒_at_L∞_position::Bool,
                                tol_𝓧_spacing,
                                #tol_h,
                                iter,
                                verbose_flag::Bool = true,
                                L∞_allowed_porpotion::Float64 = 0.95)::Bool where T
#
    @assert length(𝓟) == length(𝑒_𝓟)

    ## sanity-check: 𝑒(𝓧) must be alternating in sign.
    𝑒_𝓧_array = SharedArray{Float64}(length(𝓧))
    @sync @distributed for i = 1:length(𝓧)
        𝑒_𝓧_array[i] = convert(Float64, 𝑒(𝓧[i]) )
    end

    for i = 2:length(𝑒_𝓧_array)
        @assert !xor( sign(𝑒_𝓧_array[i-1]) > 1.0, sign(𝑒_𝓧_array[i]) > 1.0 )
    end

    ### Decide if doing multi-point exchange or one point exchange.
    # prune one entry so that the number of elements to be removed becomes even.

    if findfirst(abs.(𝑒_𝓟).>=L∞_allowed_porpotion*L∞_value) == nothing
        println("Using single exchange")
        # single exchange: replace an entry of the reference set with the
        #   L∞_position.
        # The entry should be chosen such that it is closest to L∞_position, but
        #   has the same error sign as the L∞_value

        min_dist_unused, j = findmin(abs.( 𝓧 .- L∞_position ))

        # find second minimum. Use the fact that 𝓧 is sorted in ascending order.
        j2 = j+1
        if abs( 𝓧[j-1] - L∞_position ) < abs( 𝓧[j+1] - L∞_position )
            j2 = j-1
        end

        sign_j = sign(𝑒_𝓧_array[j]) > 0
        sign_j2 = sign(𝑒_𝓧_array[j2]) > 0

        if sign_𝑒_at_L∞_position == sign_j
            𝓧[j] = L∞_position

            println("sign match as is")

        else

            # sanity check: the errors of 𝓧[j] and 𝓧[j2] should be alternating,
            #   since they are adjacent positions.
            @assert sign_𝑒_at_L∞_position == sign_j2

            𝓧[j2] = L∞_position

            println("sign match, use closest neighbour with same sign.")
        end


        sort!(𝓧)
        return true
    end


    ### Do multi-point exchange.
    println("Using multiple exchange")
    push!(𝓟, 𝓧...)

    ## discard entries that are too close in proximity.
    makenotclose!(𝓟, convert(T,1e-14))
    #unique!(𝓟)
    @assert length(𝓟) >= length(𝓧)

    if length(𝓟) > length(𝓧)
        𝓐 = removepointsAlg7!(𝓟, 𝑒, length(𝓧))
    end
    @assert length(𝓐) == length(𝓧)

    # check exchange condition.
    #@assert all(abs.(𝑒.(𝓐)).>=h)

    ### sanity-check: 𝑒(𝓧) must be alternating in sign.
    𝑒_𝓐_array = SharedArray{Float64}(length(𝓐))
    @sync @distributed for i = 1:length(𝓐)
        𝑒_𝓐_array[i] = convert(Float64, 𝑒(𝓐[i]) )
    end

    for i = 2:length(𝑒_𝓐_array)
        if xor( sign(𝑒_𝓐_array[i-1]) > 1.0, sign(𝑒_𝓐_array[i]) > 1.0 )

            # The errors of 𝓐 should be alternating.
            if verbose_flag
                Printf.@printf("Error (iter %d): 𝑒(𝓐) is not alternating! Exit.", iter )
            end

            return false
        end
    end

    ### sanity-check:
    if findfirst(abs.(𝑒_𝓐_array).>=L∞_allowed_porpotion*L∞_value) == nothing
        if verbose_flag
            Printf.@printf("Error (iter %d): none of the entries in 𝑒(𝓧_new) is near max(𝑒([-1,1])). Exit.", iter )
        end

        return false
    end

    ### safe to update reference set.
    𝓧[:] = 𝓐
    sort!(𝓧)

    return true
end

function processreferenceset(   𝓧::Vector{T},
                                h,
                                𝑒,
                                𝓟::Vector{T},
                                𝑒_𝓟::Vector{T},
                                L∞_value::Float64,
                                L∞_position::T,
                                sign_𝑒_at_L∞_position::Bool,
                                tol_𝓧_spacing,
                                #tol_h,
                                iter,
                                verbose_flag::Bool = true) where T
    𝓧_new = copy(𝓧)
    status_flag = processreferenceset!(𝓧_new, h, 𝑒, 𝓟, 𝑒_𝓟,
                                        L∞_value,
                                        L∞_position,
                                        sign_𝑒_at_L∞_position,
                                        tol_𝓧_spacing,
                                        #tol_h,
                                        iter, verbose_flag)

    return 𝓧_new, status_flag
end
