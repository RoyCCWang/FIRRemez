
function processinterval(   𝑖::Int,
                            𝓟::Vector{T},
                            sample_position_template::Vector{T},
                            𝑒::Function,
                            L::Int,
                            verbose_flag) where T
    @assert 1 < 𝑖 <= length(𝓟)

    a = 𝓟[𝑖-1]
    b = 𝓟[𝑖]

    # build approximation of f over current interval.
    X = samplepoints( sample_position_template, a, b )

    𝑒_itp_X = xx->Barycentric2nditp(xx, X, 𝑒.(X))

    # fit Chebyshevinterpolator over this approximation.
    d, ν = getChebyshevinterpolator(𝑒_itp_X, L, a, b)

    𝑒_extrema_positions_X = scalefromChebyinterval.(findChebyitpextremawithchk(d), a, b)

    # d and X can be used for plotting or debugging purposes.
    return 𝑒_extrema_positions_X, d, X
end


# sample_position_template is between [-1,1]. For each interval, rescale it to the interval at hand.
function getextremafromroots(  𝑒::Function,
                                L::Int,
                                𝓟_in::Vector{T},
                                sample_position_template::Vector{T},
                                tol_𝓧_spacing, # was defaulted to convert(T,1e-6/length(𝓟_in)
                                tol_derivative_zero,
                                verbose_flag = true) where T
    𝓟 = Vector{T}(undef,length(𝓟_in)+2)
    𝓟[1] = -one(T)
    𝓟[2:end-1] = 𝓟_in
    𝓟[end] = one(T)

    makenotclose!(𝓟, convert(T,1e-14/length(𝓟_in))) # just to be sure positions are spread out enough.

    # prepare outputs.
    out = Vector{T}(undef,0)
    d_set = Vector{Vector{T}}(undef,length(𝓟)-1)
    X_set = Vector{Vector{T}}(undef,length(𝓟)-1)

    # work on intervals.
    for i = 2:length(𝓟)

        # build approximation of f over current interval.
        X = samplepoints( sample_position_template,𝓟[i-1],𝓟[i] )

        𝑒_itp_X = xx->Barycentric2nditp(xx, X, 𝑒.(X))

        # fit Chebyshevinterpolator over this approximation.
        d, ν = getChebyshevinterpolator(𝑒_itp_X, L, 𝓟[i-1], 𝓟[i])

        if verbose_flag
            Printf.@printf("working on %d-th interval\n",i-1)
        end
        𝑒_extrema_positions_X = scalefromChebyinterval.(findChebyitpextremawithchk(d), 𝓟[i-1], 𝓟[i])

        push!(out, 𝑒_extrema_positions_X...)

        # for plotting or debugging.
        d_set[i-1] = d
        X_set[i-1] = X
    end

    # remove extrema search results that are too close to zero.
    #   These are probably incorrect results, since the local extremas
    #   should be greater or equal to the level error, which is non-zero
    #   unless the target function is in the space of Chebyshev approximations.
    while sign(𝑒(out[1])) == zero(T) && !isempty(out)
        popfirst!(out)
    end

    return out, d_set, X_set
end

function getextremafromrootsparallel(𝑒::Function,
                                L::Int,
                                𝓧::Vector{T},
                                sample_position_template::Vector{T},
                                tol_𝓧_spacing, # was defaulted to convert(T,1e-6/length(𝓟_in)
                                tol_derivative_zero,
                                verbose_flag = true) where T

    ## construct the boundary point sets that will dictate the piece-wise proxy intervals.
    𝓑 = copy(𝓧)
    if !isnumericallyclose(-one(T), 𝓧[1])
        insert!(𝓑,1,-one(T))
    end
    if !isnumericallyclose(one(T), 𝓧[end])
        push!(𝓑,one(T))
    end

    # work on intervals.
    workerfunc = xx->processinterval(xx, 𝓑, sample_position_template, 𝑒, L, verbose_flag)
    sol = pmap(workerfunc, collect(2:length(𝓑)) )

    # package up outputs.
    A = sol
    k = 1
    len = sum( length(A[i][k]) for i = 1:length(A))
    out = Vector{T}(undef,len)

    st = 1
    for i = 1:length(A)
        fin = st + length(A[i][k]) -1
        out[st:fin] = A[i][k]

        st = fin + 1
    end

    d_set = collect( A[i][2] for i = 1:length(A) )
    X_set = collect( A[i][3] for i = 1:length(A) )

    return out, d_set, X_set
end
