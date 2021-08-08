
# M is the number of positions to keep
# 𝓟 is a collection of candidate positions.
function removepointssimple!(𝓟::Vector{T}, 𝑒::Function, M::Int, h::T, tol::T) where T
    @assert length(𝓟) > M

    abs_𝑒_𝓟 = abs.(𝑒.(𝓟))

    if length(𝓟) > M
        # remove the positions associated with the smallest errors.
        ind_sort = sortperm(abs_𝑒_𝓟,rev=true)

        out = 𝓟[ind_sort[1:M]]

        resize!(𝓟,length(out))
        𝓟[:] = out
    end

    @assert length(𝓟) == M

    return nothing
end

function lg7(a::Vector{T}, N::Int) where T
    𝓟 = Vector{T}(undef,0)
    removepointsAlg7!(𝓟,a,N)

    return 𝓟
end

# Algorithm 7 from Filip's thesis.
# N is the number of points to keep.
function removepointsAlg7!(𝓟::Vector{T}, 𝑒::Function, N) where T

    ### preparation steps before removing positions.
    sort!(𝓟)

    # for each interval, find the position that has the maximum absolute error.
    𝓐 = forcealternating(𝓟, 𝑒)

    # This should be true since 𝓟 should include the previous reference set, which
    #   has N alternating error signs.
    @assert length(𝓐) >= N

    if length(𝓐) == N
        return 𝓐
    end

    ##
    N_remove = length(𝓐) - N

    # remove the first or last position if N_remove cannot be fulfilled by
    #   removing only pairs of positions.
    if isodd(N_remove)
        if abs( 𝑒(𝓐[1]) ) > abs( 𝑒(𝓐[end]) )
            deleteat!(𝓐, length(𝓐))
        else
            deleteat!(𝓐, 1)
        end
        N_remove -= 1
    end
    @assert iseven(N_remove)

    # iteratively prune pairs of values.
    while N_remove > 0

        tmp = collect( max( abs(𝑒(𝓐[i])), abs(𝑒(𝓐[i+1])) ) for i = 1:length(𝓐)-1 )
        min_e, min_ind = findmin(tmp)

        # check if the wrap-around pair has the smallest absolute error.
        if max( abs(𝑒(𝓐[1])), abs(𝑒(𝓐[end])) ) < min_e
            # remove the abs error from wrap-around positions.
            deleteat!(𝓐, (1, length(𝓐)))
        else
            # remove smallest abs error.
            deleteat!(𝓐, (min_ind, min_ind+1))
        end

        N_remove -= 2
    end

    @assert length(𝓐) == N
    return 𝓐
end



function addpointsinwithinterval(a, N)
    𝓟 = Vector{T}(undef,0)
    addpointsinwithinterval!(𝓟, a, N)

    return 𝓟
end

# N is the number of points to add to a.
# assumes a is sorted in ascending order.
function addpointsinwithinterval!(𝓟::Vector{T}, a::Vector{T}, N)::Nothing where T
    M = length(a)

    ### add at least s number of points between every pair of extrema_positions
    s = div(N,M-1)
    r = N - s*(M-1)

    # reset container.
    resize!(𝓟, 0)

    # add points in the first interval.
    push!(𝓟, a[1])

    st = a[1]

    fin = a[2]
    tmp = collect(LinRange(st,fin,s+r+2)) # +2 since we want the end points, which are the extrema.

    push!(𝓟, tmp[2:end]...)

    st = fin

    # add points for the rest of the intervals.
    for i = 2:M-1
        # add points.
        fin = a[i+1]
        tmp = collect(LinRange(st,fin,s+2)) # +2 since we want the end points, which are the extrema.
        push!(𝓟, tmp[2:end]...)

        # update.
        st = fin
    end
    @assert length(𝓟) == N+M

    return nothing
end


function isbettercandidate(a::T, b::T, f::Function, η::T = one(T)) where T
    score_a = abs(Calculus.derivative(f,a)) + one(T)/abs(f(a))
    score_b = abs(Calculus.derivative(f,b)) + one(T)/abs(f(b))

    return score_a < score_b
end
