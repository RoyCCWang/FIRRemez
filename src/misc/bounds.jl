### routines for checking proximity or intervals for floating-point numbers.

function isnumericallyclose(x::T, y::T, tol = eps(T)*2) where T
    if abs(x-y) < tol
        return true
    end

    return false
end

function isnumericallyin(x::T, a::T, b::T, tol = eps(x)*2)::Bool where T
    if x > a-tol && x < b+tol
        return true
    end

    return false
end

# remove points that are too close to the entry before it.
# also sorts the collect A.
function makenotclose!(A, tol)::Nothing
    remove_list = Vector{Int}(undef,0)
    sort!(A)

    for i = 2:length(A)
        if isnumericallyclose(A[i-1], A[i], tol)
            push!(remove_list, i)
        end
    end

    deleteat!(A,remove_list)

    return nothing
end
