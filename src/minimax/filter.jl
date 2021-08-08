function constantfunc(x)::Int
    return 1
end

function type2weightfunc(x::T, tol::T = eps(T))::T where T
    return cos(acos(x)/2)
end

function type3weightfunc(x::T, tol::T = eps(T))::T where T
    return sin(acos(x)) + 0.0001
end

function type4weightfunc(x::T, tol::T = eps(T))::T where T
    return sin(acos(x)/2) + 0.0001
end

function getweightfunc(type_number::Int)::Function
    @assert 1 <= type_number <= 4
    if type_number == 1
        return constantfunc
    elseif type_number == 2
        return type2weightfunc
    elseif type_number == 3
        return type3weightfunc
    end

    return type4weightfunc
end


# uses the Barycentric interpolation, Equation 3.17 of Filip's thesis.
function levellederror(wfunc::Function, w::Vector{T}, 𝓧, 𝓨::Vector{T})::T where T

    tmp_numerator = zero(T)
    tmp_denominator = zero(T)
    for k = 1:length(𝓧)

        tmp_numerator += 𝓨[k]*w[k]/wfunc(𝓧[k])

        tmp = w[k]/wfunc(𝓧[k])
        if iseven(k)
            tmp_denominator += tmp
        else
            tmp_denominator -= tmp
        end

    end

    @assert tmp_denominator != zero(T)
    return tmp_numerator/tmp_denominator
end

function levellederror(wfunc::Function, 𝓧, 𝓨::Vector{T})::T where T
    𝑤 = get𝑤scaled(𝓧)
    return levellederror(wfunc, w, 𝓧, 𝓨)
end


# implements || f-g ||_\infty
function discretizeL∞norm(f::Function, g::Function, xq::Vector{T})::T where T
    return maximum( abs(f(xq[i])-g(xq[i])) for i = 1:length(xq))
end

function discretizeL∞normparallel!(𝑒_xq::SharedArray{Float64},
                f::Function, g::Function, xq::Vector{T}) where T

    @assert length(𝑒_xq) == length(xq)

    @sync @distributed for i = 1:length(xq)
        𝑒_xq[i] = convert(Float64, f(xq[i])-g(xq[i]) )
    end

    k = argmax(abs.(𝑒_xq))
    return abs(𝑒_xq[k]), k, 𝑒_xq[k]
end

function discretizeL∞normparallel(f::Function, g::Function, xq::Vector{T}) where T

    𝑒_xq = SharedArray{Float64}(length(xq))

   return discretizeL∞normparallel!(𝑒_xq, f, g, xq)
end

function discretizeL∞normparallel(𝑒::Function, xq::Vector{T}) where T

    𝑒_xq = SharedArray{Float64}(length(xq))

    @sync @distributed for i = 1:length(xq)
        𝑒_xq[i] = convert(Float64, 𝑒(xq[i]) )
    end

    k = argmax(abs.(𝑒_xq))

    return abs(𝑒_xq[k]), k, 𝑒_xq[k]
end

### take 2.
function geth(  𝑤::Vector{T},
                wfunc::Function,
                dfunc::Function,
                𝓧)::T where T

    numerator = zero(T)
    denominator = zero(T)
    for i = 1:length(𝓧)
        k = i-1
        numerator += 𝑤[i]*dfunc(𝓧[i])
        denominator += (-1)^k*𝑤[i]/wfunc(𝓧[i])
    end
    out = numerator/denominator

    @assert isfinite(out)

    return out
end

function getp(x, wfunc::Function, dfunc::Function, 𝓧)
    𝑓 = get𝑓(wfunc,dfunc,𝓧)
    𝑤 = get𝑤scaled(𝓧)

    return wfunc(x)*Barycentric2nditp(x,𝑤,𝓧,𝑓)
end


function get𝑓(wfunc::Function, dfunc::Function, 𝓧::Vector{T}) where T

    𝑤 = get𝑤scaled(𝓧)
    h = geth(𝑤, wfunc, dfunc, 𝓧)

    return collect( dfunc(𝓧[i]) - (-1)^(i-1)*h/wfunc(𝓧[i]) for i = 1:length(𝓧))
end


# efficient version if 𝑤 and h were computed already.
function getp(x, 𝑤::Vector{T}, h::T, wfunc::Function, dfunc::Function, 𝓧::Vector{T}) where T

    𝑓 = collect( dfunc(𝓧[i]) - (-1)^(i-1)*h/wfunc(𝓧[i]) for i = 1:length(𝓧))

    return wfunc(x)*Barycentric2nditp(x,𝑤,𝓧,𝑓)
end

# even more efficient version, when 𝑓 is already computed.
function getp(x, 𝑓::Vector{T}, 𝑤::Vector{T}, h::T, wfunc::Function, dfunc::Function, 𝓧::Vector{T}) where T

    return wfunc(x)*Barycentric2nditp(x,𝑤,𝓧,𝑓)
end

function raisedcosinefreqrsp(x::T, β, 𝑇)::T where T
    return raisedcosinefreqrsp(x,convert(T,β),convert(T,𝑇))
end

function raisedcosinefreqrsp(x::T, β::T, 𝑇::T)::T where T
    two = one(T) + one(T)
    low = (one(T)-β)/(two*𝑇)
    high = (one(T)+β)/(two*𝑇)

    if abs(x) <= low
        return one(T)
    elseif low < abs(x) <= high
        return one(T)/two*( one(T) + cos(π*𝑇/β *(abs(x)-low)) )
    end

    return zero(T)
end


# generates L+2 positions in [-1,1]
function getinitialpositions(   L::Int,
                                candidate_multiple::Int)
    #
    N_candidates = candidate_multiple*L

    # generate canadiate positions. Take them as Chebyshev nodes.
    𝓧 = collect( Chebyshev2ndnode(k,N_candidates-1) for k = 0:N_candidates-1 )

    # select L positions.
    V = getunivariateChebyVmatrix(𝓧, L+1)

    B_unused, 𝓲 = maxvolviaqr(V)

    𝓧_chosen = 𝓧[𝓲][1:L+2]

    sort!(𝓧_chosen)

    return 𝓧_chosen, 𝓲
end
