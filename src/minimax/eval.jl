
##### This file contain routines related to evaluating polynomials or interpolators.

# the second (more stable) Barycentric interpolation formula.
function Barycentric2nditp(x::T, ð§, ð¨::Vector{T})::T where T
    w = getð¤(ð§)
    return Barycentric2nditp(x,w,ð§,ð¨)
end

function Barycentric2nditp(x::T, ð¤::Vector{T}, ð§, ð¨::Vector{T})::T where T

    tmp_numerator = zero(T)
    tmp_denominator = zero(T)
    for k = 1:length(ð§)
        tmp = x-ð§[k]
        if isnumericallyclose(tmp,zero(T))
            return ð¨[k]
        end

        tmp_numerator += ð¨[k]*ð¤[k]/(x-ð§[k])
        tmp_denominator += ð¤[k]/(x-ð§[k])
    end

    @assert tmp_denominator != zero(T)
    return tmp_numerator/tmp_denominator
end




##### other evaluations.
function evalunivariateLagrangepolynomial(  x::T,
                                            ð§::Vector{T},
                                            a::Vector{T})::T where T
    out = zero(T)
    for k = 1:length(a)
        out += a[k]*getkthLangrangebasis(x,ð§,k)
    end

    return out
end


function getkthLangrangebasis(x::T, ð§::Vector{T}, k::Int)::T where T
    Lp1 = length(ð§)

    numerator = one(T)
    denominator = one(T)
    for i = 1:Lp1
        if i != k
            numerator *= x-ð§[i]
            denominator *= ð§[k]-ð§[i]
        end
    end

    return numerator/denominator
end


##### routines related to the evaluation of Chebyshev polynomials at nodes.


function evalunivariateChebyshevpolynomial(c::Vector{T},
                                            X::Vector{T}) where T

    return collect( evalunivariateChebyshevpolynomial(c,x) for x in X )
end

# # naive implementation.
function evalunivariateChebyshevpolynomial( c::Vector{T},
                                            x::T) where T
    L = length(c) - 1

    out = zero(T)
    for k = 0:L
        out += c[k+1]*Chebyshev1st(x,k)
    end

    return out
end

function evalunivariateChebyshevpolynomialparallel( c::Vector{T},
                                                    x::T) where T
    L = length(c) - 1

    out = @distributed (+) for k = 0:L
        c[k+1]*Chebyshev1st(x,k)
    end

    return out
end


# Clenshaw implementation. does evalunivariateChebyshevpolynomial()
function clenshaweval(c::Vector{T},x::T) where T
    two = one(T) + one(T)

    L = length(c) - 1

    b_kp2 = zero(T)
    b_kp1 = c[end]

    for k = L-1:-1:1
        b_k = two*x*b_kp1 - b_kp2 + c[k+1]

        # update
        b_kp2 = b_kp1
        b_kp1 = b_k
    end

    return c[1] +x*b_kp1 - b_kp2
end

function getChebyshevinterpolator(  f::Function,
                                    L::Int)
    Lp1 = L + 1
    Î½ = collect( Chebyshev2ndnode(k,L) for k = 0:L )

    return getChebyshevinterpolator(f,L,Î½,Î½)
end

function getChebyshevinterpolator(  f::Function,
                                    L::Int, a::T, b::T) where T<:Real
    Lp1 = L + 1
    Î½ = collect( convert(T,Chebyshev2ndnode(k,L)) for k = 0:L )

    # rescale nodes to [a,b].
    Î½f = collect( scalefromChebyinterval(Î½[i],a,b) for i = 1:length(Î½) )

    return getChebyshevinterpolator(f,L,Î½f,Î½)
end

# uses Float64 SharedArrays.
function getChebyshevinterpolatorfast(  f::Function,
                                    L::Int)
    Lp1 = L + 1
    Î½ = SharedArray{Float64}(L+1)
    @sync @distributed for k = 0:L
        Î½[k+1] = convert(Float64, Chebyshev2ndnode(k,L))
    end

    return getChebyshevinterpolatorfast(f,L,Î½,Î½)
end

function getChebyshevinterpolatorfast(  f::Function,
                                    L::Int, a_inp, b_inp)
    Lp1 = L + 1
    Î½ = SharedArray{Float64}(L+1)
    @sync @distributed for k = 0:L
        Î½[k+1] = convert(Float64, Chebyshev2ndnode(k,L))
    end

    # rescale nodes to [a,b].
    a = convert(T,a_inp)
    b = convert(T,b_inp)

    Î½f = SharedArray{Float64}(L+1)
    @sync @distributed for i = 1:length(Î½)
        Î½f[i] = convert(Float64, scalefromChebyinterval(Î½[i],a,b))
    end

    return getChebyshevinterpolatorfast(f,L,Î½f,Î½)
end

# # naive implementation.
function evalunivariateChebyshevpolynomial( c::Vector{T},
                                            x::T,
                                            a::T,
                                            b::T) where T
    L = length(c) - 1

    v = scaletoChebyinterval(x,a,b)

    out = zero(T)
    for k = 0:L
        out += c[k+1]*Chebyshev1st(v,k)
    end

    return out
end


function evalunivariateChebyshevpolynomialderivative(a::Vector{T}, x::T) where T
    L = length(a) - 1
    c = collect( k*a[k+1] for k = 1:L)

    out = zero(T)
    for i = 0:L-1
        out += c[i+1]*Chebyshev2nd(x, i)
    end

    return out
end



function evalunivariateChebyshevpolynomialderivative(c::Vector{T}, X::Vector{T}) where T
    return collect( evalunivariateChebyshevpolynomialderivative(c,x) for x in X )
end



function evalChebyshevpolynomialinterval(   x::T,
                                            d_set::Vector{Vector{T}},
                                            X_set::Vector{Vector{T}}) where T

    for i = 1:length(X_set)
        X = X_set[i]
        if X[1] <= x <= X[end]
            return evalunivariateChebyshevpolynomial(d_set[i], x, X[1], X[end])
        end
    end

    return convert(T,-Inf) # when x isn't in the intervals specified by X_set
end



##### others.

# integer exponents.
# c has the format: [c_0; c_1; ...; c_{L}]
function evalunivariatepolynomial(c::Vector{T}, x::T)::T where T
    L = length(c) - 1

    out = c[1]
    for l = 1:L
        out += c[l+1]x^l
    end

    return out
end
