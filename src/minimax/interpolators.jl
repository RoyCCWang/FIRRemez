
########## Routines related to setting up interpolators.

##### Routines related to Barycentric interpolation.

# Barycentric Interpolations
function Barycentric𝓵(x::Real, 𝓧)
    return prod( x-𝓧[i] for i = 1:length(𝓧))
end

# This is l' in thesis (see Eq. 3.15).
function Barycentric𝓵prime(k::Int, 𝓧)
    return prod( 𝓧[k]-𝓧[i] for i = 1:length(𝓧) if i != k)
end

function get𝑤(k::Int, 𝓧::Vector{T})::T where T
    return one(T)/Barycentric𝓵prime(k, 𝓧)
end

function get𝑤(𝓧::Vector{T})::Vector{T} where T
    return collect( get𝑤(k,𝓧) for k = 1:length(𝓧) )
end

function get𝑤scaled!(w::Vector{T}, 𝓧::Vector{T})::Nothing where T
    L = length(𝓧) - 2
    factor::T = L*log(one(T)+one(T)) # C is 0.5 for the interval [-1,1]

    resize!(w,length(𝓧))
    for k = 0:length(𝓧)-1

        # convert to Int to speed up the multiplication used in prod().
        numerator = prod( convert(Int,sign(𝓧[k+1]-𝓧[i])) for i = 1:length(𝓧) if i != k+1)

        # log-sum-exp trick for adding ill-conditioned numbers.
        denominator = exp( factor + sum( log(abs(𝓧[k+1]-𝓧[i])) for i = 1:length(𝓧) if i != k+1 ) )

        w[k+1] = numerator/denominator
        if !isfinite(w[k+1])
            println("𝓧 is ", 𝓧)
            println("k is ", k)
            println("𝓧[k] is ", 𝓧[k])
            println("numerator is ", numerator)
            println("denominator is ", denominator)
        end
        @assert isfinite(w[k+1])
    end

    return nothing
end

function get𝑤scaled(𝓧::Vector{T})::Vector{T} where T

    w = Vector{T}(undef,length(𝓧))
    get𝑤scaled!(w,𝓧)

    return w
end

# uses SharedArray{Float64} for internal computation.
function get𝑤scaledfast(𝓧::Vector{T})::Vector{T} where T

    L = length(𝓧) - 2
    factor::T = L*log(one(T)+one(T)) # C is 0.5 for the interval [-1,1]

    w = SharedArray{Float64}(length(𝓧))
    @sync @distributed for k = 0:length(𝓧)-1

        numerator::T = prod( sign(𝓧[k+1]-𝓧[i]) for i = 1:length(𝓧) if i != k+1)
        denominator::T = exp( factor + sum( log(abs(𝓧[k+1]-𝓧[i])) for i = 1:length(𝓧) if i != k+1 ) )

        w[k+1] = convert(Float64,numerator/denominator)
    end

    # sanity check.
    for k = 0:length(𝓧)-1
        if !isfinite(w[k+1])
            println("𝓧 is ", 𝓧)
            println("k is ", k)
            println("𝓧[k] is ", 𝓧[k])
            println("numerator is ", numerator)
            println("denominator is ", denominator)
        end
        @assert isfinite(w[k+1])
    end

    return convert(Vector{T},w)
end

##### Routines related to Chebyshev interpolation at nodes.

function getChebyshevinterpolator(  f::Function,
                                    L::Int,
                                    νf::Vector{T},
                                    ν::Vector{T})::Tuple{Vector{T},Vector{T}} where T
    Lp1 = L + 1

    a = Vector{BigFloat}(undef,Lp1)
    for k = 0:L

        term1 = 1/L * f(νf[1]) * Chebyshev1st(ν[1],k)

        term2 = BigFloat("0")
        for l = 1:L-1
            term2 += f(νf[l+1]) * Chebyshev1st(ν[l+1],k)
        end
        term2 *= 2/L

        term3 = 1/L * f(νf[Lp1]) * Chebyshev1st(ν[Lp1],k)

        a[k+1] = term1 + term2 + term3
    end
    a[1] = a[1]/BigFloat(2)

    return a, ν
end

# efficient version. Less accurate.
function getChebyshevinterpolatorfast(  f::Function,
                                        L::Int,
                                        νf::SharedArray{Float64,1},
                                        ν::SharedArray{Float64,1})::Tuple{SharedArray{Float64,1},SharedArray{Float64,1}}
    Lp1 = L + 1

    ### pre-compute.

    # f(νf), length νf.
    f_νf = SharedArray{Float64}(length(νf))
    @sync @distributed for i = 1:length(f_νf)
        f_νf[i] = convert(Float64, f( convert(BigFloat,νf[i]) ))
    end

    ###

    a = SharedArray{Float64}(Lp1)
    for k = 0:L

        term1 = 1/L * f_νf[1] * Chebyshev1st(ν[1],k)

        term2 = zero(Float64)
        for l = 1:L-1
            term2 += f_νf[l+1] * Chebyshev1st(ν[l+1],k)
        end
        term2 *= 2/L

        term3 = 1/L * f_νf[Lp1] * Chebyshev1st(ν[Lp1],k)

        a[k+1] = term1 + term2 + term3
    end
    a[1] = a[1]/2

    return a, ν
end


# Barycentric interpolation coefficients to Chebyshev interpolation coefficients.
# Used in the minimax solution's conversion to impulse response in packagefiltersolution().
function BarycentrictoChebyshev(𝑋, wfunc, dfunc, type_num)
    𝑤 = get𝑤scaled(𝑋)
    h = geth(𝑤, wfunc, dfunc, 𝑋)
    𝑓 = collect( dfunc(𝑋[i]) - (-1)^(i-1)*h/wfunc(𝑋[i]) for i = 1:length(𝑋))

    f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)
    if type_num == 2
        f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)/type2weightfunc(xx)

    elseif type_num == 3
        f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)/type3weightfunc(xx)

    elseif type_num == 4
        f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)/type4weightfunc(xx)

    end

    # get Chebyshev coefficients.
    c, ν = getChebyshevinterpolator(f_itp,length(𝑋)-1)

    # sanity check.
    f_itp_Cheby = xx->evalunivariateChebyshevpolynomial(c,xx)

    xq = convert(Vector{BigFloat}, collect(LinRange(-1,1,200)) )
    yq = f_itp.(xq)
    yq2 = f_itp_Cheby.(xq)
    discrepancy = LinearAlgebra.norm(yq-yq2)
    return c, ν, discrepancy, xq, yq, yq2
end

# does not use BigFloats.
function BarycentrictoChebyshevfast(𝑋, wfunc, dfunc, type_num)

    𝑤 = get𝑤scaled(𝑋)
    h = geth(𝑤, wfunc, dfunc, 𝑋)

    𝑓 = collect( dfunc(𝑋[i]) - (-1)^(i-1)*h/wfunc(𝑋[i]) for i = 1:length(𝑋))

    f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)
    if type_num == 2
        f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)/type2weightfunc(xx)

    elseif type_num == 3
        f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)/type3weightfunc(xx)

    elseif type_num == 4
        f_itp = xx->getp(xx, 𝑓, 𝑤, h, wfunc, dfunc, 𝑋)/type4weightfunc(xx)

    end

    # get Chebyshev coefficients.
    c_SharedArray, ν_SharedArray = getChebyshevinterpolatorfast(f_itp,length(𝑋)-1)

    c = convert(Vector{BigFloat},c_SharedArray)
    ν = convert(Vector{BigFloat},ν_SharedArray)

    # sanity check.
    f_itp_Cheby = xx->clenshaweval(c,xx)

    xq = convert( Vector{BigFloat}, collect(LinRange(-1,1,200)) )

    yq = f_itp.(xq)
    yq2 = f_itp_Cheby.(xq)

    discrepancy = LinearAlgebra.norm(yq-yq2)
    return c, ν, discrepancy, xq, yq, yq2
end
