
####### Routines related to Chebyshev interpolation or polynomials that are not about
#   evaluation nor interpolation set up.



# scale x ∈ [a,b] to [-1,1].
function scaletoChebyinterval(x::T, a::T, b::T)::T where T
    @assert a <= x <= b
    @assert b > a # disallow empty interval.

    Δ::T = b-a
    return (x-a)/Δ*convert(T,2) - one(T)
end

# scale y ∈ [-1,1] to [a,b].
function scalefromChebyinterval(y::T, a::T, b::T)::T where T
    @assert -one(T) <= y <= one(T)

    Δ::T = b-a
    return (y+1)/2*Δ + a
end


function Chebyshev2nd(x::T, L::Int)::T where T
    two = one(T) + one(T)

    out = zero(T)
    for l = 0:div(L,2)
        out += (-one(Int))^l*binomial(L-l,l)*(two*x)^(L-2*l)
    end

    return out
end


function Chebyshev1stall(x::T, L::Int)::Vector{T} where T
    @assert L >= 2

    # setup.
    two = one(T)+one(T)

    out = Vector{T}(undef,L+1)
    out[1] = one(T)
    out[2] = x
    for i = 3:length(out)
        out[i] = two*x*out[i-1] - out[i-2]
    end

    return out
end

function Chebyshev1st(x::T, L::Int)::T where T
    @assert L >= 0

    @assert isnumericallyin(x,-one(T),one(T))
    clamp(x,-one(T),one(T))

    return cos(convert(T,L)*acos(x))
end

# fill array with function
# equivalent to: collect( 𝝋[j](y[i]) for j = 1:n+1 for i = 1:m+1 ), but transposed.
function getunivariateChebyVmatrix(𝓧::Vector{T}, L::Int) where T
    p = length(𝓧)
    V = Matrix{T}(undef,L+1,p)

    for j = 1:p
        V[:,j] = Chebyshev1stall(𝓧[j],L)
    end

    return V
end

function Chebyshev1stnodes(Lp1::Int)::Vector{BigFloat}
    L = Lp1 - 1
    half = BigFloat("0.5")
    return collect( cos( (BigFloat(L-i)+half)*BigFloat(π)/BigFloat(L+1) ) for i = 0:L+1 )
end


function Chebyshev2ndnodes(L::Int)::Vector{BigFloat}
    return collect( cos( (L-i)*BigFloat(π)/BigFloat(L) ) for i = 0:L+1 )
end

function Chebyshev2ndnode(i::Int, L::Int)::BigFloat
    return cos( (L-i)*BigFloat(π)/BigFloat(L) )
end


# Colleague matrix for a L-degree Chebyshev polynomial of the mth kind.
function constructcolleaguematrix(c::Vector{T}, m::Int) where T
    @assert m == 1 || m == 2

    L = length(c) - 1

    while isnumericallyclose(c[L+1], zero(T), 1e-9) && L > 2
        L -= 1
    end

    # pre-compute.
    two = convert(T,2)
    half = one(T)/two

    C = zeros(T,L,L)
    C[2,1] = half

    # fill out the band.
    for l = 2:L-1
        C[l-1,l] = half
        C[l+1,l] = half
    end
    C[L-1,L] = half

    if m == 1
        C[1,2] = one(T)
    end

    # fill out bottom row.
    for l = 1:L
        C[L,l] = -c[l]/(two*c[L+1])
    end

    C[L,L-1] += half

    return C
end



# polynomial basis is Chebyshev polynomials of the m-th kind.
# coeffs are a_k for m = 1
# coeffs are c_k-1 = k*a_k for m = 2
function findChebyitproots(a::Vector{T}, m::Int)::Vector{T} where T
    if count(!isnumericallyclose(a[i],zero(T)) for i = 1:length(a)) < 2
        # case, a is the zero vector or a has less than 2 non-zero entries.
        #   The latter case means constant term and linear term. No local extrema.
        return Vector{T}(undef,0)
    end

    C = constructcolleaguematrix(a,m)

    return findChebyitproots(C)
end

function findChebyitproots(C::Matrix{T})::Vector{T} where T
    C_eig = LinearAlgebra.eigen(convert(Matrix{Float64},C))
    𝑧 = convert(Vector{T}, real.(C_eig.values))
    unique!(𝑧)
    filter!(xx->isnumericallyin(xx,-one(T),one(T)), 𝑧)
    clamp!(𝑧,-one(T),one(T))

    return 𝑧
end

function findChebyitpextrema(a::Vector{T}) where T
    L = length(a)-1

    c = collect( k*a[k+1] for k = 1:L)

    extrema_positions = findChebyitproots(c,2)

    return extrema_positions
end

# remove roots that do not have a zero derivative.
function findChebyitpextremawithchk(a::Vector{T}, tol::T = convert(T,0.2)) where T
    X = findChebyitpextrema(a)

    𝑑itp_X = evalunivariateChebyshevpolynomialderivative(a,X)

    remove_list = Vector{Int}(undef,0)
    for i = 1:length(X)
        if !isnumericallyclose(𝑑itp_X[i],zero(T),tol)
            push!(remove_list,i)
        end
    end

    deleteat!(X,remove_list)

    return X
end

# N is the number of samples to take between each reference position.
# 𝓧 is the set of reference positions.
function findChebyitpextremainterval(  𝑒::Function,
                                    L::Int,
                                    𝓟_in::Vector{T},
                                    tol_𝓧_spacing, # was defaulted to convert(T,1e-6/length(𝓟_in)
                                    tol_derivative_zero) where T


    ### get roots.
    sample_position_template = collect( Chebyshev2ndnode(k,L) for k = 0:L )
    out, d_set, X_set = getextremafromroots(𝑒, L, 𝓟_in,
                                            sample_position_template,
                                            tol_𝓧_spacing, tol_derivative_zero)

    # check each point in out to see if its a stationary point.
    chkfunc = xx->numericalchkextrema(𝑒, xx, tol_derivative_zero)
    filter!(chkfunc,out)

    sort!(out)
    unique!(out)

    return out, d_set, X_set
end


function discreteinnerprod(f::Function, g::Function, L::Int)

    ν = collect( Chebyshev2ndnode(k,L) for k = 0:L )

    term1 = 1/2 * f(ν[1]) * g(ν[1])

    term2 = BigFloat("0")
    for l = 1:L-1
        term2 += f(ν[l+1]) * g(ν[l+1])
    end

    term3 = 1/2 * f(ν[L+1]) * g(ν[L+1])

    return term1 + term2 + term3
end
