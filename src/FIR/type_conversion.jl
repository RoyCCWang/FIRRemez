# evaluate and understand the 4 types of linear phase FIR filter responses.

# from type 1, degree L to type 2, degree L.
# 𝑎 here is b_tilde from Flilip's thesis Section 2.2.3.
# a here is b from thesis.
function converttotype2(𝑎::Vector{T}) where T
    L = length(𝑎)

    a = Vector{T}(undef,L)

    a[1] = 𝑎[1] + 𝑎[2]/2
    for k = 2:L-1
        a[k] = (𝑎[k-1+1]+𝑎[k+1])/2
    end

    a[L] = 𝑎[L]/2

    return a
end

# 𝑎 is c_tilde from Flilip's thesis Section 2.2.3.
# a here is c from thesis.
function converttotype3(𝑎::Vector{T}) where T
    L = length(𝑎)
    @assert L >= 3

    a = Vector{T}(undef,L)

    a[1] = 𝑎[1] - 𝑎[3]/2
    for k = 2:L-2
        a[k] = (𝑎[k-1+1]-𝑎[k+1+1])/2
    end

    a[L-1] = 𝑎[L-1]/2
    a[L] = 𝑎[L]/2

    return a
end

# 𝑎 is d_tilde from Flilip's thesis Section 2.2.3.
# a here is d from thesis.
function converttotype4(𝑎::Vector{T}) where T
    L = length(𝑎)

    a = Vector{T}(undef,L)

    a[1] = 𝑎[1] - 𝑎[2]/2
    for k = 2:L-1
        a[k] = (𝑎[k-1+1]-𝑎[k+1])/2
    end

    a[L] = 𝑎[L]/2

    return a
end

function getfilterfuncs(type_num::Int)
    @assert 1 <= type_num <= 4

    getafunc = getatype1
    gethfunc = gethtype1

    if type_num == 2
        getafunc = getatype2
        gethfunc = gethtype2

    elseif type_num == 3
        getafunc = getatype3
        gethfunc = gethtype3

    elseif type_num == 4
        getafunc = getatype4
        gethfunc = gethtype4
    end

    return getafunc, gethfunc
end

#### type 1
function gethtype1(a::Vector{T}) where T
    L = length(a) -1
    N = 2*L+1

    h = Vector{T}(undef, N)
    fill!(h,-1.23)

    #tmp = Vector{T}(undef, L+1)
    h[L+1] = a[1]
    for k = 1:L
        h[L-k+1] = a[k+1]/2
    end

    h[L+2:end] = reverse(h[1:L])

    return h
end

function getatype1(h::Vector{T}) where T
    N = length(h)
    @assert isodd(length(h))

    L = div(N-1,2)

    a = Vector{T}(undef,L+1)
    fill!(a,-1.23)

    a[1] = h[L+1]
    for k = 1:L
        a[k+1] = 2*h[L-k+1]
    end

    return a
end

function evalHeqn211(h::Vector{T}, ω ) where T

    a = getatype1(h)
    L = length(a)-1

    return sum( a[k+1]*cos(k*ω) for k = 0:L)
end

#### type 2

function evalHeqn212(h::Vector{T}, ω ) where T

    a = getatype2(h)
    L = length(a)

    return sum( a[k]*cos((k-0.5)*ω) for k = 1:L)
end

function getatype2(h::Vector{T}) where T

    N = length(h)
    @assert iseven(N)

    L = div(N,2)

    a = Vector{T}(undef,L)
    fill!(a,-1.23)

    for k = 1:L
        a[k] = 2*h[L-k+1]
    end

    return a
end

function gethtype2(a::Vector{T}) where T
    L = length(a)
    N = 2*L

    h = Vector{T}(undef, N)
    fill!(h,-1.23)


    for k = 1:L
        h[L-k+1] = a[k]/2
    end

    h[L+1:end] = reverse(h[1:L])

    return h
end

#### type 3

function evalHeqn213(h::Vector{T}, ω ) where T

    a = getatype3(h)
    L = length(a)

    return sum( a[k]*sin(k*ω) for k = 1:L)
end

function gethtype3(a::Vector{T}) where T
    L = length(a)
    N = 2*L+1

    h = Vector{T}(undef, N)
    fill!(h,-1.23)

    #tmp = Vector{T}(undef, L+1)
    h[L+1] = zero(T)
    for k = 1:L
        h[L-k+1] = a[k]/2
    end

    h[L+2:end] = reverse(-h[1:L])

    return h
end

function getatype3(h::Vector{T}) where T
    N = length(h)
    @assert isodd(length(h))

    L = div(N-1,2)
    @assert isnumericallyclose(h[L+1],zero(T))

    a = Vector{T}(undef,L)
    fill!(a,-1.23)

    a[1] = zero(T)
    for k = 1:L
        a[k] = 2*h[L-k+1]
    end

    return a
end

#### type 4

function evalHeqn214(h::Vector{T}, ω ) where T

    a = getatype4(h)
    L = length(a)

    return sum( a[k]*sin((k-0.5)*ω) for k = 1:L )
end

function getatype4(h::Vector{T}) where T

    N = length(h)
    @assert iseven(N)

    L = div(N,2)

    a = Vector{T}(undef,L)
    fill!(a,-1.23)

    for k = 1:L
        a[k] = 2*h[L-k+1]
    end

    return a
end

function gethtype4(a::Vector{T}) where T
    L = length(a)
    N = 2*L

    h = Vector{T}(undef, N)
    fill!(h,-1.25)


    for k = 1:L
        h[L-k+1] = a[k]/2
    end

    h[L+1:end] = reverse(-h[1:L])

    return h
end

#### generic for all 4 types.

function evalHeqn210(h,ω,type_num)
    L = 0
    if type_num == 3 || type_num == 4
        # case negative symmetric.
        L = 1
    end

    N = length(h)
    G = naiveimpulsersp(h,ω)

    return G*exp( -im*(L*π - (N-1)*ω)/2 )
end


function naiveimpulsersp(h,ω)
    N = length(h)
    return sum( h[n+1]*exp(-im*ω*n) for n = 0:N-1 )
end
