## stress test clenshaw evaluation against naive evaluate of Chebyshev polynomials.

@everywhere using Distributed
@everywhere using SharedArrays

@everywhere import Random
@everywhere import LinearAlgebra
@everywhere import Printf

@everywhere include("../src/minimax/interpolators.jl")


@everywhere include("../src/misc/utilities.jl")

Random.seed!(25)

fig_num = 1

N = 1000 # number of positions used in each test.
N_tests = 1 # number of tests.

function stresstestget𝑤scaledfast(N_tests, N)
    discrepancy = SharedArray{Float64}(N_tests)

    @sync @distributed for n = 1:N_tests
        X = collect( convert(BigFloat,rand()*2-1) for i = 1:N )

        w = get𝑤scaled(X)
        w2 = get𝑤scaledfast(X)

        discrepancy[n] = LinearAlgebra.norm(w-w2)
    end

    return convert(Vector{Float64},discrepancy)
end

discrepancy = stresstestget𝑤scaledfast(N_tests, N)
println("norm(discrepancy) = ", LinearAlgebra.norm(discrepancy))
println()


#### speed comparison.
X = collect( convert(BigFloat,rand()*2-1) for i = 1:N )

w = get𝑤scaled(X)
w2 = get𝑤scaledfast(X)

d = LinearAlgebra.norm(w-w2)

println("discrepancy between naive and parallel: ", d)
