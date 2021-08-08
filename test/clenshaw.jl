## stress test clenshaw evaluation against naive evaluate of Chebyshev polynomials.

@everywhere using Distributed
@everywhere using SharedArrays

@everywhere import Random
@everywhere import LinearAlgebra
@everywhere import Printf

@everywhere include("../src/minimax/Chebyshev.jl")
@everywhere include("../src/minimax/eval.jl")

@everywhere include("../src/misc/utilities.jl")

Random.seed!(25)

fig_num = 1

N = 10000 # number of positions used in each test.
N_tests = 10 # number of tests.
xq = convert(Vector{BigFloat}, collect(LinRange(-1,1,N_tests)) )

function stresstestclenshaw(xq, N)
    discrepancy = SharedArray{Float64}(length(xq))

    @sync @distributed for n = 1:length(xq)
        c = convert(Vector{BigFloat},randn(N))

        f_itp_Cheby = xx->evalunivariateChebyshevpolynomial(c,xx)
        f_itp_Cheby2 = xx->clenshaweval(c,xx)

        y = f_itp_Cheby.(xq)
        y2 = f_itp_Cheby2.(xq)

        discrepancy[n] = LinearAlgebra.norm(y-y2)
    end

    return convert(Vector{Float64},discrepancy)
end

discrepancy = stresstestclenshaw(xq, N)
println("norm(discrepancy) = ", LinearAlgebra.norm(discrepancy))
println()

#### speed comparison.
c = convert(Vector{BigFloat},randn(N))

f_itp_Cheby = xx->evalunivariateChebyshevpolynomial(c,xx)
f_itp_Cheby2 = xx->clenshaweval(c,xx)
f_itp_Cheby3 = xx->evalunivariateChebyshevpolynomialparallel(c,xx)


println("naive")
@time y = f_itp_Cheby.(xq)

println("clenshaw")
@time y2 = f_itp_Cheby2.(xq)

println("naive parallel")
@time y3 = f_itp_Cheby3.(xq)
println()

println("discrepancy between naive and clenshaw: ", LinearAlgebra.norm(y-y2))
println("discrepancy between naive and naive parallel: ", LinearAlgebra.norm(y-y3))


# naive parallel is slower for some reason.
