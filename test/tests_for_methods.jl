
function testnaiveimpulsersp(N_tests::Int)

    discrepancy = 0.0
    for n = 1:N_tests
        ω = rand()*900
        h = randn(5)
        Hn = h[1] + h[2]*exp(-im*ω) + h[3]*exp(-im*2*ω) + h[4]*exp(-im*3*ω) + h[5]*exp(-im*4*ω)
        Ht = naiveimpulsersp(h,ω)

        discrepancy += abs(Hn-Ht)
    end

    println("testnaiveimpulsersp: total discrepancy = ", discrepancy)
    return discrepancy
end

function testscaletoChebyinterval(N_tests::Int)

    for n = 1:N_tests
        a = abs(randn())
        b = abs(randn())
        α = rand()

        x = (1-α)*a+α*b
        y = scaletoChebyinterval(x,a,b)
        println(y)

        if !( -1 <= y <= 1)
            println("testscaletoChebyinterval failed")
            return false
        end
    end

    println("testscaletoChebyinterval passed")
    return true
end

function testscaletoChebyinterval2(N_tests::Int, N_test_locations::Int)

    function test(N_test_locations)
        discrepancy = 0.0
        for α in collect(LinRange(0.0,1.0,N_test_locations))
            x = (1-α)*a+α*b
            y = scaletoChebyinterval(x,a,b)
            xx = scalefromChebyinterval(y,a,b)

            discrepancy += abs(xx-x)
        end

        return discrepancy
    end

    discrepancy = 0.0
    for n = 1:N_tests
        a = abs(randn())
        b = abs(randn())

        discrepancy += test(a,b, N_test_locations)
    end

    return discrepancy
end


function testget𝑤scaled(𝓧::Vector{T})::Vector{T} where T

    w = get𝑤scaled(𝓧)
    w2 = get𝑤(𝓧)

    C = 0.5
    w1 = w/C^L

    return LinearAlgebra.norm(w1-w2)
end
