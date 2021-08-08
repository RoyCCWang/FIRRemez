### visualize a solution.


function visualizefiltersolution(𝓧, xq, wfunc, dfunc,f)
    𝑤_new = get𝑤scaled(𝓧)
    h_new = geth(𝑤_new, wfunc, dfunc, 𝓧)

    f_itp = xx->getp(xx, 𝑤_new, h_new, wfunc, dfunc, 𝓧)
    f_itp_xq = collect( f_itp(xq[i]) for i = 1:length(xq) )
    f_itp_𝓧 = collect( f_itp(𝓧[i]) for i = 1:length(𝓧) )

    𝑒 = xx->(f(xx)-f_itp(xx))
    𝑒_xq = collect( 𝑒(xq[i]) for i = 1:length(xq) )
    𝑒_𝓧 = collect( 𝑒(𝓧[i]) for i = 1:length(𝓧) )

    return f_itp_xq, f_itp_𝓧, 𝑒_xq, 𝑒_𝓧
end
