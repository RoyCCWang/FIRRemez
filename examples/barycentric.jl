# 1-D interpolation from irregularly sampled data.

import PyPlot
import Random
import FIRRemez

fig_num = 1
PyPlot.close("all")

Random.seed!(25)

# generate sampling locations.
N_samples = 25
𝓧 = convert( Vector{BigFloat}, randn(N_samples).*π )

# generate samples.
𝓨 = sin.(𝓧)

# query.
N_queries = 500
xq = convert(Vector{BigFloat}, collect( LinRange(minimum(𝓧), maximum(𝓧), N_queries) ) )
yq = collect( FIRRemez.Barycentric2nditp(xq[i],𝓧,𝓨) for i = 1:length(xq) )

PyPlot.figure(fig_num)
fig_num += 1
PyPlot.plot(xq, yq, label = "Barycentric")
PyPlot.plot(𝓧, 𝓨, label = "samples", "x")
PyPlot.title("Barycentric interpolation")
PyPlot.legend()
