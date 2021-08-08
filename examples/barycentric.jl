# 1-D interpolation from irregularly sampled data.

import PyPlot
import Random
import FIRRemez

fig_num = 1
PyPlot.close("all")

Random.seed!(25)

# generate sampling locations.
N_samples = 25
X = convert( Vector{BigFloat}, randn(N_samples).*π )

# generate samples.
Y = sin.(X)

# query.
N_queries = 500
xq_BF = convert(Vector{BigFloat}, collect( LinRange(minimum(X), maximum(X), N_queries) ) )
yq_BF = collect( FIRRemez.Barycentric2nditp(xq_BF[i],X,Y) for i = 1:length(xq_BF) )

xq = convert(Vector{Float64}, xq_BF)
yq = convert(Vector{Float64}, yq_BF)

X_display = convert(Vector{Float64}, X)
Y_display = convert(Vector{Float64}, Y)

PyPlot.figure(fig_num)
fig_num += 1
PyPlot.plot(xq, yq, label = "Barycentric")
PyPlot.plot(X_display, Y_display, label = "samples", "x")
PyPlot.title("Barycentric interpolation")
PyPlot.legend()
