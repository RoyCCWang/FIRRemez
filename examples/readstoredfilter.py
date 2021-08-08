import os
import h5py

filename = "./output/type1_lowpass_L200_passband_0.049530_stopband_0.068280.jld"
f = h5py.File(filename, 'r')

h = f['h'].value
X = f['X'].value
passband = f['passband'].value
stopband = f['stopband'].value
