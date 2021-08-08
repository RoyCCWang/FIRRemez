import os
import h5py

filename = "filter_passband_0.314159_stopband_0.628319_25.jld"
f = h5py.File(filename, 'r')

h = f['h'].value
X = f['X'].value
passband = f['passband'].value
stopband = f['stopband'].value
