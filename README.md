# FIRRemez

[![Build Status](https://github.com/RCCWang/FIRRemez.jl/workflows/CI/badge.svg)](https://github.com/RCCWang/FIRRemez.jl/actions)
[![Coverage](https://codecov.io/gh/RCCWang/FIRRemez.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RCCWang/FIRRemez.jl)

This Julia package lets one design type-1 mini-max FIR filters for large orders (up to order 5000 was tested). It uses local machine parallel processing (multiple cores). The current algorithm was taken from chapter 3 of Filip's 2016 thesis "Robust tools for weighted Chebyshev approximation and applications to digital filter design".

## To install this package
Install Julia, clone this repo to a local folder, say, `/home/repo/FIRRemez`

Launch Julia REPL, then press `]` on your computer keyboard.
Type `add /home/repo/FIRRemez`
Press the backspace key on your computer keyboard.

Examples are in the `./examples` folder.

## Example usage for type-1 FIR mini-max lowpass filter design
Make sure the XLSX.jl and BSON.jl packages are installed. They are used to save the resultant coefficients into a MS Excel file and binary JSON format in the `./output` folder.

Navigate to the `./examples` folder using the linux terminal. Here is an example of the syntax for running `template_lowpass.jl`:
To use 15 parallel processes to solve for a mini-max lowpass filter with cosine transition, order 200, pass band 0.025π rad, stop band 0.05π rad, run `julia -p 15 template_lowpass.jl 200 0.07853981633974483 0.15707963267948966`
This command took a few minutes to run on a 16-core AMD ThreadRipper 1950X system.

Informally speaking, the number of processes should be less than the number of physical double-precision floating-point cores on your CPU.
