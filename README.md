# FIRRemez

[![Build Status](https://github.com/RCCWang/FIRRemez.jl/workflows/CI/badge.svg)](https://github.com/RCCWang/FIRRemez.jl/actions)
[![Coverage](https://codecov.io/gh/RCCWang/FIRRemez.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RCCWang/FIRRemez.jl)

Clone this repo to a local folder, say, `/home/repo/FIRRemez`

Launch Julia REPL, then press `]` on your computer keyboard.
Type `add /home/repo/FIRRemez`
Press the backspace key on your computer keyboard.

Examples are in the `./examples` folder.

Example usage for lowpass filter design.
To use 11 parallel processes to solve for a mini-max lowpass filter with cosine transition, order 200, pass band 0.025π rad, stop band 0.05π rad, run `julia template_lowpass.jl -p 11 200 0.07853981633974483 0.15707963267948966`
