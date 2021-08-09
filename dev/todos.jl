# TODO
# - make sure running julia script in REPL, ARGS is empty.
# - what's up with computeDTFTviaformula().
# - make sure plot exports to HTML via Plotly().
# - figure out the order number vs. taps vs. L.
#   2.2.3 of https://tel.archives-ouvertes.fr/tel-01447081/document
#   make sure the Type 1 is verified.

# - export KRSurrogate ASAP. commandline pass custom (small) image,
#   or pass function in specialized file. non-adaptive kernels.
#   focus on write up the inverse solve.

## Future plans
- types 2 to 4 FIR mini-max design algorithms.
- more variety of transition bands.
- allow command line passing of the target magnitude response function.
- FIR filters by windowing.
- very large order FIR filters as implemented by Richard Lyon's book, also featured here: https://www.embedded.com/designing-high-order-fir-filters/
