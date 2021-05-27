# CorrelationTrackers.jl
[![CI](https://github.com/shamazmazum/CorrelationTrackers.jl/actions/workflows/test.yml/badge.svg)](https://github.com/shamazmazum/CorrelationTrackers.jl/actions/workflows/test.yml)

`CorrelationTrackers.jl` is a library which helps you to fastly update correlation
functions calculated by `CorrelationFunctions.jl` package when you
insignificantly change the input.

To build a documentation do the following:

1. From Julia REPL: `import Pkg; Pkg.add("Documenter")`
2. From shell, this directory being the working directory: `cd docs && julia make.jl`

Also the documentation for the most recent release is available on
[GitHub Pages](https://shamazmazum.github.io/CorrelationTrackers.jl/v0.1.0/).

## Acknowledgements

Correlation maps are written by Alexey Samarin. Special thanks to Kirill Gerke,
Efim Lavrukhin and  Marina Karsanina for their support.
