# CorrelationTrackers.jl

`CorrelationTrackers.jl` package contains a single `CorrelationTracker` type
which tracks an array and helps to recalculate correlation functions from
`CorrelationFunctions.jl` package in a fast way when the content of that array
is changed.

Currently `CorrelationTracker` can recalculate two-points(`s2`) and
lineal-path(`l2`) functions from `CorrelationFunctrions.Directional` module.

## Examples

``` jldoctest
julia> using Random

julia> using CorrelationTrackers

julia> using CorrelationFunctions

julia> system = rand(MersenneTwister(1), 0:1, (10,15))
10×15 Matrix{Int64}:
 1  0  0  0  0  1  1  0  0  0  0  0  1  1  1
 1  1  1  1  1  0  0  1  0  0  1  0  0  1  1
 0  0  0  0  1  1  1  0  1  1  1  1  1  0  1
 0  1  1  1  0  1  0  1  0  1  0  0  0  1  1
 0  0  1  0  1  0  0  1  1  0  0  1  1  0  0
 1  1  1  1  0  1  1  1  0  1  0  0  0  0  1
 1  0  0  0  1  1  0  1  1  0  1  1  0  0  0
 0  1  0  1  1  1  1  1  1  1  1  0  0  0  0
 1  1  1  0  0  0  0  1  0  1  0  1  1  1  1
 0  1  0  1  1  1  1  1  0  0  1  1  1  0  0

julia> tracker = CorrelationTracker(system, 1; len = 10)
10×15 CorrelationTracker{Int64, 2}:
 1  0  0  0  0  1  1  0  0  0  0  0  1  1  1
 1  1  1  1  1  0  0  1  0  0  1  0  0  1  1
 0  0  0  0  1  1  1  0  1  1  1  1  1  0  1
 0  1  1  1  0  1  0  1  0  1  0  0  0  1  1
 0  0  1  0  1  0  0  1  1  0  0  1  1  0  0
 1  1  1  1  0  1  1  1  0  1  0  0  0  0  1
 1  0  0  0  1  1  0  1  1  0  1  1  0  0  0
 0  1  0  1  1  1  1  1  1  1  1  0  0  0  0
 1  1  1  0  0  0  0  1  0  1  0  1  1  1  1
 0  1  0  1  1  1  1  1  0  0  1  1  1  0  0

julia> tracker[2,3] = 1 - tracker[2,3]
0

julia> Directional.l2(tracker)
┌───────────┬───────────┐
│         x │         y │
├───────────┼───────────┤
│  0.533333 │  0.533333 │
│  0.207407 │  0.307143 │
│ 0.0833333 │  0.169231 │
│  0.047619 │ 0.0916667 │
│ 0.0333333 │ 0.0545455 │
│ 0.0266667 │      0.03 │
│ 0.0166667 │ 0.0222222 │
│       0.0 │    0.0125 │
│       0.0 │       0.0 │
│       0.0 │       0.0 │
└───────────┴───────────┘


julia> Directional.s2(tracker)
┌──────────┬──────────┐
│        x │        y │
├──────────┼──────────┤
│ 0.533333 │ 0.533333 │
│ 0.207407 │ 0.307143 │
│ 0.333333 │ 0.284615 │
│ 0.238095 │    0.275 │
│ 0.277778 │      0.3 │
│     0.32 │     0.26 │
│     0.25 │      0.3 │
│ 0.288889 │   0.2875 │
│      0.3 │      0.2 │
│      0.2 │ 0.216667 │
└──────────┴──────────┘

```

You can use any arguments which you supply to `l2` or `s2` in
`CorrelationTracker` constructor.

## Caveats

Currently `l2` and `s2` return internal structures of supplied tracker. Do not
modify them. This behavior can be changed to a safer one in the future.
