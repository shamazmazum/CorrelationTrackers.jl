# CorrelationTrackers.jl

`CorrelationTrackers.jl` package has a means to do fast updates of correlation
functions calculated by `CorrelationFunctions.Directional` module in
`CorrelationFunctions.jl`. Correlation functions are recalculated when you
change an element of the underlying array. Currently, only `Directional.l2` and
`Directional.s2` functions are supported.

An update after one element of the underlying system was changed is about 20000
times faster compared to a full recalculation of the correlation functions for
that system.

## Examples

This is a step-by-step guide on how to use `CorrelationTrackers.jl`:

1. For the beginning, create your system you want to calculate correlation
   functions for. Let it be two-dimensional two-phase system with size 20x30.
2. Pick combinations of correlation function + phase which you want to track,
   e.g $S_2^1(r)$, $L_2^1(r)$ and $L_2^0(r)$. A combination of function and
   phase is represented in `TrackerData` structure,
   e.g. `TrackerData(Directional.s2, 1)` will track $S_2^1(r)$ function.
3. Create `CorrelationTracker` structure and access it as an ordinary array.
4. Obtain correlation functions at any time using functions from
   `CorrelationFunctions.Directional`. module.

```@example
using StatsBase
using PrettyTables
using Random
using CorrelationFunctions
using CorrelationTrackers

# Create our system
system = rand(MersenneTwister(348), 0:1, (20, 30))

# Suppose we want to track $S_2^1(r)$ and $L_2^0(r)$.
tracking = [TrackedData(Directional.s2, 1), 
            TrackedData(Directional.l2, 1)]

# Create the tracker. It may require some time if your system is big.
tracker = CorrelationTracker(system; tracking = tracking)

# Do some phase swapping in both the original system and the tracker
system[12, 25]  = 1 - system[12, 25]
system[11, 20]  = 1 - system[11, 20]
tracker[12, 25] = 1 - tracker[12, 25]
tracker[11, 20] = 1 - tracker[11, 20]

pretty_table(stdout, hcat(Directional.s2(system,  1) |> mean,
                          Directional.l2(system,  1) |> mean,
                          Directional.s2(tracker, 1) |> mean,
                          Directional.l2(tracker, 1) |> mean);
                     header = ["S2 recalculated", "L2 recalculated",
                               "S2 tracked",      "L2 tracked"])
```

## API

```@docs
TrackedData
default_trackers
CorrelationTracker
tracked_data
tracked_length
tracked_directions
softupdate
```

## Caveats

Currently `Directional.l2` and `Directional.s2` return internal structures of
supplied tracker. Do not modify them. This behavior can be changed to a safer
one in the future.
