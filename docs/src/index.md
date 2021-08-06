# CorrelationTrackers.jl

`CorrelationTrackers.jl` package has means to do fast updates of correlation
functions calculated by `CorrelationFunctions.Directional` module in
`CorrelationFunctions.jl`. Correlation functions are recalculated when you
change an element of the underlying array. Currently supported functions are:

* Two-point function $S_2(r)$
* Lineal-path function $L_2(r)$
* Surface-surface function $F_{ss}(r)$
* Surface-void function $F_{sv}(r)$

An update after one element of the underlying system was changed is about 10000
times faster compared to a full recalculation of the correlation functions for
that system.

This package also supports fast rollback to the previous state if the last
update must be rejected (via `AnnealingRollbackAPI`).

## Examples

This is a step-by-step guide on how to use `CorrelationTrackers.jl`:

1. For the beginning, create your system you want to calculate correlation
   functions for. Let it be two-dimensional two-phase system with size 20x30.
2. Pick combinations of correlation function + phase which you want to track,
   e.g $S_2^{(1)}(r)$, $L_2^{(1)}(r)$ and $L_2^{(0)}(r)$. A combination of
   function and phase is represented by `AbstractTracker` types,
   e.g. `S2Tracker(1)` will track $S_2^{(1)}(r)$ function.
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

# Suppose we want to track two-point function for phase 1 and lineal-path
# function for phase 0.
tracking = [S2Tracker(1), L2Tracker(0)]

# Create the tracker. It may require some time if your system is big.
tracker = CorrelationTracker(system; tracking = tracking)

# Do some phase flipping in both the original system and the tracker
indices = CartesianIndices(tracker)
state = MersenneTwister(123)
for n in 1:100
    idx = rand(state, indices)
    system[idx]  = 1 - system[idx]
    tracker[idx] = 1 - tracker[idx]
end

pretty_table(stdout, hcat(Directional.s2(system,  1) |> mean,
                          Directional.s2(tracker, 1) |> mean,
                          Directional.l2(system,  0) |> mean,
                          Directional.l2(tracker, 0) |> mean);
                     header = ["S2 recalculated", "S2 tracked",
                               "L2 recalculated", "L2 tracked"])
```

## API

### Types used to designate a correlation function

```@docs
AbstractTracker
S2Tracker
L2Tracker
SSTracker
SVTracker
```

### Constructor and accessors

```@docs
CorrelationTracker
default_trackers
tracked_data
tracked_length
tracked_directions
```

## Caveats

Currently the functions from `CorrelationFunctions.jl` package return internal
structures of supplied tracker. Do not modify returned values. This behavior can
be changed to a safer one in the future.

Phase argument to `SSTracker` and `SVTracker` is currently ignored. It should be
OK as long as tracked system is two-phase.
