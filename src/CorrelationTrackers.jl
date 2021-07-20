module CorrelationTrackers
using CorrelationFunctions
using Base.Iterators: zip, countfrom, takewhile, take
using CircularArrays: CircularArray
using LinearAlgebra: norm
using Images: imgradients, KernelFactors

include("slices.jl")
include("updaters_early.jl")
include("tracker.jl")
include("updaters_late.jl")

export
    CorrelationTracker,
    S2Tracker, L2Tracker,
    SSTracker,# SVTracker,
    AbstractTracker,
    default_trackers,

    tracked_data,
    tracked_length,
    tracked_directions
end # module
