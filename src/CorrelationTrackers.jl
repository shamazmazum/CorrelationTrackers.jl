module CorrelationTrackers
using CorrelationFunctions
using Base.Iterators: zip, countfrom, takewhile, take
using CircularArrays: CircularArray

include("slices.jl")
include("tracker.jl")

export
    CorrelationTracker,
    TrackedData,
    default_trackers,

    tracked_data,
    tracked_length,
    tracked_directions

end # module
