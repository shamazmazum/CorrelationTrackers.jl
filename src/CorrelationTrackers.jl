module CorrelationTrackers
using CorrelationFunctions
using Base.Iterators: zip, countfrom, takewhile, take
using CircularArrays: CircularArray

include("slices.jl")
include("tracker.jl")

export
    CorrelationTracker,
    TrackedData,
    tracked_data,
    tracking_by_default

end # module
