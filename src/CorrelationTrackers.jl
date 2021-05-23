module CorrelationTrackers
using CorrelationFunctions
using Base.Iterators: zip, countfrom, takewhile, take
using CircularArrays: CircularArray

include("slices.jl")
include("tracker.jl")

export CorrelationTracker

end # module
