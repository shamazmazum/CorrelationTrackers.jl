module CorrelationTrackers
using CorrelationFunctions
using Base.Iterators: zip, countfrom, takewhile

include("slices.jl")
include("tracker.jl")

export CorrelationTracker

end # module
