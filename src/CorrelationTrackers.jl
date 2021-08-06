module CorrelationTrackers
using Reexport: @reexport
@reexport using AnnealingAPI
using CorrelationFunctions
using Base.Iterators: zip, countfrom, takewhile, take
using CircularArrays: CircularArray
using LinearAlgebra: norm
using Images: imgradients, KernelFactors

include("slices.jl")
include("tracker.jl")
include("helpers.jl")

export
    CorrelationTracker,
    default_trackers
end # module
