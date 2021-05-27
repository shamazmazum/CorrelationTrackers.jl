# Tracked data is a pair of a function and a phase
struct TrackedData{T}
    func  :: Function
    phase :: T
end

struct CorrelationTracker{T, N} <: AbstractArray{T, N}
    system   :: Array{T, N}
    periodic :: Bool
    corrdata :: Dict{TrackedData{T}, Directional.CorrelationData}
end

function CorrelationTracker{T, N}(system     :: AbstractArray{T, N},
                                  tracking   :: Vector{TrackedData{T}};
                                  periodic   :: Bool = false,
                                  directions :: Vector{Symbol} =
                                      system |> Directional.default_directions,
                                  kwargs...) where {T, N}
    corrdata = Dict(data => data.func(system, data.phase;
                                      periodic   = periodic,
                                      directions = directions,
                                      kwargs...)
                    for data in tracking)
    return CorrelationTracker(copy(system), periodic, corrdata)
end

function update_corrfunc!(tracker  :: CorrelationTracker{T, N},
                          data     :: TrackedData{T},
                          val,
                          idx      :: Tuple) where {T, N}
    corrdata = tracker.corrdata[data]
    corrfunc = data.func
    phase    = data.phase
    len = length(corrdata)
    for direction in Directional.directions(corrdata)
        slice, slice_idx = get_slice(tracker.system,
                                     tracker.periodic,
                                     idx, direction)
        oldcorr = corrfunc(slice, phase;
                           periodic = tracker.periodic, len = len)
        slice[slice_idx] = val
        newcorr = corrfunc(slice, phase;
                           periodic = tracker.periodic, len = len)
        diff = newcorr.success[:x] .- oldcorr.success[:x]
        corrdata.success[direction] .+= diff
    end
end

"""
    tracked_data(x :: CorrelationTracker)

Return a vector of tracked correlation functions and phases.
"""
tracked_data(x :: CorrelationTracker) = x.corrdata |> keys |> collect

# Correlation functions interface
# ! FIXME: It should be safe to return internal structures as long as
# noone is going to modify them. I see no such scenario.
Directional.l2(x :: CorrelationTracker, phase) =
    x.corrdata[TrackedData(Directional.l2, phase)]
Directional.s2(x :: CorrelationTracker, phase) =
    x.corrdata[TrackedData(Directional.s2, phase)]

# Array interface
Base.size(x :: CorrelationTracker) = size(x.system)
Base.getindex(x :: CorrelationTracker, idx :: Vararg{Int}) = getindex(x.system, idx...)
function Base.setindex!(x   :: CorrelationTracker,
                        val,
                        idx :: Vararg{Int})
    for tracked_data in keys(x.corrdata)
        update_corrfunc!(x, tracked_data, val, idx)
    end
    x.system[idx...] = val
end
