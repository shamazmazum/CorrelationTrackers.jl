struct CorrelationTracker{T, N} <: AbstractArray{T, N}
    system   :: Array{T, N}
    phase    :: Int
    periodic :: Bool
    l2       :: Directional.CorrelationData
    s2       :: Directional.CorrelationData
end

function CorrelationTracker(system     :: AbstractArray,
                            phase      :: Int;
                            periodic   :: Bool           = false,
                            directions :: Vector{Symbol} = system |> Directional.default_directions,
                            kwargs...)
    l2 = Directional.l2(system, phase;
                        periodic   = periodic,
                        directions = directions,
                        kwargs...)
    s2 = Directional.s2(system, phase;
                        periodic   = periodic,
                        directions = directions,
                        kwargs...)
    return CorrelationTracker(copy(system), phase, periodic, l2, s2)
end

function update_corrfunc!(tracker  :: CorrelationTracker{T},
                          val      :: T,
                          corrfunc :: Function,
                          idx      :: Tuple) where T
    corrdata = corrfunc(tracker)
    len = length(corrdata)
    for (direction, _) in corrdata
        slice, slice_idx = get_slice(tracker.system, idx, direction)
        oldcorr = corrfunc(slice, tracker.phase;
                           periodic = tracker.periodic, len = len)
        slice[slice_idx] = val
        newcorr = corrfunc(slice, tracker.phase;
                           periodic = tracker.periodic, len = len)
        diff = newcorr.success[:x] .- oldcorr.success[:x]
        corrdata.success[direction] .+= diff
    end
end

# Correlation functions interface
# ! FIXME: It should be safe to return internal structures as long as
# noone is going to modify them. I see no such scenario.
Directional.l2(x :: CorrelationTracker) = x.l2
Directional.s2(x :: CorrelationTracker) = x.s2

# Array interface
Base.size(x :: CorrelationTracker) = size(x.system)
Base.getindex(x :: CorrelationTracker, idx :: Vararg{Int}) = getindex(x.system, idx...)
function Base.setindex!(x   :: CorrelationTracker{T},
                        val :: T,
                        idx :: Vararg{Int}) where T
    for corrfunc in (Directional.l2, Directional.s2)
        update_corrfunc!(x, val, corrfunc, idx)
    end
    x.system[idx...] = val
end
