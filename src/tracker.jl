struct CorrelationTracker{T, N} <: AbstractArray{T, N}
    system   :: Array{T, N}
    phase    :: Int
    periodic :: Bool
    corrdata :: Dict{Function, Directional.CorrelationData}
end

function CorrelationTracker(system     :: AbstractArray,
                            phase      :: Int;
                            periodic   :: Bool             = false,
                            directions :: Vector{Symbol}   = system |> Directional.default_directions,
                            functions  :: Vector{Function} = [Directional.s2, Directional.l2],
                            kwargs...)
    corrdata = Dict(func => func(system, phase;
                                 periodic   = periodic,
                                 directions = directions,
                                 kwargs...)
                    for func in functions)
    return CorrelationTracker(copy(system), phase, periodic, corrdata)
end

function update_corrfunc!(tracker  :: CorrelationTracker,
                          val,
                          corrfunc :: Function,
                          idx      :: Tuple)
    corrdata = tracker.corrdata[corrfunc]
    len = length(corrdata)
    for (direction, _) in corrdata
        slice, slice_idx = get_slice(tracker.system,
                                     tracker.periodic,
                                     idx, direction)
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
Directional.l2(x :: CorrelationTracker) = x.corrdata[Directional.l2]
Directional.s2(x :: CorrelationTracker) = x.corrdata[Directional.s2]

# Array interface
Base.size(x :: CorrelationTracker) = size(x.system)
Base.getindex(x :: CorrelationTracker, idx :: Vararg{Int}) = getindex(x.system, idx...)
function Base.setindex!(x   :: CorrelationTracker,
                        val,
                        idx :: Vararg{Int})
    for corrfunc in keys(x.corrdata)
        update_corrfunc!(x, val, corrfunc, idx)
    end
    x.system[idx...] = val
end
