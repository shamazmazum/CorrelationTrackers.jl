"""
    TrackedData{T}(func :: Function, phase :: T)

Construct a pair of correlation function `func` and phase `phase`
which must be tracked in `CorrelationTracker`.

# Examples
```jldoctest
julia> TrackedData(Directional.l2, 1)
TrackedData{Int64}(CorrelationFunctions.Directional.l2, 1)
```
"""
struct TrackedData{T}
    func  :: Function
    phase :: T
end

struct CorrelationTracker{T, N} <: AbstractArray{T, N}
    system     :: Array{T, N}
    periodic   :: Bool
    corrdata   :: Dict{TrackedData{T}, Directional.CorrelationData}
    softupdate :: Bool
end

"""
A vector of correlation functions and phases which are tracked by default.
"""
const tracking_by_default =
    [TrackedData(Directional.s2, 1),
     TrackedData(Directional.l2, 1),
     TrackedData(Directional.l2, 0)]

"""
    CorrelationTracker{T, N}(system   :: AbstractArray{T, N}; 
                             tracking = tracking_by_default,
                             periodic = false[, directions][, kwargs...])

Create correlation functions tracker.

Create correlation tracker for the array `system`. `tracking` is a
vector of `TrackedData` structures which specify correlation functions
you wish to track. `periodic` and `direction` have the same meaning as
in the most functions in `CorrelationFunctions.jl` package. Additional
arguments such as `len` may be passed in `kwargs`.

Returned tracker supports interface of `AbstractArray` (e.g. you can
perform element-wise read and write operations).

# Examples
```jldoctest
julia> let
       system = rand(MersenneTwister(35), 0:1, (30, 10))
       tracker = CorrelationTracker{Int,2}(system)
       end
30×10 CorrelationTracker{Int64, 2}:
 0  1  0  1  1  0  0  1  1  0
 1  1  1  0  0  0  0  0  1  1
 0  0  0  0  0  0  1  1  0  1
 1  1  1  0  1  1  1  0  1  0
 0  1  0  0  1  0  0  1  1  1
 0  0  0  0  0  0  1  0  1  1
 0  0  1  0  1  1  0  1  0  1
 1  0  0  1  0  0  1  0  1  0
 0  1  1  0  0  1  1  1  1  1
 0  0  1  1  1  1  0  0  0  0
 0  0  1  1  0  0  1  1  1  0
 0  1  0  0  0  1  0  0  1  0
 1  0  0  1  0  0  1  1  0  1
 0  1  0  1  0  0  1  1  1  0
 1  1  0  1  1  1  0  1  0  1
 1  1  1  0  0  0  0  1  0  1
 1  0  0  1  0  0  1  1  1  0
 0  0  0  1  0  0  0  1  1  0
 1  0  1  0  1  0  0  0  1  0
 1  0  0  1  0  0  0  0  0  1
 1  1  1  0  1  0  1  0  1  1
 0  1  0  1  1  0  0  1  0  1
 0  0  0  1  0  0  1  1  1  1
 0  0  1  1  1  1  0  1  1  0
 1  0  1  1  0  0  0  0  0  1
 1  1  0  1  0  1  1  0  1  0
 0  1  1  0  0  1  1  0  1  0
 0  1  0  0  1  0  0  1  0  0
 1  1  0  0  1  1  0  0  0  1
 0  0  1  1  0  1  1  1  1  0
```
"""
function CorrelationTracker{T, N}(system     :: AbstractArray{T, N};
                                  tracking   :: Vector{TrackedData{T}} = tracking_by_default,
                                  periodic   :: Bool                   = false,
                                  directions :: Vector{Symbol}         =
                                      system |> Directional.default_directions,
                                  kwargs...) where {T, N}
    corrdata = Dict{TrackedData{T},
                    Directional.CorrelationData}(data =>
        data.func(system, data.phase;
                  periodic   = periodic,
                  directions = directions,
                  kwargs...)
                                                 for data in tracking)
    return CorrelationTracker(copy(system), periodic, corrdata, false)
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

Return an iterator over `TrackedData` objects which are tracked by the
tracker.
"""
tracked_data(x :: CorrelationTracker) = x.corrdata |> keys

# Correlation functions interface
# ! FIXME: It should be safe to return internal structures as long as
# noone is going to modify them. I see no such scenario.
# Make TrackedData callable
(data :: TrackedData{T})(tracker :: CorrelationTracker{T, N}) where {T, N} =
    tracker.corrdata[data]
(data :: TrackedData{T})(array :: AbstractArray{T, N}; kwargs...) where {T, N} =
    data.func(array, data.phase; kwargs...)

@doc raw"""
    Directional.l2(x :: CorrelationTracker, phase)

Return $L_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.l2(x :: CorrelationTracker, phase) =
    TrackedData(Directional.l2, phase)(x)

@doc raw"""
    Directional.s2(x :: CorrelationTracker, phase)

Return $S_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.s2(x :: CorrelationTracker, phase) =
    TrackedData(Directional.s2, phase)(x)

# Array interface
Base.size(x :: CorrelationTracker) = size(x.system)
Base.getindex(x :: CorrelationTracker, idx :: Vararg{Int}) = getindex(x.system, idx...)
function Base.setindex!(x   :: CorrelationTracker,
                        val,
                        idx :: Vararg{Int})
    if x.softupdate
        error("Soft updates are read-only")
    end

    for tracked_data in keys(x.corrdata)
        update_corrfunc!(x, tracked_data, val, idx)
    end
    x.system[idx...] = val
end

"""
    softupdate(x :: CorrelationTracker, val, idx...)

Perform soft update of correlation functions.

This function is equivalent to `x[idx...] = val`, but instead of
modifying x it returns another `CorrelationTracker` with updated
correlation functions. The new tracker is read-only, i.e you cannot
assign values to elements of an underlying array.

# Examples
```jldoctest
julia> let
       array = rand(0:1, (200, 200))
       tracker = CorrelationTracker{Int, 2}(array)
       su = softupdate(tracker, 1 - tracker[43, 102], 43, 102)
       tracker[43, 102] = 1 - tracker[43, 102]
       Directional.s2(tracker, 1)[:x] == Directional.s2(su, 1)[:x]
       end
true
```
"""
function softupdate(x :: CorrelationTracker{T, N}, val, idx :: Vararg{Int}) where {T, N}
    soft = CorrelationTracker{T, N}(x.system, x.periodic, deepcopy(x.corrdata), true)
    for tracked_data in keys(x.corrdata)
        update_corrfunc!(soft, tracked_data, val, idx)
    end

    return soft
end
