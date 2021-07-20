function gradient(array :: AbstractArray)
    deltas = imgradients(array, KernelFactors.sobel)
    return map((x...) -> norm(x), deltas...)
end

struct CorrelationTracker{T, N, A} <: AbstractArray{T, N}
    system     :: A
    periodic   :: Bool
    corrdata   :: Dict{AbstractTracker{T},
                       Directional.CorrelationData}
    grad       :: Array{Float64, N}
    fft_plans  :: Directional.S2FTPlans

    # For quick access
    corrlen    :: Int
    directions :: Vector{Symbol}
end

@doc raw"""
    default_trackers(T)

Construct a vector of correlation functions which are tracked by
default (that is $S_2^1(x)$, $L_2^1(x)$ and $L_2^0(x)$). `T` is the
type of `x`.
"""
default_trackers(T :: Type) = [S2Tracker{T}(0), L2Tracker{T}(1), L2Tracker{T}(0)]

"""
    CorrelationTracker(system   :: AbstractArray{T, N}; 
                       tracking = default_trackers(T),
                       periodic = false[, directions][, kwargs...])

Create correlation functions tracker.

Create correlation tracker for the array `system`. `tracking` is a
vector of `AbstractTracker` structures which specify correlation
functions you wish to track. `periodic` and `direction` have the same
meaning as in the most functions in `CorrelationFunctions.jl`
package. Additional arguments such as `len` may be passed in
`kwargs`.

Returned tracker supports interface of `AbstractArray` (e.g. you can
perform element-wise read and write operations).

# Examples
```jldoctest
julia> let
       system = rand(MersenneTwister(35), 0:1, (30, 10))
       tracker = CorrelationTracker(system)
       end
30×10 CorrelationTracker{Int64, 2, Matrix{Int64}}:
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
function CorrelationTracker(system     :: AbstractArray{T, N};
                            tracking   :: Vector{<:AbstractTracker{T}} = default_trackers(T),
                            periodic   :: Bool           = false,
                            directions :: Vector{Symbol} = system |> Directional.default_directions,
                            kwargs...) where {T, N}
    corrdata = Dict{AbstractTracker{T}, Directional.CorrelationData}(
        data => data(system;
                     periodic   = periodic,
                     directions = directions,
                     kwargs...)
        for data in tracking)
    len = length(first(corrdata)[2])
    # FIXME: What about multiphase systems?
    return CorrelationTracker{T, N, typeof(system)}(
        copy(system), periodic,
        corrdata, gradient(system),
        Directional.S2FTPlans(system, periodic),
        len, directions)
end

const SimpleTracker{T}  = Union{L2Tracker{T}, S2Tracker{T}}

maybe_call_with_plans(slice :: AbstractArray{T},
                      data  :: S2Tracker{T};
                      plans :: Directional.S2FTPlans,
                      kwargs...) where T =
                          data(slice; plans = plans, kwargs...)
maybe_call_with_plans(slice :: AbstractArray{T},
                      data  :: AbstractTracker{T};
                      plans :: Directional.S2FTPlans,
                      kwargs...) where T =
                          data(slice; kwargs...)

function update_pre!(tracker  :: CorrelationTracker{T},
                     data     :: SimpleTracker{T},
                     val,
                     idx      :: Tuple) where T
    corrdata = tracker.corrdata[data]
    len = length(corrdata)

    for direction in Directional.directions(corrdata)
        slice, _ = get_slice(tracker.system,
                             tracker.periodic,
                             idx, direction)
        scorr = maybe_call_with_plans(slice, data;
                                      plans    = tracker.fft_plans,
                                      periodic = tracker.periodic,
                                      len = len)
        corrdata.success[direction] .-= scorr.success[:x]
    end

    return nothing
end

function update_pre!(tracker  :: CorrelationTracker{T},
                     data     :: SSTracker{T},
                     val,
                     index    :: Tuple) where T
    index = CartesianIndex(index)
    corrdata = tracker.corrdata[data]
    grad     = tracker.grad
    len      = length(corrdata)

    indices = CartesianIndices(grad)
    fidx, lidx = first(indices), last(indices)

    for direction in Directional.directions(corrdata)
        u = axial_index(grad, direction)
        for idx in max(index - u, fidx):min(index + u, lidx)
            slice, _ = get_slice(grad,
                                 tracker.periodic,
                                 Tuple(idx), direction)
            s2 = Directional.s2(slice, Directional.SeparableIndicator(identity);
                                periodic = tracker.periodic,
                                plans    = tracker.fft_plans,
                                len      = len)
            corrdata.success[direction] .-= s2.success[:x]
        end
    end

    return nothing
end

function update_post!(tracker  :: CorrelationTracker{T},
                      data     :: SimpleTracker{T},
                      val,
                      idx      :: Tuple) where T
    corrdata = tracker.corrdata[data]
    len = length(corrdata)

    for direction in Directional.directions(corrdata)
        slice, _ = get_slice(tracker.system,
                             tracker.periodic,
                             idx, direction)
        scorr = maybe_call_with_plans(slice, data;
                                      plans    = tracker.fft_plans,
                                      periodic = tracker.periodic,
                                      len = len)
        corrdata.success[direction] .+= scorr.success[:x]
    end

    return nothing
end

function update_post!(tracker  :: CorrelationTracker{T},
                      data     :: SSTracker{T},
                      val,
                      index    :: Tuple) where T
    index = CartesianIndex(index)
    corrdata = tracker.corrdata[data]
    grad     = tracker.grad
    len      = length(corrdata)

    indices = CartesianIndices(grad)
    fidx, lidx = first(indices), last(indices)
    uidx = oneunit(index)
    gradstart = max(index - 2uidx, fidx)
    gradstop  = min(index + 2uidx, lidx)
    subsys = tracker.system[gradstart:gradstop]
    subgrad = gradient(subsys)

    index_subgrad = index - gradstart + uidx
    indices_subgrad  = CartesianIndices(subgrad)
    fidx2, lidx2 = first(indices_subgrad), last(indices_subgrad)
    grad[max(index - uidx, fidx):min(index + uidx, lidx)] .=
        subgrad[max(index_subgrad - uidx, fidx2):min(index_subgrad + uidx, lidx2)]

    for direction in Directional.directions(corrdata)
        u = axial_index(grad, direction)
        for idx in max(index - u, fidx):min(index + u, lidx)
            slice, _ = get_slice(grad,
                                 tracker.periodic,
                                 Tuple(idx), direction)
            s2 = Directional.s2(slice, Directional.SeparableIndicator(identity);
                                periodic = tracker.periodic,
                                plans    = tracker.fft_plans,
                                len      = len)
            corrdata.success[direction] .+= s2.success[:x]
        end
    end

    return nothing
end

"""
    tracked_data(x :: CorrelationTracker)

Return an iterator over correlation function descriptors which are
tracked by the tracker.
"""
tracked_data(x :: CorrelationTracker) = x.corrdata |> keys

"""
    tracked_length(x :: CorrelationTracker)

Return maximal tracked correlation length

# Examples
```jldoctest
julia> tracked_length(CorrelationTracker{Int,2}(rand(0:1, (50, 100))))
25
```
"""
tracked_length(x :: CorrelationTracker) = x.corrlen

"""
    tracked_directions(x :: CorrelationTracker)

Return directions along which correlation functions are tracked.

# Examples
```jldoctest
julia> tracked_directions(CorrelationTracker{Int,2}(rand(0:1, (50, 100))))
2-element Vector{Symbol}:
 :x
 :y
```
"""
tracked_directions(x :: CorrelationTracker) = x.directions

# Array interface
Base.size(x :: CorrelationTracker) = size(x.system)
Base.getindex(x :: CorrelationTracker, idx :: Vararg{Int}) = getindex(x.system, idx...)
function Base.setindex!(x   :: CorrelationTracker,
                        val,
                        idx :: Vararg{Int})
    # Do everything we can before updating the value in the underlying array
    for tracked_data in keys(x.corrdata)
        update_pre!(x, tracked_data, val, idx)
    end

    # Change the value
    x.system[idx...] = val

    # Do everything else
    for tracked_data in keys(x.corrdata)
        update_post!(x, tracked_data, val, idx)
    end

    return x
end
