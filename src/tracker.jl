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

CDUpdateInfo{T} = Dict{Tuple{AbstractTracker{T}, Symbol}, Vector}
struct RollbackToken{T, N}
    val      :: T
    index    :: CartesianIndex{N}
    corrdata :: CDUpdateInfo{T}
    grad     :: Array{Float64, N}
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

function update_pre!(tracker :: CorrelationTracker{T, N},
                     uinfo   :: CDUpdateInfo{T},
                     data    :: SimpleTracker{T},
                     val,
                     index   :: CartesianIndex{N}) where {T, N}
    corrdata = tracker.corrdata[data]
    len      = length(corrdata)

    for direction in Directional.directions(corrdata)
        slice = get_slice(tracker.system,
                          tracker.periodic,
                          Tuple(index), direction)
        scorr = maybe_call_with_plans(slice, data;
                                      plans    = tracker.fft_plans,
                                      periodic = tracker.periodic,
                                      len = len)
        uinfo[(data, direction)] = -scorr.success[:x]
    end

    return nothing
end

function update_pre!(tracker :: CorrelationTracker{T, N},
                     uinfo   :: CDUpdateInfo{T},
                     data    :: SSTracker{T},
                     val,
                     index   :: CartesianIndex{N}) where {T, N}
    corrdata = tracker.corrdata[data]
    grad     = tracker.grad
    len      = length(corrdata)

    indices = CartesianIndices(grad)
    fidx, lidx = first(indices), last(indices)

    for direction in Directional.directions(corrdata)
        u = axial_index(grad, direction)
        success = mapreduce(.+, max(index - u, fidx):min(index + u, lidx)) do idx
            slice = get_slice(grad,
                              tracker.periodic,
                              Tuple(idx), direction)
            s2 = Directional.s2(slice, Directional.SeparableIndicator(identity);
                                periodic = tracker.periodic,
                                plans    = tracker.fft_plans,
                                len      = len)
            -s2.success[:x]
        end

        uinfo[(data, direction)] = success
    end

    return nothing
end

function update_post!(tracker :: CorrelationTracker{T, N},
                      uinfo   :: CDUpdateInfo{T},
                      data    :: SimpleTracker{T},
                      val,
                      index   :: CartesianIndex{N}) where {T, N}
    corrdata = tracker.corrdata[data]
    len      = length(corrdata)

    for direction in Directional.directions(corrdata)
        slice = get_slice(tracker.system,
                          tracker.periodic,
                          Tuple(index), direction)
        scorr = maybe_call_with_plans(slice, data;
                                      plans    = tracker.fft_plans,
                                      periodic = tracker.periodic,
                                      len = len)

        uinfo[(data, direction)] .+= scorr.success[:x]
        corrdata.success[direction] .+= uinfo[(data, direction)]
    end

    return nothing
end

function update_post!(tracker :: CorrelationTracker{T, N},
                      uinfo   :: CDUpdateInfo{T},
                      data    :: SSTracker{T},
                      val,
                      index   :: CartesianIndex{N}) where {T, N}
    corrdata = tracker.corrdata[data]
    grad     = tracker.grad
    len      = length(corrdata)

    indices = CartesianIndices(grad)
    fidx, lidx = first(indices), last(indices)

    for direction in Directional.directions(corrdata)
        u = axial_index(grad, direction)
        success = mapreduce(.+, max(index - u, fidx):min(index + u, lidx)) do idx
            slice = get_slice(grad,
                              tracker.periodic,
                              Tuple(idx), direction)
            s2 = Directional.s2(slice, Directional.SeparableIndicator(identity);
                                periodic = tracker.periodic,
                                plans    = tracker.fft_plans,
                                len      = len)
            s2.success[:x]
        end

        uinfo[(data, direction)] .+= success
        corrdata.success[direction] .+= uinfo[(data, direction)]
    end

    return nothing
end

function update_gradient!(tracker  :: CorrelationTracker{T, N},
                          index    :: CartesianIndex{N}) where {T, N}
    grad   = tracker.grad
    system = tracker.system

    indices    = CartesianIndices(grad)
    fidx, lidx = first(indices), last(indices)
    uidx       = oneunit(index)

    # Select 5x5 square (or 5x5x5 cube) used to update gradient
    gradstart = max(index - 2uidx, fidx)
    gradstop  = min(index + 2uidx, lidx)
    subsys    = system[gradstart:gradstop]
    subgrad   = gradient(subsys)

    # Index of the updated element in subgrad
    sindex = index - gradstart + uidx
    sindices  = CartesianIndices(subgrad)
    sfidx, slidx = first(sindices), last(sindices)

    # Update gradient
    to_update = subgrad[max(sindex - uidx, sfidx):min(sindex + uidx, slidx)]
    grad[max(index - uidx, fidx):min(index + uidx, lidx)] .= to_update

    return to_update
end

function rollback_gradient!(tracker :: CorrelationTracker{T, N},
                            index   :: CartesianIndex{N},
                            subgrad :: Array{Float64, N}) where {T, N}
    grad       = tracker.grad
    indices    = CartesianIndices(grad)
    fidx, lidx = first(indices), last(indices)
    uidx       = oneunit(index)

    grad[max(index - uidx, fidx):min(index + uidx, lidx)] .= subgrad

    return nothing
end

function update_corrfns!(tracker :: CorrelationTracker{T,N},
                         val,
                         index   :: CartesianIndex{N}) where {T, N}
    trackers = keys(tracker.corrdata)
    update_info = CDUpdateInfo{T}()

    # Do everything we can before updating the value in the underlying array
    for tr in trackers
        update_pre!(tracker, update_info,
                    tr, val, index)
    end

    # Change the value
    tracker.system[index] = val

    # Update gradient if needed
    grad = Array{Float64, N}(undef, (0 for n in 1:N)...)
    if any(update_gradient_p, trackers)
        grad = update_gradient!(tracker, index)
    end

    # Do everything else
    for tr in trackers
        update_post!(tracker, update_info,
                     tr, val, index)
    end

    return RollbackToken{T, N}(val, index, update_info, grad)
end

function rollback!(tracker  :: CorrelationTracker{T, N},
                   rollback :: RollbackToken{T, N}) where {T, N}
    trackers = keys(tracker.corrdata)
    index    = rollback.index
    val      = rollback.val

    # Roll back correlation functions
    for tr in trackers
        corrdata = tracker.corrdata[tr]

        for direction in Directional.directions(corrdata)
            corrdata.success[direction] .-= rollback.corrdata[(tr, direction)]
        end
    end

    # Roll back the value of the changed element
    tracker.system[index] = val

    # Roll back gradient
    if any(update_gradient_p, trackers)
        rollback_gradient!(tracker, index, rollback.grad)
    end

    return tracker
end

# Interface
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
Base.setindex!(x :: CorrelationTracker{T}, val, idx :: Vararg{Int}) where T =
    update_corrfns!(x, val, CartesianIndex(idx))
