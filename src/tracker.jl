function gradient(array :: AbstractArray)
    deltas = imgradients(array, KernelFactors.sobel)
    return map((x...) -> norm(x), deltas...)
end

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

struct CorrelationTracker{T, N, A} <: AbstractArray{T, N}
    system     :: A
    periodic   :: Bool
    corrdata   :: Dict{TrackedData{T}, Directional.CorrelationData}
    grad       :: Array{Float64, N}

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
default_trackers(T :: Type) = 
    [TrackedData{T}(Directional.s2, 0),
     TrackedData{T}(Directional.l2, 1),
     TrackedData{T}(Directional.l2, 0)]

"""
    CorrelationTracker(system   :: AbstractArray{T, N}; 
                       tracking = default_trackers(T),
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
                            tracking   :: Vector{TrackedData{T}} = default_trackers(T),
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
    len = length(first(corrdata)[2])
    # FIXME: What about multiphase systems?
    return CorrelationTracker{T, N, typeof(system)}(
        copy(system), periodic, corrdata, gradient(system), len, directions)
end

# Looks ugly
abstract type Updater end
struct L2Updater <: Updater end
struct SSUpdater <: Updater end

function make_updater(fn :: Function)
    if fn === Directional.s2 || fn === Directional.l2
        return L2Updater()
    elseif fn === Directional.surfsurf
        return SSUpdater()
    end
end

function update_pre!(tracker  :: CorrelationTracker{T, N},
                     data     :: TrackedData{T},
                     _        :: L2Updater,
                     val,
                     idx      :: Tuple) where {T, N}
    corrdata = tracker.corrdata[data]
    corrfunc = data.func
    phase    = data.phase
    len = length(corrdata)

    for direction in Directional.directions(corrdata)
        slice, _ = get_slice(tracker.system,
                             tracker.periodic,
                             idx, direction)
        scorr = corrfunc(slice, phase;
                         periodic = tracker.periodic, len = len)
        corrdata.success[direction] .-= scorr.success[:x]
    end

    return nothing
end

function update_pre!(tracker  :: CorrelationTracker{T, N},
                     data     :: TrackedData{T},
                     _        :: SSUpdater,
                     val,
                     index    :: Tuple) where {T, N}
    index = CartesianIndex(index)
    corrdata = tracker.corrdata[data]
    corrfunc = data.func
    grad     = tracker.grad
    phase    = data.phase
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
                                periodic = tracker.periodic, len = len)
            corrdata.success[direction] .-= s2.success[:x]
        end
    end

    return nothing
end

function update_post!(tracker  :: CorrelationTracker{T, N},
                      data     :: TrackedData{T},
                      _        :: SSUpdater,
                      val,
                      index    :: Tuple) where {T, N}
    index = CartesianIndex(index)
    corrdata = tracker.corrdata[data]
    corrfunc = data.func
    grad     = tracker.grad
    phase    = data.phase
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
                                periodic = tracker.periodic, len = len)
            corrdata.success[direction] .+= s2.success[:x]
        end
    end

    return nothing
end

function update_post!(tracker  :: CorrelationTracker{T, N},
                      data     :: TrackedData{T},
                      _        :: L2Updater,
                      val,
                      idx      :: Tuple) where {T, N}
    corrdata = tracker.corrdata[data]
    corrfunc = data.func
    phase    = data.phase
    len = length(corrdata)

    for direction in Directional.directions(corrdata)
        slice, _ = get_slice(tracker.system,
                             tracker.periodic,
                             idx, direction)
        scorr = corrfunc(slice, phase;
                         periodic = tracker.periodic, len = len)
        corrdata.success[direction] .+= scorr.success[:x]
    end

    return nothing
end

update_pre!(tracker  :: CorrelationTracker{T, N},
            data     :: TrackedData{T},
            val,
            idx      :: Tuple) where {T, N} =
                update_pre!(tracker, data, make_updater(data.func), val, idx)

update_post!(tracker  :: CorrelationTracker{T, N},
             data     :: TrackedData{T},
             val,
             idx      :: Tuple) where {T, N} =
                 update_post!(tracker, data, make_updater(data.func), val, idx)

"""
    tracked_data(x :: CorrelationTracker)

Return an iterator over `TrackedData` objects which are tracked by the
tracker.
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
Directional.l2(x :: CorrelationTracker{T, N}, phase) where {T, N} =
    TrackedData{T}(Directional.l2, phase)(x)

@doc raw"""
    Directional.s2(x :: CorrelationTracker, phase)

Return $S_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.s2(x :: CorrelationTracker{T, N}, phase) where {T, N} =
    TrackedData{T}(Directional.s2, phase)(x)

@doc raw"""
    Directional.surfsurf(x :: CorrelationTracker, phase)

Return $F_{ss}_2^{\text{phase}}$ function for an underlying system of the
tracker `x`.
"""
Directional.surfsurf(x :: CorrelationTracker{T, N}, phase) where {T, N} =
    TrackedData{T}(Directional.surfsurf, phase)(x)

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
