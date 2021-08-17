# Extrapolate correlation functions on three dimensions

## KLUDGE: We only store correlation functions for 3D image here, but
## still have to be subtype of AbstractArray to be able to tell
## dimensionality and element type of an image.

struct ExtrapolatedData{T, N} <: AbstractArray{T, N}
    periodic   :: Bool
    corrdata   :: Dict{AbstractTracker{T}, Directional.CorrelationData}
    corrlen    :: Int
    directions :: Vector{Symbol}
    shape      :: NTuple{N, Int}
end

# AnnealingAPI.jl Interface
AnnealingAPI.tracked_data(x :: ExtrapolatedData) = x.corrdata |> keys
AnnealingAPI.tracked_length(x :: ExtrapolatedData) = x.corrlen
AnnealingAPI.tracked_directions(x :: ExtrapolatedData) = x.directions

# Array interface
Base.size(x :: ExtrapolatedData) = x.shape
Base.getindex(x :: ExtrapolatedData{T}, idx :: Vararg{Int}) where T = zero(T)

function new_shape(shape :: NTuple{N, Int}, dimensions :: Int) where N
    mindim = min(N, dimensions)
    minelt = minimum(shape)
    return (shape[1:mindim]..., Tuple(minelt for n in 1:(dimensions - N))...)
end

function extrapolate_vector(data           :: Vector{Float64},
                            data_grid_step :: Float64,
                            extr_grid_step :: Float64)
    grid = range(0; length = length(data), step = data_grid_step)
    itp  = interpolate((grid,), data, Gridded(Linear()))
    ext  = extrapolate(itp, Line())

    new_grid = range(0; length = length(data), step = extr_grid_step)
    return [ext(x) for x in new_grid]
end

function extrapolate_data(data       :: Directional.CorrelationData,
                          directions :: Vector{Symbol})
    # TODO: use real extrapolation here
    corrlen = length(data)
    unit_length = Directional.unit_length

    local function statistics_gen(dict)
        Iterators.map(directions) do direction
            if direction ∈ data
                direction => dict[direction]
            else
                exdata = (extrapolate_vector(dict[dir], unit_length(dir), unit_length(direction))
                          for dir in Directional.directions(data))
                direction => reduce(.+, exdata)
            end
        end
    end

    total   = statistics_gen(data.total)
    success = statistics_gen(data.success)

    return Directional.CorrelationData(directions, Dict(success), Dict(total))
end

function ExtrapolatedData(tracker    :: CorrelationTracker{T, N},
                          dimensions :: Int,
                          directions :: Vector{Symbol}) where {T, N}
    shape = new_shape(size(tracker), dimensions)
    corrdata = CorrdataDict{T}(tracker => extrapolate_data(data, directions)
                               for (tracker, data) in tracker.corrdata)
    return ExtrapolatedData{T, dimensions}(tracker.periodic,
                                           corrdata,
                                           tracker.corrlen,
                                           directions,
                                           shape)
end

# TODO: Move to CorrelationFunctions.jl maybe?
function join_data(data1 :: Directional.CorrelationData,
                   data2 :: Directional.CorrelationData)
    @assert Directional.directions(data1) == Directional.directions(data2)
    directions = Directional.directions(data1)

    success = Dict(dir => data1.success[dir] + data2.success[dir]
                   for dir in directions)
    total   = Dict(dir => data1.total[dir] + data2.total[dir]
                   for dir in directions)
    return Directional.CorrelationData(directions, success, total)
end

function ExtrapolatedData(data1 :: ExtrapolatedData{T, N},
                          data2 :: ExtrapolatedData{T, N}) where {T, N}
    if data1.periodic             != data2.periodic        ||
        size(data1)               != size(data2)           ||
        tracked_data(data1)       != tracked_data(data2)   ||
        tracked_length(data1)     != tracked_length(data2) ||
        tracked_directions(data1) != tracked_directions(data2)
        error("Two ExtrapolationData objects must share the same properties to be joined")
    end

    joined_data = Dict(tracker => join_data(tracker(data1), tracker(data2))
                       for tracker in tracked_data(data1))

    return ExtrapolatedData(data1.periodic,
                            joined_data,
                            tracked_length(data1),
                            tracked_directions(data1),
                            size(data1))
end


Directional.l2(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[L2Tracker{T}(phase)]

Directional.s2(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[S2Tracker{T}(phase)]

Directional.surfsurf(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[SSTracker{T}(phase)]

Directional.surfvoid(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[SVTracker{T}(phase)]
