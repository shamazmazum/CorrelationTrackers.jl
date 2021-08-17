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

function extrapolate_vector(data           :: Vector{T},
                            data_grid_step :: Float64,
                            extr_grid_step :: Float64) where T
    grid = range(0; length = length(data), step = data_grid_step)
    itp  = interpolate((grid,), data, Gridded(Linear()))
    ext  = extrapolate(itp, Line())

    new_grid = range(0; length = length(data), step = extr_grid_step)
    new_data = [ext(x) for x in new_grid]

    return T <: Integer ? new_data .|> round .|> T : new_data
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


Directional.l2(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[L2Tracker{T}(phase)]

Directional.s2(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[S2Tracker{T}(phase)]

Directional.surfsurf(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[SSTracker{T}(phase)]

Directional.surfvoid(tracker :: ExtrapolatedData{T}, phase) where T =
    tracker.corrdata[SVTracker{T}(phase)]
