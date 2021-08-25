const directions_2d = [
    :x, :y, :xy, :yx
]

const directions_3d = [
    :x, :y, :z,
    :xy, :yx,
    :xz, :zx,
    :yz, :zy,
    :xyz, :yxz, :xzy, :zyx
]

const trackers = [
    S2Tracker(false),
    L2Tracker(false),
    L2Tracker(true),
    SSTracker(false),
    SVTracker(false),
    CCTracker(false, true),
    CCTracker(true, false)
]

create_tracker(array, trackers, periodic, directions) =
    CorrelationTracker(array;
                       tracking   = trackers,
                       periodic   = periodic,
                       directions = directions)

function test_tracker!(array    :: AbstractArray{T},
                       trackers :: Vector{<:AbstractTracker{T}},
                       periodic :: Bool,
                       directions) where T
    indices = CartesianIndices(array)
    tracker = create_tracker(array, trackers, periodic, directions)
    for n in 1:100
        idx = rand(indices)
        array[idx]   = 1 - array[idx]
        tracker[idx] = 1 - tracker[idx]
    end

    for data in tracked_data(tracker)
        expected = data(array;
                        periodic   = periodic,
                        directions = directions)
        got = data(tracker)
        for (direction, data) in expected
            @test data ≈ got[direction]
        end
    end
    @test tracker == array
end

function test_rollback!(array    :: AbstractArray{T},
                        trackers :: Vector{<:AbstractTracker{T}},
                        periodic :: Bool,
                        directions) where T
    indices = CartesianIndices(array)
    tracker = create_tracker(array, trackers, periodic, directions)

    for n in 1:100
        idx = rand(indices)
        token = update_corrfns!(tracker, 1 - tracker[idx], idx)
        if rand() < 0.5
            rollback!(tracker, token)
        else
            array[idx] = 1 - array[idx]
        end
    end

    for data in tracked_data(tracker)
        expected = data(array;
                        periodic   = periodic,
                        directions = directions)
        got = data(tracker)
        for (direction, data) in expected
            @test data ≈ got[direction]
        end
    end
    @test tracker == array
end

@testset "2D system" begin
    test_tracker!(rand(Bool, (100, 50)),  trackers, false, directions_2d)
    test_tracker!(rand(Bool, (100, 100)), trackers, true,  directions_2d)
end

@testset "3D system" begin
    test_tracker!(rand(Bool, (100, 50,  200)), trackers, false, directions_3d)
    test_tracker!(rand(Bool, (100, 100, 100)), trackers, true,  directions_3d)
end

@testset "2D system (rollback)" begin
    test_rollback!(rand(Bool, (100, 50)),  trackers, false, directions_2d)
    test_rollback!(rand(Bool, (100, 100)), trackers, true,  directions_2d)
end

@testset "3D system (rollback)" begin
    test_rollback!(rand(Bool, (100, 50,  200)), trackers, false, directions_3d)
    test_rollback!(rand(Bool, (100, 100, 100)), trackers, true,  directions_3d)
end

@testset "Extrapolations" begin
    # TODO: test diagonals
    array = rand(Bool, (200, 200))
    tracker = CorrelationTracker(array; periodic = true)
    extra = ExtrapolatedData(tracker, 3, [:x, :y, :z])

    s21 = Directional.s2(tracker, false)
    s22 = Directional.s2(extra, false)

    @test s22[:x] ≈ s21[:x]
    @test s22[:y] ≈ s21[:y]
    @test s22[:z] ≈ (s21[:x] + s21[:y]) / 2
end
