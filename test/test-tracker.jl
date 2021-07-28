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

const trackers_all = [S2Tracker(0), L2Tracker(0), L2Tracker(1)]
const trackers_axial = [SSTracker(0), SVTracker(0)]

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
end

function test_rollback!(array    :: AbstractArray{T},
                        trackers :: Vector{<:AbstractTracker{T}},
                        periodic :: Bool,
                        directions) where T
    indices = CartesianIndices(array)
    tracker1 = create_tracker(array, trackers, periodic, directions)
    tracker2 = create_tracker(array, trackers, periodic, directions)

    for n in 1:100
        idx = rand(indices)
        token = update_corrfns!(tracker1, 1 - tracker1[idx], idx)
        rollback!(tracker1, token)
    end

    for data in tracked_data(tracker1)
        expected = data(tracker2)
        got      = data(tracker1)
        for (direction, data) in expected
            @test data ≈ got[direction]
        end
        @test tracker1 == tracker2
    end
end

@testset "2D system" begin
    test_tracker!(rand(0:1, (100, 50)),  trackers_all,   false, directions_2d)
    test_tracker!(rand(0:1, (100, 100)), trackers_all,   true,  directions_2d)
    test_tracker!(rand(0:1, (100, 50)),  trackers_axial, false, [:x, :y])
    test_tracker!(rand(0:1, (100, 100)), trackers_axial, true,  [:x, :y])
end

@testset "3D system" begin
    test_tracker!(rand(0:1, (100, 50,  200)), trackers_all,   false, directions_3d)
    test_tracker!(rand(0:1, (100, 100, 100)), trackers_all,   true,  directions_3d)
    test_tracker!(rand(0:1, (100, 50,  200)), trackers_axial, false, [:x, :y, :z])
    test_tracker!(rand(0:1, (100, 100, 100)), trackers_axial, true,  [:x, :y, :z])
end

@testset "2D system (rollback)" begin
    test_rollback!(rand(0:1, (100, 50)),  trackers_all,   false, directions_2d)
    test_rollback!(rand(0:1, (100, 100)), trackers_all,   true,  directions_2d)
    test_rollback!(rand(0:1, (100, 50)),  trackers_axial, false, [:x, :y])
    test_rollback!(rand(0:1, (100, 100)), trackers_axial, true,  [:x, :y])
end

@testset "3D system (rollback)" begin
    test_rollback!(rand(0:1, (100, 50,  200)), trackers_all,   false, directions_3d)
    test_rollback!(rand(0:1, (100, 100, 100)), trackers_all,   true,  directions_3d)
    test_rollback!(rand(0:1, (100, 50,  200)), trackers_axial, false, [:x, :y, :z])
    test_rollback!(rand(0:1, (100, 100, 100)), trackers_axial, true,  [:x, :y, :z])
end
