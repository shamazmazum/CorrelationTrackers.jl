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
const trackers_axial = [SSTracker(0)]

function test_tracker!(array    :: AbstractArray{T},
                       trackers :: Vector{<:AbstractTracker{T}},
                       periodic :: Bool,
                       directions) where T
    indices = CartesianIndices(array)
    tracker = CorrelationTracker(array;
                                 tracking   = trackers,
                                 periodic   = periodic,
                                 directions = directions)
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
