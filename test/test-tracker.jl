directions_2d = [:x, :y, :xy_main, :xy_anti]
directions_3d = [:x, :y, :z,
                 :xy_main, :xy_anti,
                 :xz_main, :xz_anti,
                 :yz_main, :yz_anti,
                 :diag1, :diag2, :diag3, :diag4]

test_data = [TrackedData(Directional.s2, 1),
             TrackedData(Directional.l2, 1),
             TrackedData(Directional.l2, 0)]

function test_tracker!(array    :: AbstractArray,
                       periodic :: Bool,
                       directions)
    n = ndims(array)
    tracker = CorrelationTracker{Int, n}(array, test_data;
                                         periodic   = periodic,
                                         directions = directions)
    for n in 1:100
        idx = [rand(1:d) for d in size(array)]
        array[idx...]   = 1 - array[idx...]
        tracker[idx...] = 1 - tracker[idx...]
    end

    for data in test_data
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
    test_tracker!(rand(0:1, (100, 50)),  false, directions_2d)
    test_tracker!(rand(0:1, (100, 100)), true,  directions_2d)
end

@testset "3D system" begin
    test_tracker!(rand(0:1, (100, 50,  200)), false, directions_3d)
    test_tracker!(rand(0:1, (100, 100, 100)), true,  directions_3d)
end
