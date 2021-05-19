directions_2d = [:x, :y, :xy_main, :xy_anti]
directions_3d = [:x, :y, :z,
                 :xy_main, :xy_anti,
                 :xz_main, :xz_anti,
                 :yz_main, :yz_anti,
                 :diag1, :diag2, :diag3, :diag4]

function test_tracker!(array :: AbstractArray, directions)
    for p in (false, true)
        for phase in (0, 1)
            tracker = CorrelationTracker(array, phase;
                                         periodic   = p,
                                         directions = directions)
            for n in 1:100
                idx = [rand(1:d) for d in size(array)]
                array[idx...]   = 1 - array[idx...]
                tracker[idx...] = 1 - tracker[idx...]
            end

            for func in (Directional.l2, Directional.s2)
                expected = func(array, phase;
                                periodic   = p,
                                directions = directions)
                got = func(tracker)
                for (direction, data) in expected
                    @test data ≈ got[direction]
                end
            end
        end
    end
end

@testset "2D system" begin
    array = rand(0:1, (100, 50))
    test_tracker!(array, directions_2d)
end

@testset "3D system" begin
    array = rand(0:1, (100, 50, 200))
    test_tracker!(array, directions_3d)
end
