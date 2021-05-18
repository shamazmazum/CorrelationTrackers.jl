function test_tracker!(array :: AbstractArray)
    for p in (false, true)
        for phase in (0, 1)
            tracker = CorrelationTracker(array, phase; periodic = p)
            for n in 1:20
                idx = [rand(1:d) for d in size(array)]
                array[idx...] = 1 - array[idx...]
                tracker[idx...] = 1 - tracker[idx...]
            end

            for func in (Directional.l2, Directional.s2)
                expected = func(array, phase; periodic = p)
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
    test_tracker!(array)
end

@testset "3D system" begin
    array = rand(0:1, (100, 50, 200))
    test_tracker!(array)
end
