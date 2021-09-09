@testset ExtendedTestSet "checkerboard" begin
    @testset ExtendedTestSet "checkerboard" begin
        answer = [
            0.0 0.0 0.0 0.0 0.0
            0.0 0.3454915028125263 0.5590169943749475 0.5590169943749475 0.3454915028125263
            0.0 0.5590169943749475 0.9045084971874736 0.9045084971874736 0.5590169943749475
            0.0 0.5590169943749475 0.9045084971874736 0.9045084971874736 0.5590169943749475
            0.0 0.3454915028125263 0.5590169943749475 0.5590169943749475 0.3454915028125263
        ]
        test = checkerboard((5, 5))
        @test test ≈ answer
    end
end

@testset ExtendedTestSet "disk" begin
    @testset ExtendedTestSet "disk" begin
        answer = [
            0.2928932188134524 0.5 0.2928932188134524 -0.1180339887498949 -0.5811388300841898
            0.5 1.0 0.5 0.0 -0.5
            0.2928932188134524 0.5 0.2928932188134524 -0.1180339887498949 -0.5811388300841898
            -0.1180339887498949 0.0 -0.1180339887498949 -0.41421356237309515 -0.8027756377319946
            -0.5811388300841898 -0.5 -0.5811388300841898 -0.8027756377319946 -1.1213203435596424
        ]
        test = disk((5, 5); factor=1)
        @test test ≈ answer
    end

    @testset ExtendedTestSet "disk" begin
        answer = [
            -0.41421356237309515 0.0 -0.41421356237309515 -1.2360679774997898 -2.1622776601683795
            0.0 1.0 0.0 -1.0 -2.0
            -0.41421356237309515 0.0 -0.41421356237309515 -1.2360679774997898 -2.1622776601683795
            -1.2360679774997898 -1.0 -1.2360679774997898 -1.8284271247461903 -2.605551275463989
            -2.1622776601683795 -2.0 -2.1622776601683795 -2.605551275463989 -3.2426406871192848
        ]
        test = disk((5, 5); factor=2)
        @test test ≈ answer
    end
end
