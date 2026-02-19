using Test
using NURBSBOOK
using LinearAlgebra

@testset "Fitting (Ch 9)" begin

    @testset "A9.1 Global curve interpolation — collinear points" begin
        Q = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]
        crv = global_curve_interpolation(Q, 3)

        for i in 1:4
            u = (i - 1) / 3.0
            pt = curve_point(crv, u)
            @test pt ≈ Q[i] atol=1e-8
        end
    end

    @testset "A9.1 Global curve interpolation — quadratic" begin
        Q = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 0.0]]
        crv = global_curve_interpolation(Q, 2; method=:chord)

        # Curve should interpolate all data points
        for (i, qi) in enumerate(Q)
            u = (i == 1) ? 0.0 : (i == length(Q)) ? 1.0 :
                sum(norm(Q[j] - Q[j-1]) for j in 2:i) / sum(norm(Q[j] - Q[j-1]) for j in 2:length(Q))
            pt = curve_point(crv, u)
            @test norm(pt .- qi) < 1e-6
        end
    end

    @testset "A9.1 Global curve interpolation — endpoints" begin
        Q = [[0.0, 0.0], [0.5, 1.0], [1.5, 0.5], [2.0, 0.0]]
        crv = global_curve_interpolation(Q, 2)

        @test curve_point(crv, 0.0) ≈ Q[1] atol=1e-10
        @test curve_point(crv, 1.0) ≈ Q[end] atol=1e-10
    end

    @testset "Global curve approximation" begin
        Q = [[Float64(i), sin(Float64(i))] for i in 0:20]
        crv = global_curve_approximation(Q, 3, 8)

        # Approximation should pass near endpoints
        @test norm(curve_point(crv, 0.0) .- Q[1]) < 1e-6
        @test norm(curve_point(crv, 1.0) .- Q[end]) < 1e-6
    end

    @testset "Local curve interpolation" begin
        Q = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 1.0], [4.0, 0.0]]
        crv = local_curve_interpolation(Q)

        @test crv.degree == 3
        @test curve_point(crv, crv.knots[1]) ≈ Q[1] atol = 1e-8
        # Interior points should be near the data
        @test length(crv.controlpoints) >= length(Q)
    end
end
