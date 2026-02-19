using Test
using NURBSBOOK
using LinearAlgebra

@testset "B-spline Curves (Ch 3)" begin

    @testset "A3.1 CurvePoint — endpoint interpolation" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0])
        p = 3
        P = [
            [0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 1.5],
            [4.0, 0.0], [5.0, 1.0], [6.0, -1.0], [7.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)

        @test curve_point(crv, 0.0) ≈ P[1] atol=1e-14
        @test curve_point(crv, 5.0) ≈ P[end] atol=1e-14
    end

    @testset "A3.1 CurvePoint — linear curve" begin
        U = KnotVector([0.0, 0.0, 1.0, 1.0])
        p = 1
        P = [[0.0, 0.0], [1.0, 1.0]]
        crv = BSplineCurve(p, U, P)

        @test curve_point(crv, 0.0) ≈ [0.0, 0.0] atol=1e-14
        @test curve_point(crv, 0.5) ≈ [0.5, 0.5] atol=1e-14
        @test curve_point(crv, 1.0) ≈ [1.0, 1.0] atol=1e-14
    end

    @testset "A3.1 CurvePoint — quadratic Bezier" begin
        # Quadratic Bezier with no interior knots
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]
        crv = BSplineCurve(p, U, P)

        # C(0.5) for quadratic Bezier = (1/4)P0 + (1/2)P1 + (1/4)P2
        expected = 0.25 .* P[1] .+ 0.5 .* P[2] .+ 0.25 .* P[3]
        @test curve_point(crv, 0.5) ≈ expected atol=1e-14
    end

    @testset "A3.2 CurveDerivatives" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]
        crv = BSplineCurve(p, U, P)

        CK = curve_derivatives(crv, 0.0, 2)
        @test length(CK) == 3
        @test CK[1] ≈ P[1] atol=1e-14  # C(0) = P0

        # C'(0) = p * (P1 - P0) / (u_{p+1} - u_1) for Bezier = 2 * (P1 - P0)
        @test CK[2] ≈ 2 .* (P[2] .- P[1]) atol=1e-14

        CK_end = curve_derivatives(crv, 1.0, 2)
        @test CK_end[1] ≈ P[3] atol=1e-14
        @test CK_end[2] ≈ 2 .* (P[3] .- P[2]) atol=1e-14
    end

    @testset "A3.3 CurveDerivCpts" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]
        n = length(P)

        PK = curve_deriv_control_points(n, p, U, P, 1, 1, n)
        @test length(PK) == 2
        @test length(PK[1]) == 3  # 0th derivative has 3 control points
        @test length(PK[2]) == 2  # 1st derivative has 2 control points
    end

    @testset "CurvePoint — 3D curve" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
        p = 3
        P = [
            [0.0, 0.0, 0.0], [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0], [3.0, 0.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)

        pt = curve_point(crv, 0.5)
        @test length(pt) == 3
        @test curve_point(crv, 0.0) ≈ P[1] atol=1e-14
        @test curve_point(crv, 1.0) ≈ P[end] atol=1e-14
    end
end
