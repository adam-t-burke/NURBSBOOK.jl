using Test
using NURBS
using LinearAlgebra

@testset "Degree Elevation & Reduction (Ch 5)" begin

    @testset "A5.9 DegreeElevateCurve — shape preservation" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]
        crv = BSplineCurve(p, U, P)

        elevated = degree_elevate(crv, 1)
        @test elevated.degree == 3

        for u in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test curve_point(crv, u) ≈ curve_point(elevated, u) atol=1e-10
        end
    end

    @testset "A5.9 DegreeElevateCurve — degree 1 to 2" begin
        U = KnotVector([0.0, 0.0, 1.0, 1.0])
        p = 1
        P = [[0.0, 0.0], [1.0, 1.0]]
        crv = BSplineCurve(p, U, P)

        elevated = degree_elevate(crv, 1)
        @test elevated.degree == 2

        for u in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test curve_point(crv, u) ≈ curve_point(elevated, u) atol=1e-10
        end
    end

    @testset "A5.9 DegreeElevateCurve — cubic with interior knots" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0])
        p = 3
        P = [
            [0.0, 0.0], [1.0, 2.0], [2.0, 1.0],
            [3.0, 2.0], [4.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)

        elevated = degree_elevate(crv, 1)
        @test elevated.degree == 4

        for u in [0.0, 0.2, 0.5, 0.8, 1.0]
            @test curve_point(crv, u) ≈ curve_point(elevated, u) atol=1e-8
        end
    end

    @testset "NURBS degree elevation" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]
        w = [1.0, 2.0, 1.0]
        crv = NURBSCurve(p, U, P, w)

        elevated = degree_elevate(crv, 1)
        @test elevated.degree == 3

        for u in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test curve_point(crv, u) ≈ curve_point(elevated, u) atol=1e-10
        end
    end
end
