using Test
using NURBS
using LinearAlgebra

@testset "Knot Operations (Ch 5)" begin

    @testset "A5.1 CurveKnotIns — shape preservation" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0])
        p = 3
        P = [
            [0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 1.5],
            [4.0, 0.0], [5.0, 1.0], [6.0, -1.0], [7.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)

        new_crv = insert_knot(crv, 2.5)

        @test length(new_crv.controlpoints) == length(P) + 1
        @test length(new_crv.knots) == length(U) + 1

        for u in [0.0, 1.0, 2.5, 4.0, 5.0]
            @test curve_point(crv, u) ≈ curve_point(new_crv, u) atol=1e-12
        end
    end

    @testset "A5.1 CurveKnotIns — inserting existing knot" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0])
        p = 3
        P = [
            [0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 1.5],
            [4.0, 0.0], [5.0, 1.0], [6.0, -1.0], [7.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)

        new_crv = insert_knot(crv, 2.0)

        for u in [0.0, 1.5, 2.0, 3.5, 5.0]
            @test curve_point(crv, u) ≈ curve_point(new_crv, u) atol=1e-12
        end
    end

    @testset "A5.1 NURBS CurveKnotIns" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]
        w = [1.0, 2.0, 1.0]
        crv = NURBSCurve(p, U, P, w)

        new_crv = insert_knot(crv, 0.5)

        for u in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test curve_point(crv, u) ≈ curve_point(new_crv, u) atol=1e-12
        end
    end

    @testset "A5.4 RefineKnotVectCurve" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0])
        p = 3
        P = [
            [0.0, 0.0], [1.0, 1.0], [2.0, 0.5], [3.0, 1.5],
            [4.0, 0.0], [5.0, 1.0], [6.0, -1.0], [7.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)

        new_crv = refine_knots(crv, [0.5, 1.5, 2.5, 3.5, 4.5])

        for u in [0.0, 0.5, 1.0, 2.5, 4.0, 5.0]
            @test curve_point(crv, u) ≈ curve_point(new_crv, u) atol=1e-10
        end
    end

    @testset "Surface knot insertion" begin
        uknots = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 3, 2)
        P[1, 1] = [0.0, 0.0, 0.0]; P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 1] = [0.5, 0.0, 1.0]; P[2, 2] = [0.5, 1.0, 1.0]
        P[3, 1] = [1.0, 0.0, 0.0]; P[3, 2] = [1.0, 1.0, 0.0]
        surf = BSplineSurface(2, 1, uknots, vknots, P)

        new_surf = insert_knot(surf, :u, 0.5)

        for u in [0.0, 0.5, 1.0], v in [0.0, 0.5, 1.0]
            @test surface_point(surf, u, v) ≈ surface_point(new_surf, u, v) atol=1e-12
        end
    end
end
