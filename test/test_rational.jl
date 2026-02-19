using Test
using NURBSBOOK
using LinearAlgebra

@testset "NURBS Curves & Surfaces (Ch 4)" begin

    @testset "A4.1 NURBS CurvePoint — unit weights reduce to B-spline" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
        p = 3
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0]]
        w = [1.0, 1.0, 1.0, 1.0]

        bcrv = BSplineCurve(p, U, P)
        ncrv = NURBSCurve(p, U, P, w)

        for u in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test curve_point(ncrv, u) ≈ curve_point(bcrv, u) atol=1e-12
        end
    end

    @testset "A4.1 NURBS CurvePoint — weighted" begin
        # Rational quadratic with non-unit weight
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]
        w = [1.0, 2.0, 1.0]  # Weight > 1 pulls curve toward P[2]

        ncrv = NURBSCurve(p, U, P, w)
        bcrv = BSplineCurve(p, U, P)

        pt_nurbs = curve_point(ncrv, 0.5)
        pt_bspline = curve_point(bcrv, 0.5)

        # NURBS with w2 > 1 should be pulled toward P2 = [1,1]
        @test pt_nurbs[2] > pt_bspline[2]

        # Endpoints must still be interpolated
        @test curve_point(ncrv, 0.0) ≈ P[1] atol=1e-12
        @test curve_point(ncrv, 1.0) ≈ P[3] atol=1e-12
    end

    @testset "A4.2 NURBS CurveDerivatives" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]
        w = [1.0, 1.0, 1.0]

        ncrv = NURBSCurve(p, U, P, w)
        bcrv = BSplineCurve(p, U, P)

        CK_n = curve_derivatives(ncrv, 0.5, 1)
        CK_b = curve_derivatives(bcrv, 0.5, 1)

        @test CK_n[1] ≈ CK_b[1] atol=1e-12
        @test CK_n[2] ≈ CK_b[2] atol=1e-12
    end

    @testset "A4.3 NURBS SurfacePoint" begin
        uknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 2, 2)
        P[1, 1] = [0.0, 0.0, 0.0]
        P[2, 1] = [1.0, 0.0, 0.0]
        P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 2] = [1.0, 1.0, 0.0]
        w = ones(2, 2)
        surf = NURBSSurface(1, 1, uknots, vknots, P, w)

        @test surface_point(surf, 0.5, 0.5) ≈ [0.5, 0.5, 0.0] atol=1e-12
        @test surface_point(surf, 0.0, 0.0) ≈ P[1, 1] atol=1e-12
        @test surface_point(surf, 1.0, 1.0) ≈ P[2, 2] atol=1e-12
    end
end
