using Test
using NURBS
using LinearAlgebra

@testset "Advanced Algorithms (Ch 6)" begin

    @testset "Point inversion on curve" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
        p = 3
        P = [[0.0, 0.0], [1.0, 2.0], [3.0, 2.0], [4.0, 0.0]]
        crv = BSplineCurve(p, U, P)

        u_test = 0.4
        pt = curve_point(crv, u_test)
        u_found = point_inversion_curve(crv, pt; u0=0.5)
        @test abs(u_found - u_test) < 1e-6
    end

    @testset "Curve reversal" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0])
        p = 3
        P = [
            [0.0, 0.0], [1.0, 2.0], [2.0, 1.0],
            [3.0, 3.0], [4.0, 1.0], [5.0, 0.0],
        ]
        crv = BSplineCurve(p, U, P)
        rev = curve_reverse(crv)

        # C_rev(a + b - u) = C(u) where [a,b] = [0,3]
        for u in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
            @test curve_point(crv, u) ≈ curve_point(rev, 3.0 - u) atol=1e-12
        end
    end

    @testset "NURBS curve reversal" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]
        w = [1.0, 2.0, 1.0]
        crv = NURBSCurve(p, U, P, w)
        rev = curve_reverse(crv)

        for u in [0.0, 0.25, 0.5, 0.75, 1.0]
            @test curve_point(crv, u) ≈ curve_point(rev, 1.0 - u) atol=1e-12
        end
    end

    @testset "Surface reversal" begin
        uknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 2, 2)
        P[1, 1] = [0.0, 0.0, 0.0]
        P[2, 1] = [1.0, 0.0, 0.0]
        P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 2] = [1.0, 1.0, 1.0]
        surf = BSplineSurface(1, 1, uknots, vknots, P)

        rev_u = surface_reverse(surf, :u)
        for u in [0.0, 0.5, 1.0], v in [0.0, 0.5, 1.0]
            @test surface_point(surf, u, v) ≈ surface_point(rev_u, 1.0 - u, v) atol=1e-12
        end

        rev_v = surface_reverse(surf, :v)
        for u in [0.0, 0.5, 1.0], v in [0.0, 0.5, 1.0]
            @test surface_point(surf, u, v) ≈ surface_point(rev_v, u, 1.0 - v) atol=1e-12
        end
    end
end
