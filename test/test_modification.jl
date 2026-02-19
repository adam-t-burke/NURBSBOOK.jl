using Test
using NURBSBOOK
using LinearAlgebra

@testset "Shape Modification (Ch 11)" begin

    @testset "move_control_point — BSplineCurve" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [0.5, 1.0], [1.0, 0.0]]
        crv = BSplineCurve(p, U, P)

        delta = [0.0, 0.5]
        moved = move_control_point(crv, 2, delta)

        @test moved.controlpoints[2] ≈ [0.5, 1.5] atol=1e-14
        @test moved.controlpoints[1] ≈ P[1] atol=1e-14
        @test moved.controlpoints[3] ≈ P[3] atol=1e-14

        # Endpoints unchanged
        @test curve_point(moved, 0.0) ≈ P[1] atol=1e-14
        @test curve_point(moved, 1.0) ≈ P[3] atol=1e-14
    end

    @testset "modify_weight — NURBSCurve" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        p = 2
        P = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]
        w = [1.0, 1.0, 1.0]
        crv = NURBSCurve(p, U, P, w)

        pt_orig = curve_point(crv, 0.5)

        crv_high = modify_weight(crv, 2, 5.0)
        pt_high = curve_point(crv_high, 0.5)

        # Higher weight at P2 should pull the midpoint toward P2
        @test pt_high[2] > pt_orig[2]

        crv_low = modify_weight(crv, 2, 0.2)
        pt_low = curve_point(crv_low, 0.5)

        # Lower weight should push away from P2
        @test pt_low[2] < pt_orig[2]
    end

    @testset "move_control_point — BSplineSurface" begin
        uknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 2, 2)
        P[1, 1] = [0.0, 0.0, 0.0]
        P[2, 1] = [1.0, 0.0, 0.0]
        P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 2] = [1.0, 1.0, 0.0]
        surf = BSplineSurface(1, 1, uknots, vknots, P)

        moved = move_control_point(surf, 2, 2, [0.0, 0.0, 1.0])
        @test moved.controlpoints[2, 2] ≈ [1.0, 1.0, 1.0] atol=1e-14

        # Corner (0,0) unchanged
        @test surface_point(moved, 0.0, 0.0) ≈ [0.0, 0.0, 0.0] atol=1e-14
    end

    @testset "modify_weight — NURBSSurface" begin
        uknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 2, 2)
        P[1, 1] = [0.0, 0.0, 0.0]
        P[2, 1] = [1.0, 0.0, 0.0]
        P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 2] = [1.0, 1.0, 0.0]
        w = ones(2, 2)
        surf = NURBSSurface(1, 1, uknots, vknots, P, w)

        modified = modify_weight(surf, 2, 2, 3.0)
        @test modified.weights[2, 2] ≈ 3.0
        @test modified.weights[1, 1] ≈ 1.0
    end
end
