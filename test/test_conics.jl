using Test
using NURBSBOOK
using LinearAlgebra

@testset "Conics & Circles (Ch 7)" begin

    @testset "A7.1 make_arc — quarter circle" begin
        O = [0.0, 0.0]
        X = [1.0, 0.0]
        Y = [0.0, 1.0]
        r = 1.0

        arc = make_arc(O, X, Y, r, 0.0, π/2)

        pt_start = curve_point(arc, arc.knots[1])
        @test pt_start ≈ [1.0, 0.0] atol=1e-10

        pt_end = curve_point(arc, arc.knots[end])
        @test pt_end ≈ [0.0, 1.0] atol=1e-10

        # Point at u=0.5 should be on the circle
        pt_mid = curve_point(arc, 0.5)
        @test norm(pt_mid .- O) ≈ r atol=1e-8
    end

    @testset "A7.1 make_arc — semicircle" begin
        O = [0.0, 0.0]
        X = [1.0, 0.0]
        Y = [0.0, 1.0]
        r = 2.0

        arc = make_arc(O, X, Y, r, 0.0, Float64(π))

        pt_start = curve_point(arc, arc.knots[1])
        @test pt_start ≈ [2.0, 0.0] atol = 1e-10

        pt_end = curve_point(arc, arc.knots[end])
        @test pt_end ≈ [-2.0, 0.0] atol = 1e-8

        for t in range(arc.knots[1], arc.knots[end]; length=20)
            pt = curve_point(arc, t)
            @test norm(pt .- O) ≈ r atol = 1e-6
        end
    end

    @testset "make_conic_arc" begin
        P0 = [1.0, 0.0]
        T0 = [0.0, 1.0]
        P2 = [0.0, 1.0]
        T2 = [-1.0, 0.0]

        # w < 1 gives elliptic arc
        arc = make_conic_arc(P0, T0, P2, T2, cos(π/4))
        @test arc.degree == 2
        @test length(arc.controlpoints) == 3

        pt_start = curve_point(arc, 0.0)
        @test pt_start ≈ P0 atol=1e-12

        pt_end = curve_point(arc, 1.0)
        @test pt_end ≈ P2 atol=1e-12
    end
end
