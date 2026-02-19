using Test
using NURBSBOOK

@testset "B-spline Surfaces (Ch 3)" begin

    @testset "A3.5 SurfacePoint — bilinear" begin
        uknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 2, 2)
        P[1, 1] = [0.0, 0.0, 0.0]
        P[2, 1] = [1.0, 0.0, 0.0]
        P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 2] = [1.0, 1.0, 0.0]
        surf = BSplineSurface(1, 1, uknots, vknots, P)

        @test surface_point(surf, 0.0, 0.0) ≈ [0.0, 0.0, 0.0] atol=1e-14
        @test surface_point(surf, 1.0, 0.0) ≈ [1.0, 0.0, 0.0] atol=1e-14
        @test surface_point(surf, 0.0, 1.0) ≈ [0.0, 1.0, 0.0] atol=1e-14
        @test surface_point(surf, 1.0, 1.0) ≈ [1.0, 1.0, 0.0] atol=1e-14
        @test surface_point(surf, 0.5, 0.5) ≈ [0.5, 0.5, 0.0] atol=1e-14
    end

    @testset "A3.6 SurfaceDerivatives — bilinear" begin
        uknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 2, 2)
        P[1, 1] = [0.0, 0.0, 0.0]
        P[2, 1] = [1.0, 0.0, 0.0]
        P[1, 2] = [0.0, 1.0, 0.0]
        P[2, 2] = [1.0, 1.0, 0.0]
        surf = BSplineSurface(1, 1, uknots, vknots, P)

        SKL = surface_derivatives(surf, 0.5, 0.5, 1)
        @test SKL[1, 1] ≈ [0.5, 0.5, 0.0] atol=1e-14  # S(0.5, 0.5)
        @test SKL[2, 1] ≈ [1.0, 0.0, 0.0] atol=1e-14  # dS/du = [1,0,0]
        @test SKL[1, 2] ≈ [0.0, 1.0, 0.0] atol=1e-14  # dS/dv = [0,1,0]
    end

    @testset "A3.5 SurfacePoint — corner interpolation" begin
        uknots = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        vknots = KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        P = Matrix{Vector{Float64}}(undef, 3, 3)
        for i in 1:3, j in 1:3
            P[i, j] = [Float64(i-1)/2, Float64(j-1)/2, sin(Float64(i)*Float64(j))]
        end
        surf = BSplineSurface(2, 2, uknots, vknots, P)

        @test surface_point(surf, 0.0, 0.0) ≈ P[1, 1] atol=1e-14
        @test surface_point(surf, 1.0, 0.0) ≈ P[3, 1] atol=1e-14
        @test surface_point(surf, 0.0, 1.0) ≈ P[1, 3] atol=1e-14
        @test surface_point(surf, 1.0, 1.0) ≈ P[3, 3] atol=1e-14
    end
end
