using Test
using NURBSBOOK

@testset "Basis Functions (Ch 2)" begin

    @testset "A2.1 FindSpan" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0])
        p = 2
        n = length(U) - p - 1  # 8 control points

        @test find_span(n, p, 0.0, U) == p + 1
        @test find_span(n, p, 2.5, U) == 5
        @test find_span(n, p, 5.0, U) == n
    end

    @testset "A2.2 BasisFuns — partition of unity" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0])
        p = 3
        n = length(U) - p - 1

        for u in [0.0, 0.5, 1.0, 2.5, 4.99, 5.0]
            span = find_span(n, p, u, U)
            N = basis_functions(span, u, p, U)
            @test length(N) == p + 1
            @test sum(N) ≈ 1.0 atol = 1e-14
            @test all(x -> x >= 0, N)
        end
    end

    @testset "A2.2 BasisFuns — Ex2.3 values" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0])
        p = 2
        n = length(U) - p - 1
        u = 2.5
        span = find_span(n, p, u, U)
        @test span == 5
        N = basis_functions(span, u, p, U)
        @test N[1] ≈ 1 / 8 atol = 1e-14
        @test N[2] ≈ 6 / 8 atol = 1e-14
        @test N[3] ≈ 1 / 8 atol = 1e-14
    end

    @testset "A2.2 BasisFuns — endpoint interpolation" begin
        U = KnotVector([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0])
        p = 2
        n = length(U) - p - 1

        N0 = basis_functions(find_span(n, p, 0.0, U), 0.0, p, U)
        @test N0[1] ≈ 1.0 atol = 1e-14  # first basis fn = 1 at left end

        N1 = basis_functions(find_span(n, p, 1.0, U), 1.0, p, U)
        @test N1[end] ≈ 1.0 atol = 1e-14  # last basis fn = 1 at right end
    end

    @testset "A2.3 DersBasisFuns" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0])
        p = 2
        n = length(U) - p - 1
        span = find_span(n, p, 2.5, U)
        ders = basis_function_derivatives(span, 2.5, p, 2, U)

        # Function values
        @test ders[1, 1] ≈ 1 / 8 atol = 1e-14
        @test ders[1, 2] ≈ 6 / 8 atol = 1e-14
        @test ders[1, 3] ≈ 1 / 8 atol = 1e-14
        # 1st derivatives: [-0.5, 0.0, 0.5]
        @test ders[2, 1] ≈ -0.5 atol = 1e-14
        @test ders[2, 2] ≈ 0.0 atol = 1e-14
        @test ders[2, 3] ≈ 0.5 atol = 1e-14
        @test abs(sum(ders[2, :])) < 1e-12
        # 2nd derivatives: [1.0, -2.0, 1.0]
        @test ders[3, 1] ≈ 1.0 atol = 1e-14
        @test ders[3, 2] ≈ -2.0 atol = 1e-14
        @test ders[3, 3] ≈ 1.0 atol = 1e-14
        @test abs(sum(ders[3, :])) < 1e-12
    end

    @testset "A2.4 OneBasisFun" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0])
        p = 2
        n = length(U) - p - 1
        span = find_span(n, p, 2.5, U)
        N = basis_functions(span, 2.5, p, U)

        for j in 1:(p + 1)
            fn_idx = span - p + j - 1
            Nj = one_basis_function(p, U, fn_idx, 2.5)
            @test Nj ≈ N[j] atol = 1e-14
        end
        @test one_basis_function(p, U, 1, 0.0) ≈ 1.0 atol = 1e-14
        @test one_basis_function(p, U, n, 5.0) ≈ 1.0 atol = 1e-14
    end

    @testset "A2.5 DersOneBasisFun" begin
        U = KnotVector([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0])
        p = 2
        n = length(U) - p - 1
        span = find_span(n, p, 2.5, U)
        ders_all = basis_function_derivatives(span, 2.5, p, 2, U)

        for j in 1:(p + 1)
            fn_idx = span - p + j - 1
            ders_one = one_basis_function_derivatives(p, U, fn_idx, 2.5, 2)
            for k in 1:3
                @test ders_one[k] ≈ ders_all[k, j] atol = 1e-12
            end
        end
    end
end
