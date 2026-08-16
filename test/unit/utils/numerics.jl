@testitem "Utils / numeric coercion / materialized scalar boundary" tags = [:unit] begin
    using Measurements: measurement, Measurement, value, uncertainty
    const Utils = LineCableModels.Utils

    @test Utils.resolve_T(Float32(1), big"2") === Float64
    @test Utils.resolve_T(ComplexF32(1, 2), big"2") === ComplexF64
    @test Utils.resolve_T(measurement(1.0, 0.1), Float32(2)) === Measurement{Float64}
    @test Utils._hascomplex_type(ComplexF64)
    @test Utils._hascomplex_type(Tuple{Float64, ComplexF32})
    @test Utils._hascomplex_type(NamedTuple{(:value,), Tuple{ComplexF64}})
    @test !Utils._hascomplex_type(Tuple{Float64, Int})

    nested = (values = [Float32(1), Float32(2)], selector = :R, missing = nothing)
    coerced = Utils.coerce_to_T(nested, Float64)
    @test coerced.values == [1.0, 2.0]
    @test eltype(coerced.values) === Float64
    @test coerced.selector === :R
    @test coerced.missing === nothing
    @test Utils._coerce_elt_to_T(true, Float64) === true
    @test ismissing(Utils._coerce_elt_to_T(missing, Float64))
    @test Utils.coerce_to_T(1.0, Float64) === 1.0
    @test Utils.coerce_to_T(1.0 + 2.0im, ComplexF64) === 1.0 + 2.0im
    @test Utils.coerce_to_T(1.0 + 2.0im, Float64) === 1.0

    uncertain = Utils.coerce_to_T(measurement(2.0, 0.25), Measurement{Float64})
    @test value(uncertain) == 2.0
    @test uncertainty(uncertain) == 0.25
end

@testitem "Utils / matrix transforms / symmetry and circulant invariants" tags = [:unit] begin
    using LinearAlgebra
    const Utils = LineCableModels.Utils

    matrix = [1.0 2.0 9.0; 4.0 5.0 8.0; 3.0 7.0 6.0]
    symmetrized = Utils.symtrans(matrix)
    @test symmetrized == transpose(symmetrized)
    @test diag(symmetrized) == diag(matrix)
    @test matrix != symmetrized
    @test Utils.symtrans!(copy(matrix)) == symmetrized

    circulant = Utils.line_transpose!(copy(matrix))
    @test all(circulant[i, j] == circulant[mod1(i + 1, 3), mod1(j + 1, 3)]
    for i in 1:3, j in 1:3)
    @test sum(circulant) ≈ sum(matrix)
    @test_throws ArgumentError Utils.line_transpose!(ones(2, 3))

    block_source = reshape(collect(1.0:16.0), 4, 4)
    transformed = Utils.block_transform(block_source, [1, 2, 1, 2], x -> 2x)
    @test transformed[[1, 3], [1, 3]] == 2 .* block_source[[1, 3], [1, 3]]
    @test transformed[[2, 4], [2, 4]] == 2 .* block_source[[2, 4], [2, 4]]
    @test transformed[[1, 3], [2, 4]] == block_source[[1, 3], [2, 4]]
    @test_throws ArgumentError Utils.block_transform(ones(2, 3), [1, 1], identity)
    @test_throws ArgumentError Utils.block_transform(ones(2, 2), [1, 1], _ -> ones(1, 1))

    @test Utils._to_σ(Inf) == 0.0
    @test isinf(Utils._to_σ(0.0))
    @test Utils._to_σ(4.0) == 0.25
    @test Utils.offdiag_ratio(Diagonal([1.0, 2.0])) == 0.0
    @test Utils.isdiag_rel(Diagonal([1.0, 2.0]))
end

@testitem "Utils / runtime boundaries / headless, blocks, and conductivity limits" tags = [:unit] begin
    const U = LineCableModels.Utils

    @test U.to_certain(3.0) === 3.0
    @test U.to_nominal(3.0 + 4.0im) == 3.0 + 4.0im
    @test U.to_nominal([1.0, 2.0]) == [1.0, 2.0]
    @test U.uncertainty_value(3.0) == 0.0
    @test U.uncertainty_value("not numeric") == 0.0
    @test U._nudge_float(1.0) == nextfloat(1.0)
    @test U._nudge_float(1.5) == 1.5
    @test U._nudge_float(Inf) == Inf
    @test U._coerce_args_to_T(Float32(1), 2.0) === Float64
    @test U._coerce_scalar_to_T(1, Float32) === 1.0f0
    @test U._coerce_array_to_T([1, 2], Float32) == Float32[1, 2]

    @test withenv("CI" => "true", "DISPLAY" => "available", "GKSwstype" => "") do
        U.is_headless()
    end
    if Sys.islinux()
        @test withenv("CI" => "false", "DISPLAY" => nothing, "GKSwstype" => "") do
            U.is_headless()
        end
    end
    @test withenv("CI" => "false", "DISPLAY" => "available", "GKSwstype" => "100") do
        U.is_headless()
    end
    @test !withenv("CI" => "false", "DISPLAY" => "available", "GKSwstype" => "") do
        U.is_headless()
    end
    @test withenv("CI" => "true") do
        U.display_path(joinpath("some", "directory", "file.txt")) == "file.txt"
    end
    @test withenv("CI" => "false", "DISPLAY" => "available", "GKSwstype" => "") do
        U.display_path(joinpath("some", "directory", "file.txt")) ==
        relpath(joinpath("some", "directory", "file.txt"))
    end
    @test !U.is_in_testset()

    tensor = reshape(collect(1.0:18.0), 3, 3, 2)
    transformed = U.block_transform(
        tensor,
        [1, 2, 1],
        (block, weights) -> block .* reshape(weights, :, 1),
        [2.0, 3.0, 4.0];
        slice_positions = [1]
    )
    @test transformed !== tensor
    @test transformed[[1, 3], [1, 3], 1] ==
          tensor[[1, 3], [1, 3], 1] .* [2.0, 4.0]
    @test_throws ArgumentError U.block_transform!(zeros(2, 3), [1, 2], identity)
    @test_throws ArgumentError U.block_transform!(zeros(2, 2, 1, 1), [1, 2], identity)
    @test_throws ArgumentError U.block_transform!(zeros(2, 2), [1, 1], _ -> zeros(1, 1))

    @test U._to_σ(Inf) == 0.0
    @test U._to_σ(0.0) == Inf
    @test U._to_σ(4.0) == 0.25
    @test U._bessel_diff(0.0 + 0.0im, 2.0, 4.0) ≈ log(2.0)
    @test isfinite(U._bessel_diff(1.0 + 1.0im, 2.0, 4.0))
    @test !U.issymmetric_approx(zeros(2, 3))
    @test_throws ArgumentError U.offdiag_ratio(zeros(2, 3))
end
