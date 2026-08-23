@testitem "ParametricBuilder / manifests / deterministic structural identity" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    automatic=Grid((1.0, 2.0))
    named=Grid((1.0, 2.0); key = :shared)
    captured_scale=2.0
    scale=value->captured_scale*value
    first_payload=(
        automatic = automatic.key,
        repeated = automatic.key,
        named = named.key,
        type = Float64,
        transform = scale,
        dictionary = Dict(:b=>2, :a=>1),
        set = Set((3, 1, 2)),
        tuple = (1, 2.0+3.0im),
        array = reshape([1.0, 2.0], 1, 2),
        summary = SampleSummary([1.0, 2.0, 3.0]),
        sentinel = missing
    )
    second_payload=merge(
        first_payload,
        (dictionary = Dict(:a=>1, :b=>2), set = Set((2, 3, 1)))
    )

    first_tree=PB._manifest_tree(first_payload)
    second_tree=PB._manifest_tree(second_payload)
    @test isequal(first_tree, second_tree)
    @test first_tree.automatic == (automatic_grid = 1,)
    @test first_tree.repeated == first_tree.automatic
    @test first_tree.named == (named_grid = :shared,)
    @test first_tree.type == "Float64"
    @test occursin("var\"#", first_tree.transform.function_type)
    @test first_tree.array.size == (1, 2)
    @test first_tree.summary.fields.mean == 2.0
    @test PB._stable_bytes(first_tree) == PB._stable_bytes(second_tree)

    manifest_a=PB.CalculationManifest(
        (radius = 0.01,),
        first_tree,
        Formulation(),
        (policy = :full,),
        (verbosity = 0,)
    )
    manifest_b=PB.CalculationManifest(
        (radius = 0.01,),
        second_tree,
        Formulation(),
        (policy = :full,),
        (verbosity = 0,)
    )
    @test manifest_a.hash == manifest_b.hash
    @test length(manifest_a.hash) == 64
    @test manifest_a.solver == string(typeof(Formulation()))

    same_type=PB._append_result(Int[1], 2)
    @test same_type == [1, 2]
    widened=PB._append_result(Int[1], 2.5)
    @test widened == Real[1, 2.5]
    @test eltype(widened) === Real
    @test PB._append_result(nothing, 2.5) == [2.5]
    @test_throws ArgumentError PB._append_result(Int[1], "incompatible")

    configuration=PB.Configuration{Float64}(
        identity,
        (-1.0,),
        (:radius,),
        ()
    )
    failure=PB._configuration_failure(2, configuration, DomainError(-1.0))
    @test failure.index == 2
    @test failure.configuration == (radius = -1.0,)
    @test occursin("DomainError", failure.exception_type)
    @test occursin("-1.0", failure.message)
    @test all(PB._skippable_configuration_error,
        (ArgumentError("invalid"), AssertionError("invalid"),
            DimensionMismatch("invalid"), DomainError(0)))
    @test !PB._skippable_configuration_error(ErrorException("unexpected"))

    constant_histogram=UQ._histogram(fill(2.0, 4), nothing)
    @test length(constant_histogram.density) == 1
    @test first(constant_histogram.edges) < 2.0 < last(constant_histogram.edges)
    @test sum(constant_histogram.density .* diff(constant_histogram.edges)) ≈ 1.0
end
