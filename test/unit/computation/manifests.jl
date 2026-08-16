@testitem "Computation / manifests / deterministic structural identity" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const Cmp=LineCableModels.Computation
    const PB=LineCableModels.ParametricBuilder

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

    first_tree=Cmp._manifest_tree(first_payload)
    second_tree=Cmp._manifest_tree(second_payload)
    @test isequal(first_tree, second_tree)
    @test first_tree.automatic == (automatic_grid = 1,)
    @test first_tree.repeated == first_tree.automatic
    @test first_tree.named == (named_grid = :shared,)
    @test first_tree.type == "Float64"
    @test occursin("var\"#", first_tree.transform.function_type)
    @test first_tree.array.size == (1, 2)
    @test first_tree.summary.fields.mean == 2.0
    @test Cmp._stable_bytes(first_tree) == Cmp._stable_bytes(second_tree)

    manifest_a=Cmp.CalculationManifest(
        (radius = 0.01,),
        first_tree,
        Formulation(),
        (policy = :full,),
        (verbosity = 0,)
    )
    manifest_b=Cmp.CalculationManifest(
        (radius = 0.01,),
        second_tree,
        Formulation(),
        (policy = :full,),
        (verbosity = 0,)
    )
    @test manifest_a.hash == manifest_b.hash
    @test length(manifest_a.hash) == 64
    @test manifest_a.solver == string(typeof(Formulation()))

    same_type=Cmp._append_result(Int[1], 2)
    @test same_type == [1, 2]
    widened=Cmp._append_result(Int[1], 2.5)
    @test widened == Real[1, 2.5]
    @test eltype(widened) === Real
    @test Cmp._append_result(nothing, 2.5) == [2.5]
    @test_throws ArgumentError Cmp._append_result(Int[1], "incompatible")

    configuration=PB.Configuration{Float64}(
        identity,
        (-1.0,),
        (:radius,),
        ()
    )
    failure=Cmp._configuration_failure(2, configuration, DomainError(-1.0))
    @test failure.index == 2
    @test failure.configuration == (radius = -1.0,)
    @test occursin("DomainError", failure.exception_type)
    @test occursin("-1.0", failure.message)
    @test all(Cmp._skippable_configuration_error,
        (ArgumentError("invalid"), AssertionError("invalid"),
            DimensionMismatch("invalid"), DomainError(0)))
    @test !Cmp._skippable_configuration_error(ErrorException("unexpected"))

    constant_histogram=Cmp._histogram(fill(2.0, 4), nothing)
    @test length(constant_histogram.density) == 1
    @test first(constant_histogram.edges) < 2.0 < last(constant_histogram.edges)
    @test sum(constant_histogram.density .* diff(constant_histogram.edges)) ≈ 1.0
end
