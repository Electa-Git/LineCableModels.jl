@testitem "Validation / rule classes / ordering and predicate equivalence classes" tags = [:unit] begin
    import LineCableModels.Validation as V

    ordering_cases = (
        (V.Greater(:left, :right), (left = 2.0, right = 1.0), (left = 1.0, right = 2.0)),
        (V.LessEq(:left, :right), (left = 1.0, right = 1.0), (left = 2.0, right = 1.0)),
        (V.GreaterEq(:left, :right), (left = 1.0, right = 1.0), (left = 0.0, right = 1.0))
    )
    for (rule, valid, invalid) in ordering_cases
        @test V._apply(rule, valid, NamedTuple) === true
        @test_throws ArgumentError V._apply(rule, invalid, NamedTuple)
        @test_throws ArgumentError V._apply(
            rule,
            (left = 1.0 + 1.0im, right = 1.0),
            NamedTuple
        )
    end

    predicate = V.Satisfies(
        (:left, :right),
        (left, right) -> left + right ≤ 1.0,
        "sum must not exceed one"
    )
    @test V._apply(predicate, (left = 0.25, right = 0.5), NamedTuple) === nothing
    error = try
        V._apply(predicate, (left = 0.75, right = 0.5), NamedTuple)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    @test occursin("sum must not exceed one", sprint(showerror, error))
end

@testitem "Validation / traits / defaults, sanitization, and radius fallbacks" tags = [:unit] begin
    import LineCableModels.Validation as V

    struct PlainInput end
    struct TupleDefaults end
    struct NamedDefaults end
    struct BadDefaultCount end
    struct BadDefaultType end

    V.required_fields(::Type{TupleDefaults}) = (:value,)
    V.keyword_fields(::Type{TupleDefaults}) = (:scale, :label)
    V.keyword_defaults(::Type{TupleDefaults}) = (2.0, :default)
    V.keyword_fields(::Type{NamedDefaults}) = (:scale,)
    V.keyword_defaults(::Type{NamedDefaults}) = (scale = 3.0,)
    V.keyword_fields(::Type{BadDefaultCount}) = (:left, :right)
    V.keyword_defaults(::Type{BadDefaultCount}) = (1.0,)
    V.keyword_fields(::Type{BadDefaultType}) = (:value,)
    V.keyword_defaults(::Type{BadDefaultType}) = 1.0

    @test !V.has_radii(PlainInput)
    @test !V.has_temperature(PlainInput)
    @test isempty(V.extra_rules(PlainInput))
    @test isempty(V.required_fields(PlainInput))
    @test isempty(V.keyword_fields(PlainInput))
    @test isempty(V.keyword_defaults(PlainInput))
    @test V.parse(PlainInput, (value = 1,)) == (value = 1,)
    @test V.is_radius_input(PlainInput, 1.0)
    @test !V.is_radius_input(PlainInput, 1.0im)
    @test V.is_radius_input(PlainInput, Val(:other), 1.0)
    @test !V.is_radius_input(PlainInput, Val(:r_in), "radius")
    @test !V.is_radius_input(PlainInput, Val(:r_ex), "radius")

    @test V._kwdefaults_nt(PlainInput) == NamedTuple()
    @test V._kwdefaults_nt(TupleDefaults) == (scale = 2.0, label = :default)
    @test V._kwdefaults_nt(NamedDefaults) == (scale = 3.0,)
    @test_throws ErrorException V._kwdefaults_nt(BadDefaultCount)
    @test_throws ErrorException V._kwdefaults_nt(BadDefaultType)

    sanitized = V.sanitize(TupleDefaults, (4.0,), (scale = 5.0,))
    @test sanitized == (scale = 5.0, label = :default, value = 4.0)
    @test_throws ArgumentError V.sanitize(TupleDefaults, (), NamedTuple())
    @test_throws ArgumentError V.sanitize(TupleDefaults, (1.0,), (unknown = 2.0,))
end
