@testitem "Validation / rules / native failures and nonmutating evaluation" tags=[:unit] begin
    import LineCableModels.Validation as V

    struct RuleInput
        finite::Float64
        nonnegative::Float64
        positive::Float64
        count::Int
        left::Float64
        right::Float64
        mode::Symbol
        payload::Any
    end
    V.rules(::Type{RuleInput}) = (
        V.Finite(:finite),
        V.Nonnegative(:nonnegative),
        V.Positive(:positive),
        V.IntegerField(:count),
        V.Less(:left, :right),
        V.LessEqual(:left, :right),
        V.Greater(:right, :left),
        V.GreaterEqual(:right, :left),
        V.OneOf(:mode, (:first, :second)),
        V.IsA{AbstractString}(:payload)
    )

    valid=RuleInput(1.0, 0.0, 1.0, 2, 1.0, 2.0, :first, "value")
    @test V.check(RuleInput, valid) === valid

    failures=(
        RuleInput(Inf, 0.0, 1.0, 2, 1.0, 2.0, :first, "value"),
        RuleInput(1.0, -1.0, 1.0, 2, 1.0, 2.0, :first, "value"),
        RuleInput(1.0, 0.0, 0.0, 2, 1.0, 2.0, :first, "value"),
        RuleInput(1.0, 0.0, 1.0, 2, 3.0, 2.0, :first, "value")
    )
    for value in failures
        @test_throws DomainError V.check(RuleInput, value)
    end
    @test_throws ArgumentError V.check(
        RuleInput,
        RuleInput(1.0, 0.0, 1.0, 2, 1.0, 2.0, :unknown, "value")
    )
    @test_throws ArgumentError V.check(
        RuleInput,
        RuleInput(1.0, 0.0, 1.0, 2, 1.0, 2.0, :first, 7)
    )

    source=(left = 1.0, right = 2.0)
    @test V.apply(V.Less(:left, :right), source, NamedTuple) === true
    @test source == (left = 1.0, right = 2.0)
    @test_throws ArgumentError V.apply(
        V.Less(:left, :right),
        (left = 1.0 + 1.0im, right = 2.0),
        NamedTuple
    )
end
