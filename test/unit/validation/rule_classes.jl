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
    @test V.validate(valid) === valid

    failures=(
        RuleInput(Inf, 0.0, 1.0, 2, 1.0, 2.0, :first, "value"),
        RuleInput(1.0, -1.0, 1.0, 2, 1.0, 2.0, :first, "value"),
        RuleInput(1.0, 0.0, 0.0, 2, 1.0, 2.0, :first, "value"),
        RuleInput(1.0, 0.0, 1.0, 2, 3.0, 2.0, :first, "value")
    )
    for value in failures
        @test_throws DomainError V.validate(value)
    end
    @test_throws ArgumentError V.validate(
        RuleInput(1.0, 0.0, 1.0, 2, 1.0, 2.0, :unknown, "value")
    )
    @test_throws ArgumentError V.validate(
        RuleInput(1.0, 0.0, 1.0, 2, 1.0, 2.0, :first, 7)
    )

    mutable struct OrderedInput
        order::Vector{Symbol}
    end
    first_rule=value->push!(value.order, :first)
    second_rule=value->push!(value.order, :second)
    V.rules(::Type{OrderedInput}) = (
        V.OwnerRule(:first, first_rule),
        V.OwnerRule(:second, second_rule)
    )
    ordered=OrderedInput(Symbol[])
    @test V.validate(ordered) === ordered
    @test ordered.order == [:first, :second]
end
