@testitem "Engine / option grammar / owner dispatch contract" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine

    struct UnregisteredFormulation<:Grammar.AbstractFormulation end

    @test LineCableModels.FormulationOptions === Grammar.FormulationOptions === NamedTuple
    @test LineCableModels.ComputationOptions === Grammar.ComputationOptions === NamedTuple
    @test LineCableModels.formulation_options === Grammar.formulation_options
    @test LineCableModels.computation_options === Grammar.computation_options
    @test parentmodule(Grammar.formulation_options) === Grammar
    @test parentmodule(Grammar.computation_options) === Grammar

    owner=Val(AnalyticalFormulation)
    @test hasmethod(
        Grammar.formulation_options,
        Tuple{typeof(owner), NamedTuple}
    )
    @test hasmethod(
        Grammar.computation_options,
        Tuple{typeof(owner), NamedTuple}
    )
    @test_throws MethodError Grammar.formulation_options((;))
    @test_throws MethodError Grammar.computation_options((;))
    @test_throws MethodError Grammar.formulation_options(Val(:analytical), (;))
    @test_throws MethodError Grammar.computation_options(Val(:analytical), (;))
    @test_throws MethodError Grammar.formulation_options(
        Val(UnregisteredFormulation),
        (;)
    )
    @test_throws MethodError Grammar.computation_options(
        Val(UnregisteredFormulation),
        (;)
    )
    @test_throws MethodError Grammar.formulation_options(owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(owner, nothing)

    formulation=@inferred Grammar.formulation_options(owner, (;))
    @test formulation == (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = true,
        temperature_correction = true,
        output = Val(:parameters)
    )
    traced=Grammar.formulation_options(owner, (output = :trace,))
    @test traced.output == Val(:trace)
    @test_throws ArgumentError Grammar.formulation_options(owner, (unknown = true,))
    @test_throws ArgumentError Grammar.formulation_options(owner, (output = :unknown,))

    default_execution=@inferred Grammar.computation_options(owner, (;))
    @test default_execution.output_basis == Val(:pul)
    execution=Grammar.computation_options(
        owner, (
            verbosity = (default = 1, NLsolve = 0),
            output_basis = :total
        ))
    @test execution == (
        verbosity = (default = 1, NLsolve = 0),
        output_basis = Val(:total)
    )
    @test Engine.verbosity(execution, :NLsolve) == 0
    @test Engine.verbosity(execution, :unlisted) == 1
    @test_throws ArgumentError Grammar.computation_options(owner, (unknown = true,))
    @test_throws ArgumentError Grammar.computation_options(
        owner,
        (output_basis = :unknown,)
    )
    for retired_basis in (:per_length, :per_lenght, :per_unit_length)
        @test_throws ArgumentError Grammar.computation_options(
            owner,
            (output_basis = retired_basis,)
        )
    end
end
