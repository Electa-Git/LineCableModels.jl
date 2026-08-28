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

    formulation_owner=Val(LineParametersFormulation)
    computation_owner=Val(LineCableModelsEngine)
    @test hasmethod(
        Grammar.formulation_options,
        Tuple{typeof(formulation_owner), NamedTuple}
    )
    @test hasmethod(
        Grammar.computation_options,
        Tuple{typeof(computation_owner), NamedTuple}
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
    @test_throws MethodError Grammar.formulation_options(
        formulation_owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(
        computation_owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(computation_owner, nothing)

    formulation=@inferred Grammar.formulation_options(formulation_owner, (;))
    @test formulation == (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = true,
        temperature_correction = true
    )
    @test_throws ArgumentError Grammar.formulation_options(
        formulation_owner, (unknown = true,))

    default_execution=@inferred Grammar.computation_options(computation_owner, (;))
    @test default_execution.output_basis == Val(:pul)
    @test default_execution.trace == Val(false)
    execution=Grammar.computation_options(
        computation_owner, (
            verbosity = (default = 1, NLsolve = 0),
            output_basis = :total,
            trace = true
        ))
    @test execution == (
        verbosity = (default = 1, NLsolve = 0),
        output_basis = Val(:total),
        trace = Val(true)
    )
    @test Engine.verbosity(execution, :NLsolve) == 0
    @test Engine.verbosity(execution, :unlisted) == 1
    @test_throws ArgumentError Grammar.computation_options(
        computation_owner, (unknown = true,))
    @test_throws ArgumentError Grammar.computation_options(
        computation_owner,
        (output_basis = :unknown,)
    )
    for retired_basis in (:per_length, :per_lenght, :per_unit_length)
        @test_throws ArgumentError Grammar.computation_options(
            computation_owner,
            (output_basis = retired_basis,)
        )
    end
end
