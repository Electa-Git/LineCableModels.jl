@testitem "Engine / option grammar / owner dispatch contract" tags=[:unit] setup=[
    UseEngineSupport
] begin
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine

    struct UnregisteredFormulation<:Grammar.AbstractFormulation end

    @test LineCableModels.FormulationOptions === Grammar.FormulationOptions !== NamedTuple
    @test LineCableModels.ComputationOptions === Grammar.ComputationOptions !== NamedTuple
    @test LineCableModels.formulation_options === Grammar.formulation_options
    @test LineCableModels.computation_options === Grammar.computation_options
    @test parentmodule(Grammar.formulation_options) === Grammar
    @test parentmodule(Grammar.computation_options) === Grammar

    formulation_owner=LineParametersFormulation
    computation_type=LineCableModelsCoaxial
    @test hasmethod(
        Grammar.formulation_options,
        Tuple{Type{LineParametersFormulation}, FormulationOptions}
    )
    @test hasmethod(
        Grammar.computation_options,
        Tuple{Type{LineCableModelsCoaxial}, ComputationOptions}
    )
    # Passive selections cannot invoke live normalization by accident.
    @test_throws MethodError Grammar.formulation_options((;))
    retained_options = (reduce_bundle=false, kron_reduction=true,
        ideal_transposition=false)
    @test_throws MethodError Grammar.formulation_options((leaf=formulation_owner =>
        (options=retained_options,),))
    @test_throws MethodError Grammar.computation_options((;))
    @test_throws MethodError Grammar.formulation_options(:analytical, (;))
    @test_throws MethodError Grammar.computation_options(:analytical, ComputationOptions((;)))
    @test_throws MethodError Grammar.formulation_options(
        UnregisteredFormulation,
        (;)
    )
    @test_throws MethodError Grammar.computation_options(
        UnregisteredFormulation, ComputationOptions((;)))
    @test_throws MethodError Grammar.formulation_options(
        formulation_owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(
        computation_type, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(computation_type, nothing)

    formulation=@inferred Grammar.formulation_options(formulation_owner, FormulationOptions())
    @test formulation.data == (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = true
    )
    @test_throws ArgumentError Grammar.formulation_options(
        formulation_owner, FormulationOptions(unknown = true))

    default_execution=@inferred Grammar.computation_options(computation_type, ComputationOptions((;)))
    @test default_execution.data.output_basis == Val(:pul)
    @test default_execution.data.trace == Val(false)
    @test default_execution.data.on_result === nothing
    execution=Grammar.computation_options(
        computation_type, ComputationOptions((
            verbosity = (default = 1, NLsolve = 0),
            output_basis = :total,
            trace = true
        )))
    @test execution.data == (
        verbosity = (default = 1, NLsolve = 0),
        output_basis = Val(:total),
        trace = Val(true),
        on_result = nothing,
        timing = false
    )
    @test Engine.verbosity(execution, :NLsolve) == 0
    @test Engine.verbosity(execution, :unlisted) == 1
    @test_throws ArgumentError Grammar.computation_options(
        computation_type, ComputationOptions((unknown = true,)))
    @test_throws ArgumentError Grammar.computation_options(
        computation_type, ComputationOptions((output_basis = :unknown,)))
    for retired_basis in (:per_length, :per_lenght, :per_unit_length)
        @test_throws ArgumentError Grammar.computation_options(
            computation_type, ComputationOptions((output_basis = retired_basis,)))
    end
end
