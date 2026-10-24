@testitem "Engine / option grammar / owner dispatch" tags=[:unit, :engine] setup=[
    UseEngineSupport
] begin
    const Commons=LineCableModels.Commons
    const Engine=LineCableModels.Engine

    struct UnregisteredFormulation<:Commons.AbstractFormulation end

    @test LineCableModels.FormulationOptions === Commons.FormulationOptions !== NamedTuple
    @test LineCableModels.ComputationOptions === Commons.ComputationOptions !== NamedTuple
    @test LineCableModels.formulation_options === Commons.formulation_options
    @test LineCableModels.computation_options === Commons.computation_options
    @test parentmodule(Commons.formulation_options) === Commons
    @test parentmodule(Commons.computation_options) === Commons

    formulation_owner=LineParametersFormulation
    computation_type=LineCableModelsCoaxial
    @test hasmethod(
        Commons.formulation_options,
        Tuple{Type{LineParametersFormulation}, FormulationOptions}
    )
    @test hasmethod(
        Commons.computation_options,
        Tuple{Type{LineCableModelsCoaxial}, ComputationOptions}
    )
    # Passive selections cannot invoke live normalization by accident.
    @test_throws MethodError Commons.formulation_options((;))
    retained_options = (reduce_bundle=false, kron_reduction=true,
        ideal_transposition=false)
    @test_throws MethodError Commons.formulation_options((leaf=formulation_owner =>
        (options=retained_options,),))
    @test_throws MethodError Commons.computation_options((;))
    @test_throws MethodError Commons.formulation_options(:unregistered, (;))
    @test_throws MethodError Commons.computation_options(:unregistered, ComputationOptions((;)))
    @test_throws MethodError Commons.formulation_options(
        UnregisteredFormulation,
        (;)
    )
    # A formulation type without its own normalizer is refused, with the extension point.
    @test_throws r"computation_options\(::Type\{<:T\}, options\)" Commons.computation_options(
        UnregisteredFormulation, ComputationOptions((;)))
    @test_throws MethodError Commons.formulation_options(
        formulation_owner, Dict{Symbol, Any}())
    @test_throws MethodError Commons.computation_options(
        computation_type, Dict{Symbol, Any}())
    @test_throws MethodError Commons.computation_options(computation_type, nothing)

    formulation=@inferred Commons.formulation_options(formulation_owner, FormulationOptions())
    @test formulation.data == (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = false
    )
    @test_throws ArgumentError Commons.formulation_options(
        formulation_owner, FormulationOptions(unknown = true))

    default_execution=@inferred Commons.computation_options(computation_type, ComputationOptions((;)))
    @test !haskey(default_execution.data, :output_basis)
    @test default_execution.data == (trace = Val(false),)
    execution=Commons.computation_options(computation_type, ComputationOptions((trace = true,)))
    @test execution.data == (trace = Val(true),)
    # The verbosity and the callable of a computation are options of its run, which the
    # run normalizes, whatever the backend.
    for run_option in ((verbosity = (default = 1, NLsolve = 0),), (on_result = identity,))
        @test_throws ArgumentError Commons.computation_options(
            computation_type, ComputationOptions(run_option))
    end
    run=Commons.computation_options(AbstractTraversal,
        ComputationOptions((verbosity = (default = 1, NLsolve = 0),)))
    @test Engine.verbosity(run, :NLsolve) == 0
    @test Engine.verbosity(run, :unlisted) == 1
    @test_throws ArgumentError Commons.computation_options(
        computation_type, ComputationOptions((unknown = true,)))
    @test_throws ArgumentError Commons.computation_options(
        computation_type, ComputationOptions((output_basis = :total,)))
end
