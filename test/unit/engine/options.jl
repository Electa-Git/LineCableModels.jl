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

    formulation_owner=LineParametersFormulation
    computation_type=LineCableModelsCoaxial
    @test hasmethod(
        Grammar.formulation_options,
        Tuple{Type{LineParametersFormulation}, NamedTuple}
    )
    @test hasmethod(
        Grammar.computation_options,
        Tuple{Type{LineCableModelsCoaxial}, NamedTuple}
    )
    @test_throws MethodError Grammar.formulation_options((;))
    @test_throws MethodError Grammar.computation_options((;))
    @test_throws MethodError Grammar.formulation_options(:analytical, (;))
    @test_throws MethodError Grammar.computation_options(:analytical, (;))
    @test_throws MethodError Grammar.formulation_options(
        UnregisteredFormulation,
        (;)
    )
    @test_throws MethodError Grammar.computation_options(
        UnregisteredFormulation,
        (;)
    )
    @test_throws MethodError Grammar.formulation_options(
        formulation_owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(
        computation_type, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(computation_type, nothing)

    formulation=@inferred Grammar.formulation_options(formulation_owner, (;))
    @test formulation == (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = true
    )
    @test_throws ArgumentError Grammar.formulation_options(
        formulation_owner, (unknown = true,))

    default_execution=@inferred Grammar.computation_options(computation_type, (;))
    @test default_execution.output_basis == Val(:pul)
    @test default_execution.trace == Val(false)
    @test default_execution.on_result === nothing
    execution=Grammar.computation_options(
        computation_type, (
            verbosity = (default = 1, NLsolve = 0),
            output_basis = :total,
            trace = true
        ))
    @test execution == (
        verbosity = (default = 1, NLsolve = 0),
        output_basis = Val(:total),
        trace = Val(true),
        on_result = nothing
    )
    @test Engine.verbosity(execution, :NLsolve) == 0
    @test Engine.verbosity(execution, :unlisted) == 1
    @test_throws ArgumentError Grammar.computation_options(
        computation_type, (unknown = true,))
    @test_throws ArgumentError Grammar.computation_options(
        computation_type,
        (output_basis = :unknown,)
    )
    for retired_basis in (:per_length, :per_lenght, :per_unit_length)
        @test_throws ArgumentError Grammar.computation_options(
            computation_type,
            (output_basis = retired_basis,)
        )
    end
end

@testitem "Engine / FEM computation option ownership" tags=[:unit] begin
    using LineCableModels

    formulation = LineCableModelsFEM(options=(physics=:quasi_fw,))
    @test formulation.options isa FormulationOptions
    @test formulation.options.physics === Symbol("quasi-fw")
    @test fieldnames(typeof(formulation)) == (:methods, :options, :definitions)
    @test !haskey(NamedTuple(formulation), :execution)
    # Execution settings cannot enter through the formulation, or vice versa.
    @test_throws ArgumentError LineCableModelsFEM(options=(frequency_workers=4,))
    @test_throws MethodError LineCableModelsFEM(fem_options=(ui=true,))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, (physics=:quasi_fw,))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, (unknown=true,))
    defaults = @inferred computation_options(LineCableModelsFEM, (;))
    @test defaults isa ComputationOptions
    callback = (problem, index, result) -> nothing
    raw = (frequency_workers=Int32(4), solver_threads=Int16(2),
        on_result=callback, trace=true, mesh_policy=:remesh,
        getdp_executable=SubString("/tmp/getdp", 1), verbosity=(default=1,))
    configured = computation_options(LineCableModelsFEM, raw)
    @test keys(configured) == keys(defaults)
    @test isconcretetype(typeof(configured))
    @test configured.frequency_workers isa Int
    @test configured.solver_threads isa Int
    @test configured.getdp_executable isa String
    @test configured.on_result === callback
    @test configured.trace === Val(true)
    resources(options) = options.frequency_workers * options.solver_threads
    @test (@inferred resources(configured)) == 8
    @test formulation.options.physics === Symbol("quasi-fw")
    for invalid in ((ui=1,), (plot_field_maps=:yes,), (keep_run_directory=1,),
        (mesh_path="",), (getdp_executable=1,), (log_file="",),
        (frequency_workers=true,), (solver_threads=1.5,), (gmsh_verbosity=true,),
        (getdp_verbosity=6,), (resume_run_directory="",), (resume_run_directory=:bad,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM, invalid)
    end
end
