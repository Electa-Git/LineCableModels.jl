@testmodule GmshCleanupFailures begin
    const FEMGmshSession = NamedTuple
    const FEMResolvedModel = NamedTuple
    struct LineCableModelsFEMError <: Exception end
    const primary = LineCableModelsFEMError()
    const calls = Any[]
    const fail_inspection = Ref(false)
    const gmsh = (
        view=(get_tags=() -> [1, 2], remove=tag -> begin
            push!(calls, (:remove_view, tag))
            error("view removal failed")
        end),
        model=(list=() -> ["caller", "temporary"], get_current=() -> "caller",
            add=name -> push!(calls, (:add_model, name)),
            set_current=name -> begin
                push!(calls, (:set_current, name))
                name == "caller" && error("model restoration failed")
            end,
            remove=() -> begin
                push!(calls, :remove_model)
                error("model removal failed")
            end),
        option=(set_number=(name, value) -> push!(calls, (name, value)),),
        merge=path -> push!(calls, (:merge, path)))
    _restore_onelab(snapshot) = push!(calls, (:onelab, snapshot))
    function _inspect_loaded_mesh(model, path)
        push!(calls, :inspect)
        fail_inspection[] && throw(primary)
    end

    # Load the actual cleanup methods against a local Gmsh API fixture.
    for (file, name) in (("compute.jl", :_finish_gmsh), ("mesh.jl", :_validate_mesh_file))
        Base.include(@__MODULE__, joinpath(@__DIR__, "../../ext/LineCableModelsGmshExt", file)) do expression
            expression isa Expr && expression.head === :function &&
                expression.args[1] isa Expr && expression.args[1].head === :call &&
                expression.args[1].args[1] === name ? expression : nothing
        end
    end
end

@testitem "Gmsh FEM / cleanup reports failures and preserves the primary exception" tags=[:extension, :fem] default_imports=false setup=[GmshCleanupFailures] begin
    using Test, Logging
    const G = GmshCleanupFailures
    session = (owned=false, previous_model="caller", initial_models=Set(["caller"]),
        initial_views=Set([1]), terminal_option=1.0, verbosity_option=2.0,
        geometry_tolerance=1e-7, boolean_tolerance=2e-7, onelab=:snapshot)
    logger = Test.TestLogger(min_level=Logging.Warn)
    empty!(G.calls)
    caught = try
        try
            throw(G.primary)
        finally
            with_logger(logger) do
                G._finish_gmsh(session)
            end
        end
    catch exception
        exception
    end
    @test caught === G.primary
    @test G.calls == [(:remove_view, 2), (:set_current, "temporary"), :remove_model,
        (:set_current, "caller"), ("General.Terminal", 1.0), ("General.Verbosity", 2.0),
        ("Geometry.Tolerance", 1e-7), ("Geometry.ToleranceBoolean", 2e-7), (:onelab, :snapshot)]
    @test length(logger.logs) == 3
    @test logger.logs[1].kwargs[:view] == 2
    @test logger.logs[2].kwargs[:model] == "temporary"
    @test logger.logs[3].kwargs[:model] == "caller"
    @test all(record -> first(record.kwargs[:exception]) isa ErrorException, logger.logs)

    mktemp() do path, io
        close(io)
        for fails in (false, true)
            G.fail_inspection[] = fails
            empty!(G.calls)
            logger = Test.TestLogger(min_level=Logging.Warn)
            result = with_logger(logger) do
                try
                    G._validate_mesh_file((;), path)
                catch exception
                    exception
                end
            end
            @test result === (fails ? G.primary : nothing)
            @test G.calls[end-1:end] == [:remove_model, (:set_current, "caller")]
            @test length(logger.logs) == 2
            @test startswith(logger.logs[1].kwargs[:model], "LineCableModelsFEM-validation-")
            @test logger.logs[2].kwargs[:model] == "caller"
            @test all(record -> record.kwargs[:mesh] == path, logger.logs)
            @test all(record -> first(record.kwargs[:exception]) isa ErrorException, logger.logs)
        end
    end
end
