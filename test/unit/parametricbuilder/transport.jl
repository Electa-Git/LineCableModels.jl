@testitem "ParametricBuilder / result transport / owned result spaces" tags=[:unit] begin
    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct TransportResult <: Grammar.AbstractCoreResult
        value::Float64
    end

    struct DownstreamProblem <: Grammar.AbstractProblemDefinition
        source::TransportResult
    end
    LineCableModels.validate(problem::DownstreamProblem) = problem

    struct TransportFormulation <: Grammar.AbstractFormulation end

    values=TransportResult[TransportResult(1.0), TransportResult(2.0)]
    combinatorial=PB.Combinatorial(TransportFormulation())
    source=PB.ParametricResult(combinatorial, values)
    transported=PB.Gridspace{DownstreamProblem}(source)

    @test transported.grids === (source,)
    @test length(transported) == length(source)
    @test eltype(transported) === DownstreamProblem
    @test collect(transported) == DownstreamProblem.(values)

    linear=UQ.LinearErrorResult(
        UQ.LinearError(TransportFormulation()),
        values
    )
    uncertain_transport=PB.Gridspace{DownstreamProblem}(linear)
    @test uncertain_transport.grids === (linear,)
    @test collect(uncertain_transport) == DownstreamProblem.(values)

    unsupported=PB.ParametricResult(TransportFormulation(), values)
    error=try
        PB.Gridspace{DownstreamProblem}(unsupported)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    message=sprint(showerror, error)
    @test occursin("Gridspace transport from ParametricResult", message)
    @test occursin("to problem DownstreamProblem", message)
    @test occursin("?Gridspace", message)
    @test occursin("electa-git.github.io/LineCableModels.jl/dev/gridspace/", message)

    monte_carlo=UQ.MonteCarloResult(
        UQ.MonteCarlo(TransportFormulation(); trials = 1, seed = 1),
        TransportResult[TransportResult(1.0)],
        [1],
        nothing,
        nothing,
        UInt64(1),
        UInt64[1],
        [1]
    )
    monte_error=try
        PB.Gridspace{DownstreamProblem}(monte_carlo)
        nothing
    catch exception
        exception
    end
    @test monte_error isa ArgumentError
    @test occursin("requires a reconstruction for result type", sprint(showerror,
        monte_error))

    @test !isdefined(LineCableModels, :project)
    @test !isdefined(PB, :project)
    @test !isdefined(PB, :AbstractProjectionDefinition)
    @test !isdefined(PB, :entitle)
    @test !isdefined(PB, :select)
    @test !isdefined(PB, :derive)
    @test !isdefined(PB, :finish)
end

@testitem "Architecture / Gridspace exclusively owns result transport" tags=[:quality] begin
    using LineCableModels

    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const Engine=LineCableModels.Engine

    root=pkgdir(LineCableModels)

    # The retired transaction and its staged protocol must not regain a public
    # binding under their old names.
    for name in (
        :project, :AbstractProjectionDefinition, :ProjectionResult,
        :entitle, :select, :derive, :finish
    )
        @test !isdefined(LineCableModels, name)
        @test !isdefined(PB, name)
    end

    for path in (
        joinpath("src", "parametricbuilder", "project.jl"),
        joinpath("test", "unit", "parametricbuilder", "projection.jl")
    )
        @test !ispath(joinpath(root, path))
    end

    # A typed result consumer added to ParametricBuilder is an architectural
    # action. Keep Gridspace's own transport error; reject a second
    # result-to-problem protocol under a novel verb.
    function consumes_result_space(method)
        signature=Base.unwrap_unionall(method.sig)
        return any(signature.parameters[2:end]) do argument
            try
                argument <: Grammar.AbstractResultSpace
            catch
                false
            end
        end
    end

    consumers=Set{Symbol}()
    for name in names(PB; all = true, imported = false)
        isdefined(PB, name) || continue
        methods_for_binding=try
            methods(getfield(PB, name))
        catch
            continue
        end
        any(consumes_result_space, methods_for_binding) && push!(consumers, name)
    end
    @test consumers == Set((:_transport_error,))

    # Result-bearing Gridspace constructors may only be implemented by the
    # result family that owns the transport. Any new owner must amend this gate
    # explicitly instead of introducing a parallel conversion layer.
    representative=PB.Gridspace{Engine.LineParametersProblem}
    transport_methods=filter(consumes_result_space, collect(methods(representative)))
    transport_owners=Set(normpath(abspath(String(method.file)))
    for method in transport_methods)
    core_owners=Set((
        normpath(joinpath(root, "src", "parametricbuilder", "results.jl")),
        normpath(joinpath(root, "src", "uq", "results.jl"))
    ))
    allowed_owners=union(core_owners,
        Set((
            normpath(joinpath(root, "ext", "LineCableModelsMeasurementsExt.jl")),
        )))
    @test core_owners ⊆ transport_owners
    @test transport_owners ⊆ allowed_owners

    # Keep the retired terminology out of the exact source and documentation
    # surfaces that own result transport. Other uses of "project" elsewhere in
    # the repository (package projects, CAD projections) are unrelated.
    surfaces=String[
        joinpath(root, "src", "LineCableModels.jl"),
        joinpath(
            root, "ext", "LineCableModelsMeasurementsExt.jl"),
        joinpath(
            root, "docs", "src", "gridspace.md"),
        joinpath(
            root, "docs", "src", "developers.md"),
        joinpath(
            root, "docs", "src", "extensions.md")
    ]
    for directory in (
        joinpath(root, "src", "parametricbuilder"),
        joinpath(root, "src", "uq")
    )
        append!(surfaces,
            [joinpath(path, file)
             for (path, _, files) in walkdir(directory)
             for file in files
             if endswith(file, ".jl")])
    end
    retired_pattern=r"(?i)\bproject(?:ion|ed|ing)?\b|AbstractProjectionDefinition|ProjectionResult"
    offenders=filter(path->occursin(retired_pattern, read(path, String)), surfaces)
    @test isempty(offenders)
end
