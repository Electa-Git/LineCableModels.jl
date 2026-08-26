@testitem "ParametricBuilder / projection / template method and result space" tags=[:unit] begin
    using Statistics

    const Grammar=LineCableModels.Grammar
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct ProjectionPayload
        x::Float64
        y::Float64
    end

    struct ExternalResultSpace{T} <: Grammar.AbstractResultSpace{T}
        values::Vector{T}
    end
    Base.length(source::ExternalResultSpace) = length(source.values)
    Base.size(source::ExternalResultSpace) = (length(source),)
    Base.getindex(source::ExternalResultSpace, index::Integer) = source.values[index]
    Base.iterate(source::ExternalResultSpace, state...) = iterate(source.values, state...)
    Base.firstindex(source::ExternalResultSpace) = firstindex(source.values)
    Base.lastindex(source::ExternalResultSpace) = lastindex(source.values)

    struct DownstreamProblem <: Grammar.AbstractProblemDefinition
        x::Float64
        y::Float64
    end

    struct StageProjection <: PB.AbstractProjectionDefinition
        calls::Vector{Symbol}
    end
    function PB.entitle(definition::StageProjection, source::ExternalResultSpace)
        push!(definition.calls, :entitle)
        return source
    end
    function PB.select(definition::StageProjection, source::ExternalResultSpace)
        push!(definition.calls, :select)
        return (results = collect(source),)
    end
    function PB.derive(definition::StageProjection, selected::NamedTuple)
        push!(definition.calls, :derive)
        return selected.results
    end
    function PB.materialize(definition::StageProjection, value::ProjectionPayload)
        push!(definition.calls, :materialize)
        return DownstreamProblem(value.x, value.y)
    end
    function PB.finish(
            definition::StageProjection,
            problems::Vector{DownstreamProblem}
    )
        push!(definition.calls, :finish)
        return invoke(
            PB.finish,
            Tuple{PB.AbstractProjectionDefinition, Vector{DownstreamProblem}},
            definition,
            problems
        )
    end

    calls=Symbol[]
    source=ExternalResultSpace(ProjectionPayload[ProjectionPayload(1.0, 10.0)])
    projected=PB.project(StageProjection(calls), source)
    @test calls == [:entitle, :select, :derive, :materialize, :finish]
    @test projected isa PB.Gridspace{DownstreamProblem}
    @test eltype(projected) === DownstreamProblem
    @test collect(projected) == [DownstreamProblem(1.0, 10.0)]

    unsupported_calls=Symbol[]
    unsupported=PB.ParametricResult(nothing, ProjectionPayload[ProjectionPayload(1.0, 1.0)])
    @test_throws MethodError PB.project(StageProjection(unsupported_calls), unsupported)
    @test isempty(unsupported_calls)

    struct PreserveProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::PreserveProjection, source::ExternalResultSpace) = source
    PB.select(::PreserveProjection, source::ExternalResultSpace) = (results = collect(source),)
    PB.derive(::PreserveProjection, selected::NamedTuple) = selected.results
    PB.materialize(::PreserveProjection, value::ProjectionPayload) = DownstreamProblem(value.x, value.y)
    @test parentmodule(which(
        PB.finish,
        (PreserveProjection, Vector{DownstreamProblem})
    )) === PB

    multiple=ExternalResultSpace(ProjectionPayload[
        ProjectionPayload(1.0, 10.0),
        ProjectionPayload(2.0, 20.0)
    ])
    selected=PB.select(PreserveProjection(), multiple)
    @test selected isa @NamedTuple{results::Vector{ProjectionPayload}}
    @test collect(PB.project(PreserveProjection(), multiple)) == [
        DownstreamProblem(1.0, 10.0),
        DownstreamProblem(2.0, 20.0)
    ]

    struct ReduceProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::ReduceProjection, source::ExternalResultSpace) = source
    PB.select(::ReduceProjection, source::ExternalResultSpace) = (results = collect(source),)
    function PB.derive(::ReduceProjection, selected::NamedTuple)
        return (ProjectionPayload(
            mean(value.x for value in selected.results),
            mean(value.y for value in selected.results)
        ),)
    end
    PB.materialize(::ReduceProjection, value::ProjectionPayload) = DownstreamProblem(value.x, value.y)
    @test collect(PB.project(ReduceProjection(), multiple)) == [
        DownstreamProblem(1.5, 15.0),
    ]

    struct ExpandProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::ExpandProjection, source::ExternalResultSpace) = source
    PB.select(::ExpandProjection, source::ExternalResultSpace) = (results = collect(source),)
    function PB.derive(::ExpandProjection, selected::NamedTuple)
        return (
            ProjectionPayload(value.x, sign * value.y)
        for value in selected.results for sign in (-1.0, 1.0)
        )
    end
    PB.materialize(::ExpandProjection, value::ProjectionPayload) = DownstreamProblem(value.x, value.y)
    @test collect(PB.project(ExpandProjection(), multiple)) == [
        DownstreamProblem(1.0, -10.0),
        DownstreamProblem(1.0, 10.0),
        DownstreamProblem(2.0, -20.0),
        DownstreamProblem(2.0, 20.0)
    ]

    struct EmptyProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::EmptyProjection, source::ExternalResultSpace) = source
    PB.select(::EmptyProjection, source::ExternalResultSpace) = source
    PB.derive(::EmptyProjection, source::ExternalResultSpace) = ()
    @test_throws ArgumentError PB.project(EmptyProjection(), source)

    struct InconsistentProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::InconsistentProjection, source::ExternalResultSpace) = source
    PB.select(::InconsistentProjection, source::ExternalResultSpace) = source
    PB.derive(::InconsistentProjection, source::ExternalResultSpace) = (1, 2)
    PB.materialize(::InconsistentProjection, value::Int) = value == 1 ?
                                                           DownstreamProblem(1.0, 1.0) :
                                                           :wrong_type
    @test_throws ArgumentError PB.project(InconsistentProjection(), source)

    struct EnvelopeProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::EnvelopeProjection, source::ExternalResultSpace) = source
    PB.select(::EnvelopeProjection, source::ExternalResultSpace) = source
    PB.derive(::EnvelopeProjection, source::ExternalResultSpace) = (only(source),)
    PB.materialize(::EnvelopeProjection, value::ProjectionPayload) = PB.ParametricResult(
        nothing, ProjectionPayload[value])
    @test_throws ArgumentError PB.project(EnvelopeProjection(), source)

    struct JointStatistics
        x_mean::Float64
        y_mean::Float64
    end
    statistical=UQ.MonteCarloResult(
        nothing,
        ProjectionPayload[ProjectionPayload(0.0, 0.0)],
        JointStatistics[JointStatistics(1.5, 150.0)],
        [[(1.0, 100.0), (2.0, 200.0), (3.0, 300.0)]],
        nothing,
        UInt64(1),
        UInt64[2],
        Int[3]
    )

    struct SummaryProjection <: PB.AbstractProjectionDefinition end
    PB.entitle(::SummaryProjection, source::UQ.MonteCarloResult{ProjectionPayload}) = source
    PB.select(::SummaryProjection, source::UQ.MonteCarloResult) = (statistics = UQ.statistics(source),)
    PB.derive(::SummaryProjection, selected::NamedTuple) = (
        ProjectionPayload(value.x_mean, value.y_mean) for value in selected.statistics
    )
    PB.materialize(::SummaryProjection, value::ProjectionPayload) = DownstreamProblem(value.x, value.y)
    @test collect(PB.project(SummaryProjection(), statistical)) == [
        DownstreamProblem(1.5, 150.0),
    ]

    struct RetainedTrialProjection{I <: Tuple} <: PB.AbstractProjectionDefinition
        indices::I
    end
    function PB.entitle(
            ::RetainedTrialProjection,
            source::UQ.MonteCarloResult{ProjectionPayload}
    )
        UQ.samples(source) === nothing && throw(ArgumentError("retained samples required"))
        length(UQ.samples(source)) == length(source) || throw(DimensionMismatch(
            "retained samples must align with stored core results",
        ))
        return source
    end
    PB.select(::RetainedTrialProjection, source::UQ.MonteCarloResult) = (samples = UQ.samples(source),)
    function PB.derive(definition::RetainedTrialProjection, selected::NamedTuple)
        joint_samples = only(selected.samples)
        return (joint_samples[index] for index in definition.indices)
    end
    PB.materialize(::RetainedTrialProjection, value::Tuple{
        Float64, Float64}) = DownstreamProblem(value...)

    retained=collect(PB.project(RetainedTrialProjection((1, 3)), statistical))
    @test retained == [
        DownstreamProblem(1.0, 100.0),
        DownstreamProblem(3.0, 300.0)
    ]
    @test all(problem -> problem.y == 100 * problem.x, retained)

    @test LineCableModels.project === PB.project
    @test parentmodule(LineCableModels.project) === PB
    @test Base.ispublic(PB, :AbstractProjectionDefinition)
    @test Base.ispublic(PB, :entitle)
    @test Base.ispublic(PB, :select)
    @test Base.ispublic(PB, :derive)
    @test Base.ispublic(PB, :materialize)
    @test !Base.ispublic(PB, :finish)
    @test !isdefined(LineCableModels, :AbstractProjectionDefinition)
    @test !isdefined(LineCableModels, :entitle)
    @test !isdefined(LineCableModels, :select)
    @test !isdefined(LineCableModels, :derive)
    @test !isdefined(LineCableModels, :materialize)
    @test !isdefined(LineCableModels, :finish)
    @test !isdefined(LineCableModels, :ProjectionResult)
    @test !isdefined(PB, :ProjectionResult)
end
