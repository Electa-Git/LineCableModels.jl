@testitem "ParametricBuilder / result transport / owned result spaces" tags=[:unit, :uq] begin
    const Commons=LineCableModels.Commons
    const PB=LineCableModels.ParametricBuilder
    const UQ=LineCableModels.UQ

    struct TransportResult <: Commons.AbstractCoreResult
        value::Float64
    end

    struct DownstreamProblem <: Commons.AbstractProblemDefinition
        source::TransportResult
    end
    LineCableModels.validate(problem::DownstreamProblem) = problem

    struct TransportFormulation <: Commons.AbstractFormulation end

    values=TransportResult[TransportResult(1.0), TransportResult(2.0)]
    source=Commons.ParametricResult(values,
        fieldtype(Commons.ParametricResult, :failures)(), CartesianIndices((2, 1)),
        Commons.ComputationDetails())
    transported=PB.Gridspace{DownstreamProblem}(source)

    @test transported.grids === (source,)
    @test length(transported) == length(source)
    @test eltype(transported) === DownstreamProblem
    @test collect(transported) == DownstreamProblem.(values)

    linear=UQ.LinearErrorResult(values,
        fieldtype(UQ.LinearErrorResult, :failures)(), Commons.ComputationDetails())
    uncertain_transport=PB.Gridspace{DownstreamProblem}(linear)
    @test uncertain_transport.grids === (linear,)
    @test collect(uncertain_transport) == DownstreamProblem.(values)

    # A result space of another owner transports only through its owner's method.
    struct UnsupportedResults <: Commons.AbstractResultSpace{TransportResult} end
    unsupported=UnsupportedResults()
    error=try
        PB.Gridspace{DownstreamProblem}(unsupported)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError
    message=sprint(showerror, error)
    @test occursin("Gridspace transport from UnsupportedResults", message)
    @test occursin("to problem DownstreamProblem", message)
    @test occursin("?Gridspace", message)
    @test occursin("electa-git.github.io/LineCableModels.jl/dev/gridspace/", message)

    monte_carlo=UQ.MonteCarloResult(
        (confidence = 0.95, cdf_tol = 0.02, distribution = :normal),
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
    @test occursin("requires a transport for result type", sprint(showerror,
        monte_error))

end

# Current result-family transport is exercised above. Retired terminology,
# filenames and exact consumer-method inventories do not define the API.
