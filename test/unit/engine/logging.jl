@testitem "Logging / verbosity respects module ancestry and caller policy" tags=[:unit] begin
    using Logging
    const LCM = LineCableModels
    @test LCM.verbosity((default = 0, progress = 1, Engine = 2)) ==
          (default = 0, progress = 1, Engine = 2)
    for invalid in ((progress = 1,), (default = 3,), (default = 1.5,), 1)
        @test_throws ArgumentError LCM.verbosity(invalid)
    end
    log = Test.TestLogger(min_level = Logging.Debug)
    with_logger(LCM.VerbosityLogger(log, (
        default = 0, progress = 1, Engine = 2, EarthAdmittance = 0))) do
        @info "progress" _group=:progress
        @info "quiet"
        @debug "ancestor" _module=LCM.Engine.EarthImpedance
        @info "engine quiet" _module=LCM.Engine.EarthAdmittance
        @warn "warning"
        @info "module-less quiet" _module=nothing
        @warn "module-less warning" _module=nothing
        @info "module-less progress" _module=nothing _group=:progress
    end
    @test getproperty.(log.logs, :message) == ["progress", "ancestor", "warning",
        "module-less warning", "module-less progress"]
    # This parent deliberately leaves level filtering to min_enabled_level.
    struct ThresholdLogger <: AbstractLogger
        records::Vector{String}
    end
    Logging.min_enabled_level(::ThresholdLogger) = Logging.Warn
    Logging.catch_exceptions(::ThresholdLogger) = false
    Logging.shouldlog(::ThresholdLogger, level, source, group, id) = group !== :excluded
    Logging.handle_message(logger::ThresholdLogger, level, message, args...;
        kwargs...) = push!(logger.records, string(message))
    parent = ThresholdLogger(String[])
    filter = LCM.VerbosityLogger(parent, (default = 2, progress = 2))
    @test !Logging.catch_exceptions(filter)
    @test !Logging.shouldlog(filter, Logging.Info, @__MODULE__, :progress, :test)
    with_logger(filter) do
        @info "blocked by minimum" _group=:progress
        @warn "blocked by predicate" _group=:excluded
        @warn "accepted"
    end
    @test parent.records == ["accepted"]
    evaluations = Ref(0)
    with_logger(LCM.VerbosityLogger(NullLogger(), (default = 2, progress = 2))) do
        @info "discarded" expensive=(evaluations[]+=1)
    end
    @test evaluations[] == 0
    @test_throws ErrorException with_logger(filter) do
        @warn error("caller exception policy")
    end
end

@testitem "Logging / FEM file destination does not bypass the parent threshold" tags=[:extension] begin
    using Gmsh, Logging
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    parent = Test.TestLogger(min_level = Logging.Warn, respect_maxlog = false)
    file = Test.TestLogger(min_level = Logging.Debug, respect_maxlog = false)
    filter = LineCableModels.VerbosityLogger(parent, (default = 2, progress = 1))
    with_logger(FEM.FEMTeeLogger(filter, file)) do
        @info "file only" _group=:progress
        @debug "file debug"
        @warn "both"
    end
    @test getproperty.(parent.logs, :message) == ["both"]
    @test getproperty.(file.logs, :message) == ["file only", "file debug", "both"]
end
