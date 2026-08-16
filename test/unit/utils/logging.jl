@testitem "Utils / logging / timestamp, level, and file contracts" tags = [:unit] begin
    using Logging
    const U = LineCableModels.Utils

    buffer = IOBuffer()
    delegate = SimpleLogger(buffer, Logging.Info)
    logger = U.TimestampLogger(delegate)

    @test Logging.min_enabled_level(logger) === Logging.Info
    for level in (Logging.Debug, Logging.Info)
        @test Logging.shouldlog(logger, level, @__MODULE__, :contract, :message) ==
              Logging.shouldlog(delegate, level, @__MODULE__, :contract, :message)
    end

    with_logger(logger) do
        @info "timestamp contract" quantity = 3
    end
    rendered = String(take!(buffer))
    @test occursin(r"\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\] timestamp contract", rendered)
    @test occursin("quantity = 3", rendered)

    previous_logger = global_logger()
    mktempdir() do directory
        logfile = joinpath(directory, "linecablemodels.log")
        console_logfile = joinpath(directory, "console.log")
        open(console_logfile, "w+") do console_output
            redirect_stderr(console_output) do
                try
                    U.set_verbosity!(1, logfile)
                    @debug "filtered file contract"
                    @info "persisted file contract"

                    contents = read(logfile, String)
                    @test occursin("persisted file contract", contents)
                    @test !occursin("filtered file contract", contents)
                finally
                    global_logger(previous_logger)
                end
            end
            seekstart(console_output)
            @test occursin("persisted file contract", read(console_output, String))
        end

        missing_logfile = joinpath(directory, "missing", "linecablemodels.log")
        try
            @test_logs (:warn, r"Failed to set up file logging") U.set_verbosity!(
                1, missing_logfile)
        finally
            global_logger(previous_logger)
        end
        @test !isfile(missing_logfile)
    end
end
