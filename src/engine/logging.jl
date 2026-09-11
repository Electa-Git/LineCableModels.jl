struct ConsoleVerbosityLogger{L <: AbstractLogger, V <: NamedTuple} <: AbstractLogger
    console::L
    levels::V
end

# AbstractLogger's required methods are not marked public by the Logging stdlib.
Logging.min_enabled_level(::ConsoleVerbosityLogger) = Logging.Debug
Logging.catch_exceptions(::ConsoleVerbosityLogger) = false

function _log_level(level::Integer)
    level == 0 ? Logging.Warn :
    level == 1 ? Logging.Info : Logging.Debug
end

function _log_key(source::Module)
    path = Base.fullname(source)
    return isempty(path) ? :default : first(path)
end

# AbstractLogger's required methods are not marked public by the Logging stdlib.
function Logging.shouldlog(
        logger::ConsoleVerbosityLogger,
        level,
        source,
        group,
        id
)
    performance_sample_active() && return false
    selected = get(logger.levels, _log_key(source), logger.levels.default)
    return level >= _log_level(selected) &&
           Logging.shouldlog(logger.console, level, source, group, id)
end

function Logging.handle_message(
        logger::ConsoleVerbosityLogger,
        level,
        message,
        source,
        group,
        id,
        file,
        line;
        kwargs...
)
    return with_progress_output() do
        Logging.handle_message(
            logger.console, level, message, source, group, id, file, line; kwargs...
        )
    end
end
