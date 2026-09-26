"""
$(TYPEDSIGNATURES)

Validate verbosity levels, or select a level from execution options. Levels
0, 1, and 2 permit warnings, information, and debug messages respectively.
The `progress` group uses its explicit level or `default`; other messages use
the nearest explicitly configured module ancestor before falling back to `default`.
"""
function verbosity(levels::NamedTuple)
    haskey(levels, :default) ||
        throw(ArgumentError("verbosity must define a default level"))
    all(value -> value isa Integer && value in 0:2, values(levels)) ||
        throw(ArgumentError("verbosity levels must be integers from 0 to 2"))
    return NamedTuple{keys(levels)}(Int.(values(levels)))
end
verbosity(levels) = throw(ArgumentError("verbosity must be a named tuple"))
function verbosity(record::ComputationOptions, key::Symbol)
    levels = get(record.data, :verbosity, (default = 0,))
    return get(levels, key, levels.default)
end
function verbosity(levels::NamedTuple, source::Union{Module, Nothing}, group)
    group === :progress && return get(levels, :progress, levels.default)
    source === nothing && return levels.default
    while true
        haskey(levels, nameof(source)) && return levels[nameof(source)]
        ancestor = parentmodule(source)
        ancestor === source && return levels.default
        source = ancestor
    end
end

"""
$(TYPEDEF)

Filter ordinary Julia log records by execution verbosity and forward accepted
records to the caller's logger. The parent logger retains its filtering and
exception policy. This filter owns no progress counters or output resources.

$(TYPEDFIELDS)
"""
struct VerbosityLogger{L <: Logging.AbstractLogger, V <: NamedTuple} <:
       Logging.AbstractLogger
    "Caller-owned destination."
    parent::L
    "Validated verbosity levels."
    levels::V
end

# Julia's logger protocol methods are not public stdlib bindings.
function Logging.min_enabled_level(logger::VerbosityLogger)
    Logging.min_enabled_level(logger.parent)
end
Logging.catch_exceptions(logger::VerbosityLogger) = Logging.catch_exceptions(logger.parent)
function Logging.shouldlog(logger::VerbosityLogger, level, source, group, id)
    selected = verbosity(logger.levels, source, group)
    threshold = selected == 0 ? Logging.Warn : selected == 1 ? Logging.Info : Logging.Debug
    return level >= threshold && level >= Logging.min_enabled_level(logger.parent) &&
           Logging.shouldlog(logger.parent, level, source, group, id)
end
function Logging.handle_message(logger::VerbosityLogger, level, message, source, group,
        id, file, line; kwargs...)
    return Logging.handle_message(logger.parent, level, message, source, group, id,
        file, line; kwargs...)
end

public verbosity, VerbosityLogger
