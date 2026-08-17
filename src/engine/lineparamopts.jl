"""
$(TYPEDEF)

Select the physical reductions and corrections used by an EMT formulation.

$(TYPEDFIELDS)
"""
Base.@kwdef struct EMTOptions <: AbstractFormulationOptions
    "Whether bundled conductors are reduced to an equivalent conductor."
    reduce_bundle::Bool=true
    "Whether unassigned conductors are eliminated by Kron reduction."
    kron_reduction::Bool=true
    "Whether ideal transposition is applied."
    ideal_transposition::Bool=true
    "Whether conductor resistivity is corrected to operating temperature."
    temperature_correction::Bool=true
end

"""
$(TYPEDEF)

Execution settings for ordinary and parametric calculations.

Verbosity is a named tuple whose `default` value applies when no component
key is present. Levels are `0` for warnings, `1` for information, and `2` for
debugging.

$(TYPEDFIELDS)
"""
struct ComputeOptions{V <: NamedTuple}
    "Per-component verbosity settings."
    verbosity::V
    "Output basis: `:per_length` or `:total`."
    output_basis::Symbol

    function ComputeOptions(;
            verbosity::NamedTuple = (default = 0,),
            output_basis::Symbol = :per_length
    )
        haskey(verbosity, :default) || throw(ArgumentError(
            "verbosity must define a default level",
        ))
        all(value -> value isa Integer && value in 0:2, values(verbosity)) ||
            throw(ArgumentError("verbosity levels must be integers from 0 to 2"))
        output_basis in (:per_length, :total) || throw(ArgumentError(
            "output_basis must be :per_length or :total; got :$output_basis",
        ))
        levels = NamedTuple{keys(verbosity)}(Int.(values(verbosity)))
        return new{typeof(levels)}(levels, output_basis)
    end
end

function verbosity(options::ComputeOptions, key::Symbol)
    get(options.verbosity, key, options.verbosity.default)
end

_to_namedtuple(options::NamedTuple) = options
_to_namedtuple(options::Base.Pairs) = (; options...)
_to_namedtuple(options::AbstractDict) = (; options...)
_to_namedtuple(::Nothing) = (;)

function _strict_options(::Type{Options}, options) where {Options}
    options isa Options && return options
    values = _to_namedtuple(options)
    allowed = Set(fieldnames(Options))
    unknown = setdiff(Set(keys(values)), allowed)
    isempty(unknown) || throw(ArgumentError(
        "unknown options for $Options: $(sort!(collect(unknown)))",
    ))
    return Options(; values...)
end

formulation_options(options) = _strict_options(EMTOptions, options)
compute_options(options) = _strict_options(ComputeOptions, options)
