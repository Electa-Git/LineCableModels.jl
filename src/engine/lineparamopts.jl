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

    "Whether material resistivity is corrected to the operating temperature."
    temperature_correction::Bool=true
end

"""
$(TYPEDEF)

Execution options shared by ordinary, full-parametric, and Monte Carlo runs.
`output_basis=:total` scales computed line impedances and admittances by the
materialized system length; cable constants have no system length and therefore
support only `:per_length`.

$(TYPEDFIELDS)
"""
struct ComputeOptions
    "Whether unreduced impedance and admittance matrices are retained."
    store_primitive_matrices::Bool

    "Logging verbosity."
    verbosity::Int

    "Optional path for calculation logs."
    logfile::Union{Nothing,String}

    "Output basis: `:per_length` or `:total`."
    output_basis::Symbol

    function ComputeOptions(;
        store_primitive_matrices::Bool=true,
        verbosity::Integer=0,
        logfile::Union{Nothing,AbstractString}=nothing,
        output_basis::Symbol=:per_length,
    )
        output_basis in (:per_length, :total) || throw(ArgumentError(
            "output_basis must be :per_length or :total; got :$output_basis",
        ))
        return new(
            store_primitive_matrices,
            Int(verbosity),
            logfile === nothing ? nothing : String(logfile),
            output_basis,
        )
    end
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
