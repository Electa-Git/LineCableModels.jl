"""
$(TYPEDSIGNATURES)

Normalize the line-parameter reductions. `reduce_bundle=true` merges the
conductors of each active phase, `kron_reduction=true` eliminates conductors with
phase zero, and `ideal_transposition=false` leaves the retained matrices
untransposed. With `ideal_transposition=true`, the retained ``Z`` and
potential-coefficient matrix ``P`` are averaged over cyclic transposition before
the inversion of ``P`` to ``Y``.
"""
function formulation_options(
        ::Type{LineParametersFormulation},
        record::FormulationOptions
)::FormulationOptions
    options = record.data
    allowed = (
        :reduce_bundle,
        :kron_reduction,
        :ideal_transposition
    )
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown line-parameter formulation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = false
        ),
        options
    )
    all(name -> getproperty(normalized, name) isa Bool,
        (:reduce_bundle, :kron_reduction, :ideal_transposition)) || throw(ArgumentError(
        "reduction and transposition options must be Bool",
    ))
    return FormulationOptions(;
        reduce_bundle = normalized.reduce_bundle,
        kron_reduction = normalized.kron_reduction,
        ideal_transposition = normalized.ideal_transposition
    )
end

# The coaxial backend's one computation option is `trace`. The verbosity and the callable of
# a computation are options of its run.
function computation_options(
        ::Type{LineCableModelsCoaxial},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    unknown = filter(key -> key !== :trace, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsCoaxial computation options: $(sort!(collect(unknown)))",
    ))
    trace = get(options, :trace, false)
    trace isa Bool || throw(ArgumentError("trace must be Bool"))
    trace && return ComputationOptions((trace = Val(true),))
    return ComputationOptions((trace = Val(false),))
end

# The coaxial backend normalizes the computation options of a line-parameter formulation.
computation_options(::Type{<:LineParametersFormulation}, record::ComputationOptions) =
    computation_options(LineCableModelsCoaxial, record)

# A normalized record stores the trace switch as a type. This method returns it as it is.
computation_options(::Type{LineCableModelsCoaxial},
    record::ComputationOptions{<:NamedTuple{(:trace,), <:Tuple{Val}}}) = record


description(::Type{<:LineParametersFormulation},
    ::Val{:reduce_bundle},value::Bool;compact::Bool=false) = "bundle reduction="*string(value)
description(::Type{<:LineParametersFormulation},
    ::Val{:kron_reduction},value::Bool;compact::Bool=false) = "Kron reduction="*string(value)
description(::Type{<:LineParametersFormulation},
    ::Val{:ideal_transposition},value::Bool;compact::Bool=false) = "ideal transposition="*string(value)
