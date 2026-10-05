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
            ideal_transposition = true
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

function computation_options(
        ::Type{LineCableModelsCoaxial},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    allowed = (:verbosity, :output_basis, :trace, :on_result, :timing)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsCoaxial computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (verbosity = (default = 0,), output_basis = :pul,
            trace = false, on_result = nothing, timing = false),
        options
    )
    levels = verbosity(normalized.verbosity)
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    normalized.trace isa Bool || throw(ArgumentError("trace must be Bool"))
    normalized.timing isa Bool || throw(ArgumentError("timing must be Bool"))
    return ComputationOptions(;
        verbosity = levels,
        output_basis = Val(basis_value),
        trace = Val(normalized.trace),
        on_result = normalized.on_result,
        timing = normalized.timing
    )
end


description(::Type{<:LineParametersFormulation},
    ::Val{:reduce_bundle},value::Bool;compact::Bool=false) = "bundle reduction="*string(value)
description(::Type{<:LineParametersFormulation},
    ::Val{:kron_reduction},value::Bool;compact::Bool=false) = "Kron reduction="*string(value)
description(::Type{<:LineParametersFormulation},
    ::Val{:ideal_transposition},value::Bool;compact::Bool=false) = "ideal transposition="*string(value)
