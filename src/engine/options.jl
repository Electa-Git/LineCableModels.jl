function formulation_options(
        ::Val{LineParametersFormulation},
        options::NamedTuple
)::FormulationOptions
    allowed = (
        :reduce_bundle,
        :kron_reduction,
        :ideal_transposition,
        :temperature_correction
    )
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown line-parameter formulation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = true,
            temperature_correction = true
        ),
        options
    )
    all(name -> getproperty(normalized, name) isa Bool,
        (:reduce_bundle, :kron_reduction, :ideal_transposition,
            :temperature_correction)) || throw(ArgumentError(
        "reduction, transposition, and temperature-correction options must be Bool",
    ))
    return (
        reduce_bundle = normalized.reduce_bundle,
        kron_reduction = normalized.kron_reduction,
        ideal_transposition = normalized.ideal_transposition,
        temperature_correction = normalized.temperature_correction
    )
end

function computation_options(
        ::Val{LineCableModelsEngine},
        options::NamedTuple
)::ComputationOptions
    allowed = (:verbosity, :output_basis, :trace)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LineCableModelsEngine computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (verbosity = (default = 0,), output_basis = :pul, trace = false),
        options
    )
    verbosity_values = normalized.verbosity
    verbosity_values isa NamedTuple || throw(ArgumentError(
        "verbosity must be a named tuple",
    ))
    haskey(verbosity_values, :default) || throw(ArgumentError(
        "verbosity must define a default level",
    ))
    all(value -> value isa Integer && value in 0:2, values(verbosity_values)) ||
        throw(ArgumentError("verbosity levels must be integers from 0 to 2"))
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    normalized.trace isa Bool || throw(ArgumentError("trace must be Bool"))
    levels = NamedTuple{keys(verbosity_values)}(Int.(values(verbosity_values)))
    return (
        verbosity = levels,
        output_basis = Val(basis_value),
        trace = Val(normalized.trace)
    )
end

function verbosity(options::NamedTuple, key::Symbol)
    haskey(options, :verbosity) || throw(ArgumentError(
        "computation options do not define verbosity",
    ))
    return get(options.verbosity, key, options.verbosity.default)
end
