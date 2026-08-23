const _FORMULATION_DEFAULTS = (
    reduce_bundle = true,
    kron_reduction = true,
    ideal_transposition = true,
    temperature_correction = true,
    output = :parameters
)

const _COMPUTATION_DEFAULTS = (
    verbosity = (default = 0,),
    output_basis = :per_length
)

_to_namedtuple(options::NamedTuple) = options
_to_namedtuple(options::Base.Pairs) = (; options...)
_to_namedtuple(options::AbstractDict) = (; options...)
_to_namedtuple(::Nothing) = (;)

function _strict_namedtuple(options, defaults::NamedTuple, owner::AbstractString)
    values = _to_namedtuple(options)
    unknown = setdiff(Set(keys(values)), Set(keys(defaults)))
    isempty(unknown) || throw(ArgumentError(
        "unknown $owner options: $(sort!(collect(unknown)))",
    ))
    return merge(defaults, values)
end

function formulation_options(options)::FormulationOptions
    normalized = _strict_namedtuple(options, _FORMULATION_DEFAULTS, "formulation")
    all(name -> getproperty(normalized, name) isa Bool,
        (:reduce_bundle, :kron_reduction, :ideal_transposition,
            :temperature_correction)) || throw(ArgumentError(
        "reduction, transposition, and temperature-correction options must be Bool",
    ))
    output = normalized.output
    output in (:parameters, :trace) || throw(ArgumentError(
        "formulation output must be :parameters or :trace; got $(repr(output))",
    ))
    return merge(normalized, (output = Val(output),))
end

function computation_options(options)::ComputationOptions
    normalized = _strict_namedtuple(options, _COMPUTATION_DEFAULTS, "computation")
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
    basis_value in (:per_length, :total) || throw(ArgumentError(
        "output_basis must be :per_length or :total; got $(repr(basis_value))",
    ))
    levels = NamedTuple{keys(verbosity_values)}(Int.(values(verbosity_values)))
    return (verbosity = levels, output_basis = Val(basis_value))
end

@inline _output_basis(options::NamedTuple) = typeof(options.output_basis).parameters[1]

function verbosity(options::NamedTuple, key::Symbol)
    haskey(options, :verbosity) || throw(ArgumentError(
        "computation options do not define verbosity",
    ))
    return get(options.verbosity, key, options.verbosity.default)
end
