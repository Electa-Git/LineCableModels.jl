"""
    validate(subject, context...)

Reject an invalid input before its consumer uses it. Return `subject`, the first
positional argument, unchanged.

`subject` is user input or an input derived from it. The optional context says what the
subject is checked for: its consumer, such as the problem, formulation, physical model or
result being built, or, for the result of a law or record, the function that produced it.
Data that the check needs follows the context. A method returns the subject without
converting or resolving it or building a new record, and without a check when the
consumer does not restrict it. `validate` is the only verb for input checks.

# Errors

- Throws a native exception identifying the invalid field, value, and required
  condition.
- Throws `RequiredInterfaces.NotImplementedError` when a concrete problem type
  does not implement input validation.
"""
function validate end

@required AbstractProblemDefinition begin
    validate(::AbstractProblemDefinition)
end

"""
$(SIGNATURES)

Validate and normalize the options owned by a formulation type.

The implementation that owns `FormulationType` defines a method for the type
itself. Dispatch requires an explicitly supported type. An unregistered formulation
raises `MethodError`.

# Arguments

- `owner`: formulation-type dispatch token.
- `options`: `FormulationOptions` supplied to the defining normalizer. Public
  constructors also accept `options=(...)` shorthand and wrap it at entry.

# Returns

- A formulation-owned [`FormulationOptions`](@ref) record with a fixed set
  of keys for the selected owner.
"""
function formulation_options end

"""
$(TYPEDSIGNATURES)

Allocate a selected formula's reusable arrays during computation initialization.
The arguments are the resolved selection, scalar type, completed numerical input,
plan of fixed indices and geometry, and existing buffer record. Return the extended
record without replacing another owner's storage. Array blocks contain no copied
geometry, material model, selection or validity state. The default uses the existing storage. No material law or integrand is evaluated here.
"""
function initialize_buffers end

function initialize_buffers(
        ::Union{AbstractFormulation, Nothing}, ::Type, input, plan, buffers)
    buffers
end

function initialize_buffers(
        selections::Union{NamedTuple, Tuple}, ::Type{T}, input, plan, buffers) where {T}
    return foldl(values(selections); init = buffers) do accumulated, selected
        initialized = initialize_buffers(selected, T, input, plan, accumulated)
        # Extension methods may append arrays but cannot replace another owner's
        # storage.
        retained = map(values(accumulated),
            values(initialized[keys(accumulated)])) do before, after
            before === after
        end
        all(retained) || throw(ArgumentError(
            "buffer initialization replaced existing storage :$(keys(accumulated)[findfirst(!, retained)])"))
        initialized
    end
end

"""
$(SIGNATURES)

Validate and normalize the options owned by one computation.

The implementation that owns `OwnerType` defines a method for the type itself.
`OwnerType` may identify a core solver or another composite calculation owner.
Dispatch requires an explicitly supported type. An unregistered owner raises `MethodError`.

# Arguments

- `owner`: computation-owner dispatch token.
- `options`: `ComputationOptions` supplied to the defining normalizer. Public
  compute calls also accept `options=(...)` shorthand and wrap it at entry.

# Returns

- A computation-owned [`ComputationOptions`](@ref) record with a fixed set
  of outer keys for the selected owner.
"""
function computation_options end

"""
$(SIGNATURES)

Return supplemental output owned by one core or composite computation.

The formulation type defines a method for itself and returns a
fixed-key [`ComputationDetails`](@ref) record. Dispatch requires an explicitly supported type.
An unregistered formulation raises `MethodError`.
"""
function computation_details end

"""
$(SIGNATURES)

Return the typed supplemental output retained by a completed
result. Result owners define narrow methods beside their result containers.
"""
function details end

"""
$(SIGNATURES)

Calculate a completed result from an explicit problem and formulation.

Concrete solver methods validate and normalize execution `options` through
[`computation_options`](@ref) for their owner. Composite computations may
forward caller options to that solver for validation. Scientific choices use [`formulation_options`](@ref), and supplemental
results use [`computation_details`](@ref). Unsupported problem-formulation
pairs fail through ordinary Julia dispatch.
"""
function compute end

"""
$(SIGNATURES)

Return native numerical values selected from a completed scientific result.

The selector and optional transform are function objects. Result owners define
the supported combinations beside their result representations.
"""
function observe end

function _observe_macro_parts(request)
    valid_request = request isa Expr &&
                    request.head === :ref &&
                    length(request.args) >= 2
    valid_request || throw(ArgumentError(
        "@observe expects indexed `accessor[...]` or `(accessor, transform)[...]`; " *
        "got `$(request)`.",
    ))

    selector = first(request.args)
    indices = request.args[2:end]
    if selector isa Expr && selector.head === :tuple
        length(selector.args) in (2, 3) || throw(ArgumentError(
            "@observe selector tuples require two or three functions; " *
            "got `$(selector)`.",
        ))
        return Tuple(selector.args), indices
    end
    return (selector,), indices
end

"""
    @observe accessor[indices...]
    @observe (accessor, transform)[indices...]
    @observe (statistics, quantity, statistic)[indices...]

Construct a plain observable-request tuple without reading a result.
"""
macro observe(request)
    selectors, indices = _observe_macro_parts(request)
    parts = (selectors..., indices...)
    return Expr(:tuple, map(esc, parts)...)
end

"""
    @observe source accessor[indices...]
    @observe source (accessor, transform)[indices...]
    @observe source (statistics, quantity, statistic)[indices...]

Expand indexed observable syntax into an immediate [`observe`](@ref) call. The
indices follow the selected observable. Ordinary line quantities use row,
column, and sample indices. Diagonal transforms use mode and sample indices.
"""
macro observe(source, request)
    selectors, indices = _observe_macro_parts(request)
    parts = (source, selectors..., indices...)
    escaped = map(esc, parts)
    return :(observe($(escaped...)))
end

"""
$(SIGNATURES)

Publish explicitly requested scientific values for presentation or reporting.

`observables(::Type{T})` declares the selectors supported by `T`.
`observables(source, requests; units, length_unit, frequency_unit,
quantity_units, clip, atol, frequencies)` returns one [`ObservedResult`](@ref),
or an ordinary vector for a result collection. Its four sections retain gridpoint
descriptions, quantities, completed errors, and recorded timings. Quantity-wise
tables are materialized by `ReportBuilder.tabulate` from these detached records. `units` is empty or
positionally aligned with `requests`. With `clip=true` (default), the result
owner's declared native-unit reporting resolution is applied before conversion.
`atol` optionally overrides that resolution. Multiple quantities require keyed
cutoffs. `frequencies` supplies standalone tensor context \\[Hz\\]. These cutoffs
are not certified floating-point error bounds. Each clipped value becomes exact
zero, including its uncertainty. Quantities lacking both a declared resolution and an explicit cutoff remain unchanged. `clip=false` retains raw values in
the requested display units. Absolute/relative error products are never clipped
using their operands' physical cutoffs.
"""
function observables end
