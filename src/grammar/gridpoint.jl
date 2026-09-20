"""
$(TYPEDSIGNATURES)

Copy retained records and arrays without replacing uncertainty-source identities.
Numerical leaves keep their scalar representation; no problem or formulation is
reconstructed. Structured scientific owners supply their own conversion.
"""
detach(value::Union{Number,Symbol,Nothing,Missing,Type,Val,UUIDs.UUID}) = value
detach(value::AbstractString) = String(value)
detach(value::NamedTuple) = map(detach, value)
detach(value::Tuple) = map(detach, value)
detach(value::AbstractArray) = map(detach, value)
detach(value::Pair) = detach(first(value)) => detach(last(value))
detach(value::AbstractDict) = Dict(detach(key) => detach(item) for (key, item) in value)
detach(value::Union{Units.UnitExpr,Units.Unit,Units.Quantity}) = value
function detach(value::Function)
    Base.issingletontype(typeof(value)) || throw(ArgumentError("observation records cannot retain a closure"))
    return value
end
detach(value::Base.Fix2) = Base.Fix2(detach(value.f),detach(value.x))

"""
$(TYPEDSIGNATURES)

Read the retained physical description and original identity of a completed
gridpoint. Result owners extend this method; consumers never inspect a lazy
computation axis or materialize a problem to recover its description.

The fallback identifies an externally supplied result with unspecified inputs.
Such a record cannot establish physical equivalence with another observation.
"""
observation_gridpoint(source) = (id=nothing, inputs=nothing, formulations=nothing,
    coordinates=nothing, uncertainty=nothing, missing_reason=:physical_inputs_not_supplied)

"""
$(TYPEDSIGNATURES)

Identify a completed physical point and formulation within one calculation.
The default source UUID uses system randomness independently of scientific RNGs.
Collection owners share `source_id` and supply the original one-based indices.
"""
gridpoint_id(; source_id=UUIDs.uuid4(Random.RandomDevice()),
    problem_index::Integer=1, formulation_index::Integer=1) =
    (; source_id, problem_index=Int(problem_index), formulation_index=Int(formulation_index))
