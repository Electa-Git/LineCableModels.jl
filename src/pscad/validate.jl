"""
$(TYPEDSIGNATURES)

Check the physical earth inventory accepted by PSCAD model export and execution.
The supported model consists of air and one infinite horizontal soil half-space.
"""
function validate(model::EarthModel, ::Type{<:PSCADFormulation})
    _pscad_deterministic(eltype(model))
    !model.vertical_layers && length(model.layers) == 2 &&
    all(layer -> isinf(layer.thickness), model.layers) || throw(ArgumentError(
        "PSCAD supports physical air and one homogeneous horizontal soil half-space"))
    return model
end

"""
$(TYPEDSIGNATURES)

Check that a line-parameter problem uses the real frequency axis accepted by PSCAD.
"""
function validate(problem::LineParametersProblem, ::Type{<:PSCADFormulation})
    iszero(problem.sigma) || throw(ArgumentError(
        "PSCAD line constants are computed on the real frequency axis. The Laplace abscissa σ must be zero."))
    return problem
end

function _pscad_deterministic(types::Type...)
    any(Engine.has_uncertainty_type, types) && throw(ArgumentError(
        "PSCAD does not support Measurement values. Supply deterministic inputs."))
    return nothing
end

function _pscad_blueprints(system::LineCableSystem)
    _pscad_deterministic(eltype(system))
    T = eltype(system)
    return DataModel.CableBlueprint{T}[DataModel.flatten(design, T)
                                       for design in system.designs]
end

# Computation frequencies of the PSCAD phase-scan adapter.
function validate(values::AbstractVector, ::Type{<:PSCADFormulation})
    _pscad_deterministic(eltype(values))
    isempty(values) && throw(ArgumentError("PSCAD requires computation frequencies"))
    invalid = findall(value -> !isfinite(value) || value < 0.1, values)
    isempty(invalid) || throw(DomainError(values[invalid],
        "PSCAD computation frequencies must be finite and at least 0.1 Hz; invalid indices: $invalid"))
    issorted(values) && allunique(values) || throw(ArgumentError(
        "PSCAD frequencies must be strictly increasing"))
    length(values) - 1 in (100, 200, 500, 1000) || throw(ArgumentError(
        "PSCAD phase-scan adapter supports 101, 201, 501 or 1001 logarithmic samples; requested $(length(values))"))
    expected = 10.0 .^ range(log10(Float64(first(values))),
        log10(Float64(last(values))); length = length(values))
    tolerance = max(32eps(Float64), 8Float64(eps(float(one(eltype(values))))))
    all(isapprox(a, b; rtol = tolerance, atol = 0) for (a, b) in zip(values, expected)) ||
        throw(ArgumentError("PSCAD phase-scan adapter requires logarithmic samples; arbitrary-grid execution is not implemented"))
    return values
end
