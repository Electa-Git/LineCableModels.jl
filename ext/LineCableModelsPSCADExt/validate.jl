"""
$(TYPEDSIGNATURES)

Check the physical earth inventory accepted by PSCAD model export and execution.
The supported model consists of air and one infinite horizontal soil half-space.
"""
function validate(model::EarthModel, ::Val{:pscad})
    validate(model)
    !model.vertical_layers && length(model.layers) == 2 &&
    all(layer -> isinf(layer.thickness), model.layers) || throw(ArgumentError(
        "PSCAD supports physical air and one homogeneous horizontal soil half-space"))
    return model
end

function _validate_frequencies(values::AbstractVector)
    isempty(values) && throw(ArgumentError("PSCAD requires calculation frequencies"))
    invalid = findall(value -> !isfinite(value) || value < 0.1, values)
    isempty(invalid) || throw(DomainError(values[invalid],
        "PSCAD calculation frequencies must be finite and at least 0.1 Hz; invalid indices: $invalid"))
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

function validate(problem::LineParametersProblem, formulation::PSCADFormulation)
    pscad_setting(formulation, problem)
    _validate_frequencies(problem.frequencies)
    _pscad_size(problem)
    return problem
end
