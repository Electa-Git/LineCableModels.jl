# Computation frequencies of the PSCAD phase-scan adapter.
function validate(values::AbstractVector, ::Type{<:PSCADFormulation})
    Engine.has_uncertainty_type(eltype(values)) && throw(ArgumentError(
        "PSCAD does not support Measurement values. Supply deterministic inputs."))
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
