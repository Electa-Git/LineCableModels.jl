function UniversalLineModel.LineFrequencyData(segment::PropagationParameters)
    diagnostics = details(segment).data.modal.diagnostics
    coupling = (z_coupling=maximum(diagnostics.z_coupling),
                y_coupling=maximum(diagnostics.y_coupling))
    largest, index = findmax(values(coupling))
    if largest > 1e-6
        @warn "Phase admittance and current propagation use the diagonal modal approximation" diagnostic=keys(coupling)[index] coupling=largest
    end
    f = frequencies(segment)
    ℓ = line_length(segment)
    admittance = Yc(segment, PhaseDomain)
    propagation = H(segment, PhaseDomain; field=:current)
    transformation = Ti(segment)
    return LineFrequencyData(f, ComplexF64.(2π * im .* f),
        [Matrix(admittance[:, :, k]) for k in eachindex(f)],
        [Matrix(propagation[:, :, k]) for k in eachindex(f)],
        permutedims(exp.(-gamma(segment) .* ℓ)), permutedims(gamma(segment)),
        [Matrix(transformation[:, :, k]) for k in eachindex(f)], ℓ)
end
