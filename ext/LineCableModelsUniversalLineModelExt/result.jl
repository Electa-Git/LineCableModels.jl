function WidebandModel(model::ULMModel)
    admittance = (poles=copy(model.pY), residues=stack(model.RY), direct=copy(model.DY))
    propagation = [(delay=model.tau[g], poles=copy(model.pH[g]),
        residues=stack(model.RH[g])) for g in eachindex(model.tau)]
    return WidebandModel(admittance, propagation)
end

function UniversalLineModel.ULMModel(result::LineModel{<:WidebandModel})
    admittance = result.model.admittance
    propagation = result.model.propagation
    return ULMModel(size(admittance.direct, 1), result.line_length,
        copy(admittance.poles), [Matrix(slice) for slice in eachslice(admittance.residues; dims=3)],
        copy(admittance.direct), [group.delay for group in propagation],
        [copy(group.poles) for group in propagation],
        [[Matrix(slice) for slice in eachslice(group.residues; dims=3)] for group in propagation],
        nothing, Dict(:f=>copy(result.frequencies)))
end

function measure_fitting_error(model::ULMModel, data::LineFrequencyData)
    records = map((eval_yc, eval_h), (data.Yc, data.H)) do evaluate, samples
        fitted = [evaluate(model, s) for s in data.s]
        n = model.n
        rms = zeros(n, n)
        maxima = zeros(n, n)
        f_max = zeros(n, n)
        for j in 1:n, i in 1:n
            scale = maximum(abs(sample[i, j]) for sample in samples)
            errors = [abs(fit[i, j] - sample[i, j]) for (fit, sample) in zip(fitted, samples)]
            relative = iszero(scale) ? ifelse.(iszero.(errors), 0.0, Inf) : errors ./ scale
            rms[i, j] = sqrt(sum(abs2, relative) / length(relative))
            maxima[i, j], index = findmax(relative)
            f_max[i, j] = data.f[index]
        end
        return (rms=rms, max=maxima, f_max=f_max)
    end
    return (admittance=records[1], propagation=records[2])
end
