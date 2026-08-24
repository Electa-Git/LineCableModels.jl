function _mc_unit(quantity::Symbol, result_basis, length_unit, quantity_units)
    return Units.line_component_unit(
        quantity,
        result_basis;
        length_unit,
        quantity_units
    )
end

function _mc_entry(result::MonteCarloResult)
    length(result) == 1 || throw(ArgumentError(
        "DataFrame requires one Monte Carlo point; select one result explicitly",
    ))
    observed = observables(result)
    return (
        representation = only(observed.result),
        statistics = only(observed.statistics),
        samples = samples(result) === nothing ? nothing : only(samples(result)),
        histograms = histograms(result) === nothing ? nothing : only(histograms(result)),
        trials = only(result.trial_counts),
        confidence = result.formulation.confidence,
        cdf_tol = result.formulation.cdf_tol,
        distribution = result.formulation.distribution,
        root_seed = result.root_seed,
        seed = only(result.point_seeds)
    )
end

function _mc_summary_frame(
        entry,
        quantities::Tuple,
        summaries::Tuple,
        result_basis::Symbol;
        length_unit::Symbol,
        quantity_units
)
    resolved_units = map(
        quantity -> _mc_unit(quantity, result_basis, length_unit, quantity_units),
        quantities
    )
    scales = getproperty.(resolved_units, :scale)
    frame = DataFrame(
        quantity = collect(String.(quantities)),
        mean = collect(map((summary, scale) -> summary.mean * scale, summaries, scales)),
        std = collect(map((summary, scale) -> summary.std * abs(scale), summaries, scales)),
        min = collect(map((summary, scale) -> summary.min * scale, summaries, scales)),
        q05 = collect(map((summary, scale) -> summary.q05 * scale, summaries, scales)),
        median = collect(map((summary, scale) -> summary.median * scale, summaries, scales)),
        q95 = collect(map((summary, scale) -> summary.q95 * scale, summaries, scales)),
        max = collect(map((summary, scale) -> summary.max * scale, summaries, scales)),
        n = getproperty.(summaries, :n),
        unit = collect(Units.label.(getproperty.(resolved_units, :units))),
        trials = fill(entry.trials, length(quantities)),
        confidence = fill(entry.confidence, length(quantities)),
        cdf_tol = fill(entry.cdf_tol, length(quantities))
    )
    metadata!(
        frame,
        "monte_carlo",
        (
            trials = entry.trials,
            confidence = entry.confidence,
            cdf_tol = entry.cdf_tol,
            distribution = string(entry.distribution),
            seed = entry.seed,
            root_seed = entry.root_seed
        );
        style = :note
    )
    return frame
end

function _montecarlo_dataframe(entry, ::DataModel.CableConstants;
        length_unit::Symbol, quantity_units)
    quantities = (:R, :L, :C)
    summaries = map(quantity -> getproperty(entry.statistics, quantity), quantities)
    return _mc_summary_frame(
        entry,
        quantities,
        summaries,
        :per_length;
        length_unit,
        quantity_units
    )
end

function _montecarlo_dataframe(entry, representation::Engine.LineParameters;
        length_unit::Symbol, quantity_units)
    quantities = (:R, :L, :C, :G)
    shape = size(entry.statistics.R)
    frames = Array{DataFrame, 3}(undef, shape)
    for index in CartesianIndices(frames)
        summaries = map(
            quantity -> getproperty(entry.statistics, quantity)[index],
            quantities
        )
        frames[index] = _mc_summary_frame(
            entry,
            quantities,
            summaries,
            basis(representation);
            length_unit,
            quantity_units
        )
    end
    return frames
end

function _montecarlo_dataframe(entry, representation; kwargs...)
    throw(ArgumentError(
        "DataFrame presentation is not defined for MonteCarloResult{$(typeof(representation))}",
    ))
end

function DataFrame(
        result::MonteCarloResult;
        length_unit::Symbol = :kilo,
        quantity_units = nothing
)
    entry = _mc_entry(result)
    return _montecarlo_dataframe(
        entry,
        entry.representation;
        length_unit,
        quantity_units
    )
end
