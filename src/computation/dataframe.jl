import DataFrames: DataFrame, metadata!

function _mc_unit(quantity::Symbol, result_basis, length_unit, quantity_units)
    return UnitHandler.line_component_unit(
        quantity,
        result_basis;
        length_unit,
        quantity_units
    )
end

function _mc_summary_frame(
        result::MonteCarloResult,
        quantities::Tuple,
        summaries::Tuple,
        result_basis::Symbol;
        length_unit::Symbol,
        quantity_units
)
    resolved_units = map(
        quantity -> _mc_unit(
            quantity,
            result_basis,
            length_unit,
            quantity_units
        ),
        quantities
    )
    scales = getproperty.(resolved_units, :scale)
    frame = DataFrame(
        quantity = collect(String.(quantities)),
        mean = collect(map((summary, scale) -> summary.mean * scale, summaries, scales)),
        std = collect(map((summary, scale) -> summary.std * abs(scale), summaries, scales)),
        min = collect(map((summary, scale) -> summary.min * scale, summaries, scales)),
        q05 = collect(map((summary, scale) -> summary.q05 * scale, summaries, scales)),
        q50 = collect(map((summary, scale) -> summary.q50 * scale, summaries, scales)),
        q95 = collect(map((summary, scale) -> summary.q95 * scale, summaries, scales)),
        max = collect(map((summary, scale) -> summary.max * scale, summaries, scales)),
        unit = collect(UnitHandler.get_label.(getproperty.(resolved_units, :units))),
        trials = fill(result.trials, length(quantities)),
        confidence = fill(result.confidence, length(quantities)),
        cdf_tol = fill(result.cdf_tol, length(quantities))
    )
    metadata!(
        frame,
        "monte_carlo",
        (
            trials = result.trials,
            confidence = result.confidence,
            cdf_tol = result.cdf_tol,
            distribution = string(result.distribution),
            seed = result.seed,
            manifest_hash = result.manifest.hash
        );
        style = :note
    )
    return frame
end

function _montecarlo_dataframe(
        result::MonteCarloResult,
        ::DataModel.CableConstants;
        length_unit::Symbol,
        quantity_units
)
    quantities = (:R, :L, :C)
    summaries = map(quantity -> getproperty(result.statistics, quantity), quantities)
    return _mc_summary_frame(
        result,
        quantities,
        summaries,
        :per_length;
        length_unit,
        quantity_units
    )
end

function _montecarlo_dataframe(
        result::MonteCarloResult,
        representation::Engine.LineParameters;
        length_unit::Symbol,
        quantity_units
)
    quantities = (:R, :L, :C, :G)
    shape = size(result.statistics.R)
    frames = Array{DataFrame, 3}(undef, shape)
    for index in CartesianIndices(frames)
        summaries = map(
            quantity -> getproperty(result.statistics, quantity)[index],
            quantities
        )
        frames[index] = _mc_summary_frame(
            result,
            quantities,
            summaries,
            basis(representation);
            length_unit,
            quantity_units
        )
    end
    return frames
end

function _montecarlo_dataframe(
        result::MonteCarloResult,
        representation;
        length_unit::Symbol,
        quantity_units
)
    throw(ArgumentError(
        "DataFrame presentation is not defined for MonteCarloResult{$(typeof(representation))}",
    ))
end

"""
$(TYPEDSIGNATURES)

Render the marginal summaries of a [`MonteCarloResult`](@ref) without
performing a calculation.

Cable-constant results produce one table with R, L, and C rows. Line-parameter
results produce an `n × n × n_frequency` array of tables with R, L, C, and G
rows. The DKW confidence and `cdf_tol` columns describe the simultaneous
marginal empirical-CDF bound used for automatic trial sizing; they are not
confidence intervals for the sample mean.

# Keywords

- `length_unit`: Metric prefix for the denominator of per-length quantities.
  Default: `:kilo`.
- `quantity_units`: Optional numerator-prefix override accepted by
  `UnitHandler.line_component_unit`.

# Returns

- A `DataFrame` for cable constants, or an array of `DataFrame`s for line
  parameters.
"""
function DataFrame(
        result::MonteCarloResult;
        length_unit::Symbol = :kilo,
        quantity_units = nothing
)
    return _montecarlo_dataframe(
        result,
        result.representation;
        length_unit,
        quantity_units
    )
end
