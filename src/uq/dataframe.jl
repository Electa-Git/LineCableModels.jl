import DataFrames: DataFrame

function _confidence_columns(summary::SampleSummary, result)
    z = Distributions.quantile(
        Distributions.Normal(),
        0.5 + confidence(result) / 2
    )
    half_width = z * summary.std / sqrt(ntrials(result))
    relative = half_width / max(abs(summary.mean), eps(typeof(summary.mean)))
    return half_width, relative
end

function DataFrame(result::LineParametersMC)
    shape = size(result.statistics.R)
    frames = Array{DataFrame, 3}(undef, shape)
    for k in axes(frames, 3), j in axes(frames, 2), i in axes(frames, 1)
        summaries = (
            statistics(result, :R, i, j, k),
            statistics(result, :L, i, j, k),
            statistics(result, :C, i, j, k),
            statistics(result, :G, i, j, k)
        )
        confidence_values = _confidence_columns.(summaries, Ref(result))
        frames[i, j, k] = DataFrame(
            quantity = ["R", "L", "C", "G"],
            mean = getproperty.(summaries, :mean),
            std = getproperty.(summaries, :std),
            min = getproperty.(summaries, :min),
            q05 = getproperty.(summaries, :q05),
            q50 = getproperty.(summaries, :q50),
            q95 = getproperty.(summaries, :q95),
            max = getproperty.(summaries, :max),
            n = fill(ntrials(result), 4),
            conf = fill(confidence(result), 4),
            ci_half = first.(confidence_values),
            ci_rel = last.(confidence_values)
        )
    end
    return frames
end

function DataFrame(result::CableConstantsMC)
    summaries = (
        statistics(result, :R),
        statistics(result, :L),
        statistics(result, :C)
    )
    confidence_values = _confidence_columns.(summaries, Ref(result))
    return DataFrame(
        variable = ["R", "L", "C"],
        mean = getproperty.(summaries, :mean),
        std = getproperty.(summaries, :std),
        min = getproperty.(summaries, :min),
        q05 = getproperty.(summaries, :q05),
        q50 = getproperty.(summaries, :q50),
        q95 = getproperty.(summaries, :q95),
        max = getproperty.(summaries, :max),
        ntrials = fill(ntrials(result), 3),
        confidence = fill(confidence(result), 3),
        ci_half = first.(confidence_values),
        ci_rel = last.(confidence_values)
    )
end
