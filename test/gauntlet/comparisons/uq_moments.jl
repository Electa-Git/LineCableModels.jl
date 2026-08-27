struct MomentResult{V <: NamedTuple, F <: AbstractVector, D}
    values::V
    frequencies::F
    basis::Symbol
    domain::D
    port_order::Vector{String}
end

struct MomentBenchmark{E <: NamedTuple}
    errors::E
end

const _MOMENT_QUANTITIES = (
    R = LineCableModels.R,
    L = LineCableModels.L,
    C = LineCableModels.C,
    G = LineCableModels.G
)

function _moment_metadata(parameters, port_order)
    return (
        frequencies = Float64.(
            LineCableModels.nominal.(LineCableModels.frequencies(parameters))
        ),
        basis = LineCableModels.basis(parameters),
        domain = LineCableModels.domain(parameters),
        port_order = String.(port_order)
    )
end

function extract_moments(
        result::LineCableModels.LinearErrorResult,
        port_order::AbstractVector{<:AbstractString}
)
    length(result) == 1 || throw(DimensionMismatch(
        "UQ Gauntlet moments require exactly one Gridspace point",
    ))
    parameters = only(result)
    products = map(_MOMENT_QUANTITIES) do selector
        uncertain = LineCableModels.observe(parameters, selector)
        (
            mean = Float64.(LineCableModels.nominal.(uncertain)),
            std = Float64.(LineCableModels.standard_uncertainty.(uncertain))
        )
    end
    metadata = _moment_metadata(parameters, port_order)
    return MomentResult(
        products,
        metadata.frequencies,
        metadata.basis,
        metadata.domain,
        metadata.port_order
    )
end

function extract_moments(
        result::LineCableModels.MonteCarloResult,
        port_order::AbstractVector{<:AbstractString}
)
    length(result) == 1 || throw(DimensionMismatch(
        "UQ Gauntlet moments require exactly one Gridspace point",
    ))
    parameters = only(result)
    summaries = only(LineCableModels.statistics(result))
    keys(summaries) == keys(_MOMENT_QUANTITIES) || throw(ArgumentError(
        "Monte Carlo summaries must contain exactly R, L, C, and G",
    ))
    products = map(summaries) do field
        (
            mean = Float64.(getproperty.(field, :mean)),
            std = Float64.(getproperty.(field, :std))
        )
    end
    metadata = _moment_metadata(parameters, port_order)
    return MomentResult(
        products,
        metadata.frequencies,
        metadata.basis,
        metadata.domain,
        metadata.port_order
    )
end

function _validate_moment_result(value::MomentResult)
    keys(value.values) == keys(_MOMENT_QUANTITIES) || throw(ArgumentError(
        "UQ moment products must contain exactly R, L, C, and G",
    ))
    isempty(value.frequencies) && throw(ArgumentError(
        "UQ moment frequency axis cannot be empty",
    ))
    terminal_count = length(value.port_order)
    terminal_count > 0 || throw(ArgumentError(
        "UQ moment terminal order cannot be empty",
    ))
    for quantity in keys(value.values)
        product = getproperty(value.values, quantity)
        keys(product) == (:mean, :std) || throw(ArgumentError(
            "UQ $quantity moment product must contain mean and std",
        ))
        expected = (terminal_count, terminal_count, length(value.frequencies))
        size(product.mean) == expected || throw(DimensionMismatch(
            "UQ $quantity mean has size $(size(product.mean)); expected $expected",
        ))
        size(product.std) == expected || throw(DimensionMismatch(
            "UQ $quantity std has size $(size(product.std)); expected $expected",
        ))
        all(isfinite, product.mean) || throw(ArgumentError(
            "UQ $quantity means must be finite",
        ))
        all(value -> isfinite(value) && value >= 0, product.std) ||
            throw(ArgumentError("UQ $quantity standard deviations must be finite and nonnegative"))
    end
    return value
end

function _moment_rms(reference::AbstractArray{<:Real, 3}, candidate::AbstractArray{<:Real, 3})
    errors = [begin
        left = @view reference[row, column, :]
        right = @view candidate[row, column, :]
        difference_norm = sum(abs2, left .- right)
        reference_norm = sum(abs2, left)
        absolute = sqrt(difference_norm / length(left))
        relative = iszero(reference_norm) ?
                   (iszero(difference_norm) ? zero(absolute) : oftype(absolute, Inf)) :
                   sqrt(difference_norm / reference_norm)
        (; absolute, relative)
    end for row in axes(reference, 1), column in axes(reference, 2)]
    absolute = getproperty.(errors, :absolute)
    relative = getproperty.(errors, :relative)
    T = promote_type(eltype(absolute), eltype(relative))
    return RMSError{T}(Matrix{T}(absolute), Matrix{T}(relative))
end

function Base.isapprox(reference::MomentResult, candidate::MomentResult; kwargs...)
    reference.frequencies == candidate.frequencies || return false
    reference.basis === candidate.basis || return false
    reference.domain === candidate.domain || return false
    reference.port_order == candidate.port_order || return false
    return all(keys(reference.values)) do quantity
        left = getproperty(reference.values, quantity)
        right = getproperty(candidate.values, quantity)
        isapprox(left.mean, right.mean; kwargs...) &&
            isapprox(left.std, right.std; kwargs...)
    end
end

function LineCableModels.Engine.compare(
        reference::MomentResult,
        candidate::MomentResult
)
    _validate_moment_result(reference)
    _validate_moment_result(candidate)
    reference.frequencies == candidate.frequencies || throw(ArgumentError(
        "UQ moment frequencies must match exactly and in order",
    ))
    reference.basis === candidate.basis || throw(ArgumentError(
        "UQ moment bases must match",
    ))
    reference.domain === candidate.domain || throw(ArgumentError(
        "UQ moment domains must match",
    ))
    reference.port_order == candidate.port_order || throw(ArgumentError(
        "UQ moment terminal orders must match",
    ))
    errors = map(reference.values, candidate.values) do left, right
        size(left.mean) == size(right.mean) || throw(DimensionMismatch(
            "UQ moment mean dimensions must match",
        ))
        size(left.std) == size(right.std) || throw(DimensionMismatch(
            "UQ moment standard-deviation dimensions must match",
        ))
        (
            mean = _moment_rms(left.mean, right.mean),
            std = _moment_rms(left.std, right.std)
        )
    end
    return MomentBenchmark(errors)
end

function moment_comparison_passes(comparison::MomentBenchmark, tolerances::NamedTuple)
    keys(tolerances) == (:mean, :std) || throw(ArgumentError(
        "UQ moment tolerances must contain exactly mean and std",
    ))
    for statistic in (:mean, :std), quantity in keys(_MOMENT_QUANTITIES)
        error = getproperty(getproperty(comparison.errors, quantity), statistic)
        statistic_tolerances = getproperty(tolerances, statistic)
        haskey(statistic_tolerances, quantity) || throw(ArgumentError(
            "UQ $statistic tolerances are missing $quantity",
        ))
        limit = getproperty(statistic_tolerances, quantity)
        all((error.absolute .<= limit.absolute) .|
            (error.relative .<= limit.relative)) || return false
    end
    return true
end

function moment_error_summary(comparison::MomentBenchmark)
    return map(comparison.errors) do product
        map(product) do error
            absolute_index = argmax(error.absolute)
            relative_index = argmax(error.relative)
            (
                maximum_absolute = error.absolute[absolute_index],
                absolute_location = Tuple(absolute_index),
                maximum_relative = error.relative[relative_index],
                relative_location = Tuple(relative_index)
            )
        end
    end
end


function moment_error_summary(
        comparison::MomentBenchmark,
        tolerances::NamedTuple
)
    mean_limits = tolerances.mean
    std_limits = tolerances.std
    return map(comparison.errors, mean_limits, std_limits) do product, mean_limit, std_limit
        limits = (mean = mean_limit, std = std_limit)
        map(product, limits) do error, limit
            gating_ratio = min.(
                error.absolute ./ limit.absolute,
                error.relative ./ limit.relative
            )
            index = argmax(gating_ratio)
            visible = findall(error.absolute .> limit.absolute)
            relative = isempty(visible) ? missing : maximum(error.relative[visible])
            (
                maximum_absolute = maximum(error.absolute),
                maximum_relative_above_absolute_floor = relative,
                maximum_gating_ratio = gating_ratio[index],
                gating_location = Tuple(index)
            )
        end
    end
end
