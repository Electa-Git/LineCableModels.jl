"""
$(TYPEDEF)

Select direct linear uncertainty propagation with `inner`.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation} <: AbstractFormulation
    "Formulation used for each materialised problem."
    inner::F
end

"""
$(TYPEDEF)

Select conditional Monte Carlo propagation over a
[`ParametricProblem`](@ref). Randomness is local and reproducible when `seed`
is supplied.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, D, S} <: AbstractFormulation
    "Formulation used for each sampled problem."
    inner::F
    "Requested trials, or `nothing` for DKW sizing."
    trials::Union{Nothing, Int}
    "Simultaneous empirical-CDF confidence [dimensionless]."
    confidence::Float64
    "Maximum empirical-CDF deviation used for DKW sizing [dimensionless]."
    cdf_tol::Float64
    "Sampling family or extension-provided univariate distribution."
    distribution::D
    "Optional root random seed."
    seed::S
    "Whether joint primitive samples are retained."
    return_samples::Bool
    "Whether marginal histogram densities are retained."
    return_histograms::Bool
    "Optional histogram bin count."
    bins::Union{Nothing, Int}
    function MonteCarlo(
            inner::F;
            trials::Union{Nothing, Integer} = nothing,
            confidence::Real = 0.95,
            cdf_tol::Real = 0.02,
            distribution = :normal,
            seed::Union{Nothing, Integer} = nothing,
            return_samples::Bool = false,
            return_histograms::Bool = false,
            bins::Union{Nothing, Integer} = nothing
    ) where {F <: AbstractFormulation}
        trials === nothing || trials > 0 || throw(ArgumentError("trials must be positive"))
        0 < confidence < 1 ||
            throw(ArgumentError("confidence must lie between zero and one"))
        0 < cdf_tol < 1 || throw(ArgumentError("cdf_tol must lie between zero and one"))
        bins === nothing || bins > 0 || throw(ArgumentError("bins must be positive"))
        distribution isa Symbol && distribution ∉ (:normal, :uniform) &&
            throw(ArgumentError(
                "unsupported distribution $(repr(distribution)); expected :normal, :uniform, a sampler function, or an extension-supported distribution",
            ))
        actual_seed = seed === nothing ? nothing : UInt64(seed)
        return new{F, typeof(distribution), typeof(actual_seed)}(
            inner,
            trials === nothing ? nothing : Int(trials),
            Float64(confidence),
            Float64(cdf_tol),
            distribution,
            actual_seed,
            return_samples,
            return_histograms,
            bins === nothing ? nothing : Int(bins)
        )
    end
end
