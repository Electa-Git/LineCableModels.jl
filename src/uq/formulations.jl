"""
$(TYPEDEF)

Select direct linear uncertainty propagation with `inner`.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation, O <: ComputationOptions} <: AbstractFormulation
    "Formulation used for each materialised problem."
    inner::F
    "Supplemental-output retention options owned by this propagation."
    options::O
end

function computation_options(
        ::Type{LinearError},
        options::NamedTuple
)::ComputationOptions
    unknown = filter(key -> key !== :retain_details, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LinearError computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((retain_details = false,), options)
    normalized.retain_details isa Bool || throw(ArgumentError(
        "LinearError retain_details must be Bool",
    ))
    return (retain_details = normalized.retain_details,)
end

function LinearError(
        inner::F;
        options::NamedTuple = (;)
) where {F <: AbstractFormulation}
    normalized = computation_options(LinearError, options)
    return LinearError{F, typeof(normalized)}(inner, normalized)
end

"""
$(TYPEDEF)

Select non-intrusive polynomial-chaos propagation with `inner`. The optional
PolyChaos extension owns basis construction, coefficient fitting, and
independent surrogate validation.

$(TYPEDFIELDS)
"""
struct PolynomialChaos{F, O <: ComputationOptions} <: AbstractFormulation
    "Formulation used for each collocation and validation problem."
    inner::F
    "Maximum total polynomial degree."
    degree::Int
    "Univariate Gaussian quadrature order."
    quadrature_order::Int
    "Independent standard-measure family."
    distribution::Symbol
    "Independent core solves used to validate each fitted expansion."
    validation_points::Int
    "Maximum relative RMS and maximum validation error."
    validation_rtol::Float64
    "Maximum collocation-plus-validation evaluations per outer point."
    max_evaluations::Int
    "Root seed used only for independent surrogate validation."
    validation_seed::UInt64
    "Supplemental-output retention options owned by this propagation."
    options::O
end

function computation_options(
        ::Type{PolynomialChaos},
        options::NamedTuple
)::ComputationOptions
    unknown = filter(key -> key !== :retain_details, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown PolynomialChaos computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((retain_details = false,), options)
    normalized.retain_details isa Bool || throw(ArgumentError(
        "PolynomialChaos retain_details must be Bool",
    ))
    return (retain_details = normalized.retain_details,)
end

function PolynomialChaos(
        inner;
        degree = 3,
        quadrature_order = nothing,
        distribution = :normal,
        validation_points = 16,
        validation_rtol = 1.0e-3,
        max_evaluations = 10_000,
        validation_seed = 0,
        options::NamedTuple = (;)
)
    degree isa Integer && !(degree isa Bool) && degree >= 1 || throw(
        ArgumentError("degree must be an integer greater than or equal to one"),
    )
    normalized_degree = Int(degree)
    order = quadrature_order === nothing ? normalized_degree + 1 : quadrature_order
    order isa Integer && !(order isa Bool) && order >= normalized_degree + 1 || throw(
        ArgumentError("quadrature_order must be at least degree + 1"),
    )
    distribution isa Symbol && distribution in (:normal, :uniform) || throw(
        ArgumentError("distribution must be :normal or :uniform"),
    )
    validation_points isa Integer && !(validation_points isa Bool) &&
        validation_points >= 1 || throw(ArgumentError(
        "validation_points must be a positive integer",
    ))
    validation_tolerance = Float64(validation_rtol)
    isfinite(validation_tolerance) && validation_tolerance > 0 || throw(
        ArgumentError("validation_rtol must be positive and finite"),
    )
    max_evaluations isa Integer && !(max_evaluations isa Bool) &&
        max_evaluations >= 1 || throw(ArgumentError(
        "max_evaluations must be a positive integer",
    ))
    seed = UInt64(validation_seed)
    normalized_options = computation_options(PolynomialChaos, options)
    return PolynomialChaos{typeof(inner), typeof(normalized_options)}(
        inner,
        normalized_degree,
        Int(order),
        distribution,
        Int(validation_points),
        validation_tolerance,
        Int(max_evaluations),
        seed,
        normalized_options
    )
end

"""
$(TYPEDEF)

Select conditional Monte Carlo propagation over a
[`ParametricProblem`](@ref). Randomness is local and reproducible when `seed`
is supplied. Computation option `on_error=:fail` propagates every exception.
`on_error=:retry` rejects only realisations that raise `DomainError`, retains
their sampled arguments and error summaries, and continues until the requested
number of successful trials is obtained or `max_failures` is reached. Retry
mode requires `retain_details=true` and estimates the output distribution
conditional on successful problem construction and computation.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, D, S, O <: ComputationOptions} <:
       AbstractFormulation
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
    "Whether joint samples are retained."
    return_samples::Bool
    "Whether marginal histogram densities are retained."
    return_histograms::Bool
    "Optional histogram bin count."
    bins::Union{Nothing, Int}
    "Supplemental-output retention options owned by this propagation."
    options::O
    function MonteCarlo(
            inner::F;
            trials::Union{Nothing, Integer} = nothing,
            confidence::Real = 0.95,
            cdf_tol::Real = 0.02,
            distribution = :normal,
            seed::Union{Nothing, Integer} = nothing,
            return_samples::Bool = false,
            return_histograms::Bool = false,
            bins::Union{Nothing, Integer} = nothing,
            options::NamedTuple = (;)
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
        normalized_options = computation_options(MonteCarlo, options)
        return new{F, typeof(distribution), typeof(actual_seed), typeof(normalized_options)}(
            inner,
            trials === nothing ? nothing : Int(trials),
            Float64(confidence),
            Float64(cdf_tol),
            distribution,
            actual_seed,
            return_samples,
            return_histograms,
            bins === nothing ? nothing : Int(bins),
            normalized_options
        )
    end
end

function computation_options(
        ::Type{MonteCarlo},
        options::NamedTuple
)::ComputationOptions
    supported = (:retain_details, :on_error, :max_failures)
    unknown = filter(key -> key ∉ supported, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown MonteCarlo computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(
        (retain_details = false, on_error = :fail, max_failures = 100),
        options
    )
    normalized.retain_details isa Bool || throw(ArgumentError(
        "MonteCarlo retain_details must be Bool",
    ))
    normalized.on_error in (:fail, :retry) || throw(ArgumentError(
        "MonteCarlo on_error must be :fail or :retry",
    ))
    normalized.max_failures isa Integer && !(normalized.max_failures isa Bool) &&
    normalized.max_failures > 0 || throw(ArgumentError(
        "MonteCarlo max_failures must be a positive integer",
    ))
    normalized.on_error === :retry && !normalized.retain_details &&
        throw(
            ArgumentError(
            "MonteCarlo on_error=:retry requires retain_details=true",
        ),
        )
    return (
        retain_details = normalized.retain_details,
        on_error = normalized.on_error,
        max_failures = Int(normalized.max_failures)
    )
end
