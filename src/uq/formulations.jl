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

Select variance-based sensitivity attribution for explicit native observable
requests. Dependency-owned method and sampler values determine the optional
integration that implements [`compute`](@ref).

$(TYPEDFIELDS)
"""
struct Sensitivity{
    F <: AbstractFormulation,
    M,
    R <: Tuple,
    S,
    L,
    O <: ComputationOptions
} <: AbstractFormulation
    "Formulation used for every concrete core problem."
    inner::F
    "Dependency-owned global-sensitivity method."
    method::M
    "Native observable requests in declaration order."
    requests::R
    "Base design sample count."
    samples::Int
    "Dependency-owned quasi-random sampler."
    sampler::S
    "Physical uncertainty mapping, `:normal` or `:uniform`."
    distribution::Symbol
    "Dependency-owned first-order estimator name."
    estimator::Symbol
    "Optional labels aligned with every uncertainty descriptor."
    input_labels::L
    "Maximum total number of core evaluations."
    max_evaluations::Int
    "Maximum flattened output values retained per outer point."
    max_output_values::Int
    "Supplemental-output retention options owned by this calculation."
    options::O
end

function computation_options(
        ::Type{Sensitivity},
        options::NamedTuple
)::ComputationOptions
    unknown = filter(key -> key !== :retain_details, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown Sensitivity computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((retain_details = false,), options)
    normalized.retain_details isa Bool || throw(ArgumentError(
        "Sensitivity retain_details must be Bool",
    ))
    return (retain_details = normalized.retain_details,)
end

function _validate_sensitivity_request(request)
    identity = request_identity(request)
    identity isa Function ||
        identity isa Tuple &&
        length(identity) == 2 &&
        all(value -> value isa Function, identity) ||
        throw(ArgumentError(
            "sensitivity requests must use observable selector functions",
        ))
    for index in request_indices(request)
        index isa Union{Integer, AbstractRange, AbstractVector, Colon} || throw(
            ArgumentError(
            "observable indices must be integers, ranges, vectors, or `:`",
        ),
        )
    end
    return nothing
end

function Sensitivity(
        inner::F,
        method,
        requests::Tuple;
        samples,
        sampler,
        distribution = :normal,
        estimator = :Jansen1999,
        input_labels = nothing,
        max_evaluations = 200_000,
        max_output_values = 20_000_000,
        options::NamedTuple = (;)
) where {F <: AbstractFormulation}
    isempty(requests) && throw(ArgumentError(
        "Sensitivity requires at least one observable request",
    ))
    foreach(_validate_sensitivity_request, requests)
    samples isa Integer && !(samples isa Bool) && samples >= 2 || throw(
        ArgumentError("Sensitivity samples must be an integer of at least two"),
    )
    distribution isa Symbol && distribution in (:normal, :uniform) || throw(
        ArgumentError("Sensitivity distribution must be :normal or :uniform"),
    )
    estimator isa Symbol || throw(ArgumentError(
        "Sensitivity estimator must be a Symbol",
    ))
    max_evaluations isa Integer && !(max_evaluations isa Bool) &&
    max_evaluations >= 1 || throw(ArgumentError(
        "Sensitivity max_evaluations must be a positive integer",
    ))
    max_output_values isa Integer && !(max_output_values isa Bool) &&
    max_output_values >= 1 || throw(ArgumentError(
        "Sensitivity max_output_values must be a positive integer",
    ))
    input_labels === nothing ||
        input_labels isa Tuple &&
        all(label -> label isa Union{AbstractString, Symbol}, input_labels) ||
        throw(
            ArgumentError("Sensitivity input_labels must be nothing or a tuple of strings or symbols"),
        )
    labels = input_labels === nothing ? nothing : map(String, input_labels)
    normalized = computation_options(Sensitivity, options)
    return Sensitivity{
        F,
        typeof(method),
        typeof(requests),
        typeof(sampler),
        typeof(labels),
        typeof(normalized)
    }(
        inner,
        method,
        requests,
        Int(samples),
        sampler,
        distribution,
        estimator,
        labels,
        Int(max_evaluations),
        Int(max_output_values),
        normalized
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
