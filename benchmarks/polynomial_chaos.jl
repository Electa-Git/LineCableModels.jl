using DataFrames
using LineCableModels
using LinearAlgebra: norm
using PolyChaos
using Random
using Statistics

const QUANTITIES = (:R, :L, :C, :G)
const MONTE_CARLO_TRIALS = max(2, parse(Int, get(
    ENV, "LINECABLEMODELS_PCE_BENCHMARK_TRIALS", "5000")))
const SURROGATE_DRAWS = max(2, parse(Int, get(
    ENV, "LINECABLEMODELS_PCE_SURROGATE_DRAWS", "100000")))

conductor = Material(
    kind = :conductor,
    rho = 1.7241e-8,
    mu_r = Grid(0.999994, AbsoluteError(1.0e-4)),
    alpha = 0.00393
)
dielectric = Material(
    kind = :insulator,
    rho = 1.0e14,
    eps_r = Grid(2.3, AbsoluteError(0.01)),
    tan_delta = Grid(0.01, AbsoluteError(1.0e-4))
)
design = build(
    CableDesign,
    "polynomial-chaos-benchmark",
    terminal(:core, core(conductor; r = 0.01)),
    insulation(dielectric; t = 0.005)
)
problem = ParametricProblem(CableConstantsProblem(design; frequency = 50.0))

function moment_vector(product, field::Symbol)
    return reduce(vcat, (
        let summaries = getproperty(product, quantity)
            values = summaries isa NamedTuple ? getproperty(summaries, field) :
                     getproperty.(summaries, field)
            vec(values)
        end
        for quantity in QUANTITIES
    ))
end

function relative_error(candidate, reference)
    return norm(candidate - reference) / max(norm(reference), eps(Float64))
end

monte_carlo = nothing
monte_carlo_seconds = @elapsed global monte_carlo = compute(
    problem,
    MonteCarlo(
        CableConstantsFormulation(); trials = MONTE_CARLO_TRIALS, seed = 0x5eed)
)
reference = only(statistics(monte_carlo))
reference_mean = moment_vector(reference, :mean)
reference_std = moment_vector(reference, :std)
reference_R = only(reference.R)

rows = NamedTuple[(
    method = "Monte Carlo reference",
    core_solves = trial_count(monte_carlo, 1),
    wall_seconds = monte_carlo_seconds,
    mean_relative_error = 0.0,
    std_relative_error = 0.0,
    validation_error = missing,
    R_q05 = reference_R.q05,
    R_q50 = reference_R.median,
    R_q95 = reference_R.q95,
)]

for degree in 1:4
    completed = nothing
    seconds = @elapsed global completed = compute(
        problem,
        PolynomialChaos(
            CableConstantsFormulation();
            degree,
            validation_points = 32,
            validation_rtol = 1.0e-5,
            validation_seed = 0x5eed
        )
    )
    moments = only(statistics(completed))
    diagnostic = only(validation(completed))
    expansion = only(expansions(completed))

    # These are cheap offline surrogate evaluations, not core-solver trials
    # and therefore never enter samples(completed).
    dimension = diagnostic.stochastic_dimension
    latent = randn(
        Random.Xoshiro(0x5eed + UInt64(degree)), SURROGATE_DRAWS, dimension)
    R_surrogate = PolyChaos.evaluatePCE(
        vec(expansion.coefficients.R[1, :]), latent, expansion.basis)
    selected_quantiles = quantile(R_surrogate, [0.05, 0.5, 0.95])
    validation_error = maximum((
        Base.values(diagnostic.relative_rms_error)...,
        Base.values(diagnostic.relative_max_error)...,
    ))

    push!(rows, (
        method = "PolynomialChaos degree $degree",
        core_solves = diagnostic.collocation_evaluations +
                      diagnostic.validation_evaluations,
        wall_seconds = seconds,
        mean_relative_error = relative_error(
            moment_vector(moments, :mean), reference_mean),
        std_relative_error = relative_error(
            moment_vector(moments, :std), reference_std),
        validation_error,
        R_q05 = selected_quantiles[1],
        R_q50 = selected_quantiles[2],
        R_q95 = selected_quantiles[3],
    ))
end

display(DataFrame(rows))
