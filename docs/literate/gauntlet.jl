# # Benchmarks
#
# Gauntlet validates LineCableModels against external reference calculations
# and against independently implemented uncertainty methods. Every benchmark
# file records the compared scientific results, numerical differences, timing,
# software environment, and content hashes.
#
# ## External-reference benchmarks
#
# No PSCAD-specific result tensor survives past parsing.
#
# The complete recorded suite covers seven cable arrangements, each at 101
# frequencies from 1 Hz to 1 MHz. Relative values below the configured absolute
# floor are reported as `missing`; the absolute difference remains stored.

using LineCableModels #hide
using LineCableModels.Engine #hide
using CairoMakie #hide
using DataFrames #hide
using JLD2 #hide
using Measurements #hide

pscad_suite = DataFrame( #hide
    case = [ #hide
        "132 kV flat horizontal", #hide
        "18 kV trifoil", #hide
        "380 kV flat vertical", #hide
        "525 kV bipole", #hide
        "640 kV bipole", #hide
        "1000 mm² solid single", #hide
        "two bare wires", #hide
    ], #hide
    conductors = [9, 9, 9, 6, 6, 1, 2], #hide
    maximum_Z_rms_difference_percent = [ #hide
        2.1625, 1.7168, 1.8092, 1.9195, 1.9258, 1.1860, 14.9325, #hide
    ], #hide
    maximum_Y_rms_difference_percent = Union{Missing, Float64}[ #hide
        0.2469, 0.2469, 0.2469, 0.0005, missing, missing, missing, #hide
    ], #hide
    linecablemodels_median_ms = [ #hide
        27.828, 27.496, 27.171, 13.174, 13.362, 2.834, 5.356, #hide
    ], #hide
    pscad_wall_ms = [50.000, 74.000, 53.000, 71.000, 54.000, 73.000, 65.000], #hide
) #hide
pscad_suite

# These comparisons retain the method choices in each snapshot. The recorded
# PSCAD cases use Wedepohl earth return, while LineCableModels uses Pollaczek;
# the table therefore reports validation differences rather than bitwise
# equivalence. The two-bare-wire case has the largest relative Z difference.
# Its full absolute and relative matrices remain available in the snapshot.
#
# The code creates a small hidden snapshot with the fields used below. Set
# `snapshot_path` to a recorded case path when inspecting stored results.

frequencies_value = [1.0, 10.0, 100.0] #hide
reference_Z = zeros(ComplexF64, 2, 2, 3) #hide
reference_Y = zeros(ComplexF64, 2, 2, 3) #hide
for index in eachindex(frequencies_value) #hide
    scale = frequencies_value[index] #hide
    reference_Z[:, :, index] .= [0.1 + 0.2im 0.02 + 0.03im; 0.02 + 0.03im 0.11 + 0.21im] .* #hide
                                scale #hide
    reference_Y[:, :, index] .= [1.0e-6im 4.0e-20im; 2.0e-8im 1.2e-6im] .* scale #hide
end #hide
candidate_Z = reference_Z .* 1.01 #hide
candidate_Y = copy(reference_Y) #hide
candidate_Y[1, 2, :] .*= 1.5 #hide
candidate_Y[2, 2, :] .*= 1.00001 #hide
reference = LineParameters(PhaseDomain, reference_Z, reference_Y, frequencies_value) #hide
candidate = LineParameters(PhaseDomain, candidate_Z, candidate_Y, frequencies_value) #hide
fixture_comparison = compare(reference, candidate) #hide
fixture_benchmark = ( #hide
    minimum_seconds = 0.0129, #hide
    median_seconds = 0.0135, #hide
    bytes = 1_705_184, #hide
    allocations = 9_005, #hide
    samples = 10, #hide
    environment = ( #hide
        julia_version = string(VERSION), #hide
        kernel = string(Sys.KERNEL), #hide
        architecture = string(Sys.ARCH), #hide
        threads = Threads.nthreads(), #hide
        blas = "documentation fixture" #hide
    ) #hide
) #hide
fixture_directory = mktempdir() #hide
fixture_path = joinpath(fixture_directory, "snapshot.jld2") #hide
JLD2.jldsave( #hide
    fixture_path; #hide
    reference_comparison = fixture_comparison, #hide
    julia_benchmark = fixture_benchmark, #hide
    reference_execution = ( #hide
        backend = :fixture, #hide
        version = "1.0", #hide
        elapsed_seconds = 0.112 #hide
    ) #hide
) #hide

snapshot_path = joinpath(
    "test", "gauntlet", ".artifacts",
    "staging", "pscad", "benchmarks",
    "benchmark_525kv_1600mm2_bipole_pscad", "snapshot.jld2"
)
snapshot_path = fixture_path #hide
snapshot = JLD2.load(snapshot_path);

# ## Element-wise RMS errors
#
# `reference_comparison` compares the recorded external result against the
# accepted LineCableModels result. Each matrix entry is one RMS calculation
# over the complete frequency vector. The original matrices remain available
# without filtering.

comparison = snapshot["reference_comparison"]
comparison.Z.absolute

# The managed comparison table is wide: matrix coordinates identify each row,
# while absolute and relative Z and Y errors occupy their own columns. No
# display threshold suppresses the comparison values.

errors = DataFrame(comparison)
errors

# The modified Y mutual term appears directly in the corresponding wide-table
# column.

noise_row = only(eachrow(errors[(errors.row .== 1) .& (errors.column .== 2), :]))
(
    raw_relative = comparison.Y.relative[1, 2],
    displayed_relative = noise_row.εY,
)

# ## Impedance plots
#
# The comparison recipe accepts two or more `LineParameters` results and one
# legend entry for each result. Selecting `Z` produces separate real and
# imaginary impedance pages. Every matrix position keeps its own panel and the
# corresponding series are overlaid across frequency.

impedance_plots = CairoMakie.plot(
    reference,
    candidate;
    legend = ("Reference", "LineCableModels"),
    requests = (Z,),
    xscale = :log10,
    fig_size = (1000, 650),
    display_plot = false, #hide
    controls = false #hide
)

# The real part of the impedance matrix:
impedance_plots[1].figure #hide

# The imaginary part of the impedance matrix:
impedance_plots[2].figure #hide

# ## Execution measurements
#
# The snapshot keeps the accepted local benchmark and reference wall time. Timing
# gates apply only when the current and recorded environment identities match.

benchmark = snapshot["julia_benchmark"]
DataFrame(
    metric = [
        "Julia minimum", "Julia median", "allocated bytes",
        "allocations", "samples", "reference wall time"
    ],
    value = [
        1.0e3 * benchmark.minimum_seconds,
        1.0e3 * benchmark.median_seconds,
        benchmark.bytes,
        benchmark.allocations,
        benchmark.samples,
        snapshot["reference_execution"].elapsed_seconds
    ],
    unit = ["ms", "ms", "bytes", "count", "count", "s"]
)

# `benchmark.environment` contains the recorded Julia, operating-system, thread,
# and BLAS identifiers.

benchmark.environment

# ## Uncertainty quantification benchmarks
#
# Owned UQ benchmarks apply a 10% relative standard uncertainty to every
# continuous cable-layer geometry parameter and compare linear error
# propagation (LEP) with fixed-seed Monte Carlo. All seven recorded cases pass
# the engineering gates of 5% for means, 10% for standard deviations, and 2×
# for the Monte Carlo-to-LEP timing ratio.

uq_suite = DataFrame( #hide
    case = [ #hide
        "132 kV flat horizontal", #hide
        "18 kV trifoil", #hide
        "380 kV flat vertical", #hide
        "525 kV bipole", #hide
        "640 kV bipole", #hide
        "1000 mm² solid single", #hide
        "two bare wires", #hide
    ], #hide
    monte_carlo_trials = [512, 512, 512, 2048, 512, 512, 2048], #hide
    maximum_mean_difference_percent = [ #hide
        1.2295, 1.8823, 1.6797, 1.3101, 0.8628, 0.6833, 1.7548, #hide
    ], #hide
    maximum_std_difference_percent = [ #hide
        3.5812, 8.5139, 9.1603, 9.2480, 6.6836, 7.3302, 6.2435, #hide
    ], #hide
    monte_carlo_over_lep_time = [ #hide
        23.37, 6.25, 24.99, 80.76, 21.28, 34.28, 165.28, #hide
    ], #hide
) #hide
uq_suite

# Across the suite, the largest mean difference is 1.88% and the largest
# propagated standard-deviation difference is 9.25%. The methods therefore
# give consistent propagated uncertainties for the stated engineering gates.
# LEP is faster in every recorded case; the smallest observed timing ratio is
# 6.25×. Timing ratios are machine-specific.
#
# The plots below use the 18 kV, 1000 mm² trifoil case. It uses 512 accepted
# trials. Fixed wire counts and 20% formation clearance keep every accepted
# realisation within the supported geometry.

function _expand_uq_fixture(anchors) #hide
    return [index == 101 ? last(anchors) : begin #hide
                segment = div(index - 1, 10) + 1 #hide
                fraction = rem(index - 1, 10) / 10 #hide
                anchors[segment] + fraction * (anchors[segment + 1] - anchors[segment]) #hide
            end for index in 1:101] #hide
end #hide
function _uq_fixture_product(mean_anchors, std_anchors) #hide
    mean_values = _expand_uq_fixture(mean_anchors) #hide
    std_values = _expand_uq_fixture(std_anchors) #hide
    return (mean = reshape(mean_values, 1, 1, :), std = reshape(std_values, 1, 1, :)) #hide
end #hide
uq_fixture_frequencies = collect(10.0 .^ range(0, stop = 6, length = 101)) #hide
uq_fixture_zeros = zeros(11) #hide
uq_fixture_reference = ( #hide
    values = ( #hide
        R = _uq_fixture_product( #hide
            [3.567535312787254e-5, 3.863485615199304e-5, 5.06033512042166e-5, 1.0082637926629017e-4, 3.118557925966963e-4, 1.114080052889108e-3, 4.256781417674817e-3, 1.7093995507515593e-2, 7.018708008715781e-2, 2.8435793675690574e-1, 1.1765132611611668], #hide
            [6.9383330411604535e-6, 6.935326475949944e-6, 6.888462026428324e-6, 6.32268200539055e-6, 6.974233659660597e-6, 1.253562344924284e-5, 2.3859472114688702e-5, 7.788234170570126e-5, 1.8094904800021496e-4, 6.4400264661532e-4, 1.361271011375702e-3] #hide
        ), #hide
        L = _uq_fixture_product(uq_fixture_zeros, uq_fixture_zeros), #hide
        C = _uq_fixture_product( #hide
            fill(3.993576520461312e-10, 11), #hide
            fill(4.6452250314918605e-11, 11) #hide
        ), #hide
        G = _uq_fixture_product(uq_fixture_zeros, uq_fixture_zeros) #hide
    ), #hide
    frequencies = uq_fixture_frequencies, #hide
    basis = :pul, #hide
    domain = :PhaseDomain, #hide
    port_order = ["cable:1:core"] #hide
) #hide
uq_fixture_candidate = ( #hide
    values = ( #hide
        R = _uq_fixture_product( #hide
            [3.6463171459741566e-5, 3.942298920616265e-5, 5.13961983736761e-5, 1.0164292326951122e-4, 3.1228689010569656e-4, 1.114921051794952e-3, 4.258403698251689e-3, 1.7095719538089773e-2, 7.017683913341152e-2, 2.844724943615983e-1, 1.1766624247259188], #hide
            [7.707927333377946e-6, 7.704936435727461e-6, 7.658390905200914e-6, 7.1042357389742805e-6, 7.508153512661742e-6, 1.3438623738111746e-5, 2.5432295943057146e-5, 7.789377463360494e-5, 1.6949950245540816e-4, 6.489996351688275e-4, 1.3453048905541728e-3] #hide
        ), #hide
        L = _uq_fixture_product(uq_fixture_zeros, uq_fixture_zeros), #hide
        C = _uq_fixture_product( #hide
            fill(4.053770800971695e-10, 11), #hide
            fill(4.819929218973482e-11, 11) #hide
        ), #hide
        G = _uq_fixture_product(uq_fixture_zeros, uq_fixture_zeros) #hide
    ), #hide
    frequencies = uq_fixture_frequencies, #hide
    basis = :pul, #hide
    domain = :PhaseDomain, #hide
    port_order = ["cable:1:core"] #hide
) #hide
uq_fixture_directory = mktempdir() #hide
uq_fixture_path = joinpath(uq_fixture_directory, "snapshot.jld2") #hide
JLD2.jldsave( #hide
    uq_fixture_path; #hide
    accepted_reference = uq_fixture_reference, #hide
    accepted_candidate = uq_fixture_candidate #hide
) #hide

uq_snapshot_path = joinpath(
    "test", "gauntlet", ".artifacts", "staging", "uq", "benchmarks",
    "benchmark_18kv_1000mm2_trifoil_lep_montecarlo", "snapshot.jld2"
)
uq_snapshot_path = uq_fixture_path #hide
uq_snapshot = JLD2.load(uq_snapshot_path);

# The snapshot stores marginal moments without inventing another result type.
# For this comparison plot, an explicit presentation adapter combines each
# mean and standard deviation into a `Measurement` and reconstructs a
# `LineParameters` value. This does not claim Monte Carlo covariance: it only
# supplies the marginal uncertainty bars consumed by the existing line recipe.

function moment_line_parameters(record; row = 1, column = 1, samples = :)
    record.domain === :PhaseDomain || error("only phase-domain moments are supported")
    indices = samples isa Colon ? eachindex(record.frequencies) : samples
    moments = map(record.values) do product
        measurement.(
            product.mean[[row], [column], indices],
            product.std[[row], [column], indices]
        )
    end
    selected_frequencies = record.frequencies[indices]
    angular_frequency = reshape(2π .* selected_frequencies, 1, 1, :)
    return LineParameters(
        PhaseDomain,
        complex.(moments.R, moments.L .* angular_frequency),
        complex.(moments.G, moments.C .* angular_frequency),
        selected_frequencies;
        basis = record.basis
    )
end
nothing #hide

# Every tenth recorded frequency is enough to keep overlapping uncertainty
# whiskers legible in the documentation figure.

frequency_indices = 1:10:length(uq_snapshot["accepted_reference"].frequencies)
uq_plot_sources = (
    LEP = moment_line_parameters(
        uq_snapshot["accepted_reference"]; samples = frequency_indices),
    Monte_Carlo = moment_line_parameters(
        uq_snapshot["accepted_candidate"]; samples = frequency_indices)
)

resistance_uncertainty_plot = CairoMakie.plot(
    uq_plot_sources,
    @observe(R[1, 1, :]);
    legend = ("LEP", "Monte Carlo"),
    xscale = :log10,
    yscale = :log10,
    fig_size = (950, 420),
    display_plot = false, #hide
    controls = false #hide
)

capacitance_uncertainty_plot = CairoMakie.plot(
    uq_plot_sources,
    @observe(C[1, 1, :]);
    legend = ("LEP", "Monte Carlo"),
    xscale = :log10,
    fig_size = (950, 420),
    display_plot = false, #hide
    controls = false #hide
)

# Core self-resistance means and one-standard-deviation error bars:
resistance_uncertainty_plot[1].figure #hide

# Core self-capacitance means and one-standard-deviation error bars:
capacitance_uncertainty_plot[1].figure #hide

# ### Engineering agreement
#
# The recorded comparison calculates RMS differences over all 101 frequencies
# for every matrix term. Relative differences are reported only when the
# corresponding absolute RMS exceeds the quantity-specific numerical floor.

uq_summary = DataFrame(
    quantity = ["R", "L", "C", "G"],
    maximum_mean_difference_percent = Union{Missing, Float64}[
        0.0212, missing, 1.8823, missing
    ],
    maximum_std_difference_percent = Union{Missing, Float64}[
        4.8633, 4.6896, 8.5139, missing
    ]
)

# The largest meaningful mean difference is 1.88%, and the largest propagated
# standard-deviation difference is 8.51%. Both remain inside the benchmark's
# practical-engineering gates of 5% and 10%, respectively. `missing` denotes a
# comparison below its absolute numerical floor; the lossless formulation gives
# exactly zero shunt conductance in both calculations.

# ### Performance
#
# Timing is a separate warmed benchmark of the same calculations. These values
# are machine-specific and the snapshot records the complete environment.

DataFrame(
    method = ["LEP", "Monte Carlo"],
    execution_seconds = [6.407, 16.718],
    warmed_median_seconds = [2.992, 18.711],
    timing_samples = [3, 2]
)

# On the recorded machine, the warmed Monte Carlo median was 6.25× the LEP
# median, comfortably above the required 2× advantage. Monte Carlo completed
# two warmed timing samples within the 20-second sampling budget, so the
# ratio is an engineering performance indicator rather than a portable timing
# constant.

# ## Benchmark data API
#
# These result and observable types support owned line-parameter comparisons.
#
# ```@docs
# LineCableModels.Engine.RMSError
# LineCableModels.Engine.LineParametersBenchmark
# LineCableModels.Engine.compare
# LineCableModels.Engine.absolute_error
# LineCableModels.Engine.relative_error
# ```
