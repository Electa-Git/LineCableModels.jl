```@meta
EditURL = "../literate/gauntlet.jl"
```

# Inspecting PSCAD gauntlet results

A recorded gauntlet snapshot contains the external reference, the accepted
LineCableModels result, their element-wise comparison, and execution timing.
The snapshot uses the package's `LineParameters` and
`LineParametersComparison` types. No PSCAD-specific result tensor survives
past parsing.

This example creates a small hidden snapshot with the same fields as a
recorded case. Replace `snapshot_path` with the path of the case that you want
to inspect.

````@example gauntlet
using LineCableModels
using LineCableModels.Engine
using DataFrames
using JLD2

frequencies_value = [1.0, 10.0, 100.0] #hide
reference_Z = zeros(ComplexF64, 2, 2, 3) #hide
reference_Y = zeros(ComplexF64, 2, 2, 3) #hide
for index in eachindex(frequencies_value) #hide
    scale = frequencies_value[index] #hide
    reference_Z[:, :, index] .= [0.1 + 0.2im 0.02 + 0.03im; 0.02 + 0.03im 0.11 + 0.21im] .*
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
    pscad_elapsed_seconds = 0.112 #hide
) #hide

snapshot_path = joinpath(
    "test", "gauntlet", ".artifacts", "example3", "snapshot.jld2"
)
snapshot_path = fixture_path #hide
snapshot = JLD2.load(snapshot_path)
````

## Element-wise RMS errors

`reference_comparison` compares the recorded external result against the
accepted LineCableModels result. Each matrix entry is one RMS calculation
over the complete frequency vector. The original matrices remain available
without filtering.

````@example gauntlet
comparison = snapshot["reference_comparison"]
comparison.Z.absolute
````

A relative error is not useful when its absolute difference is already below
the chosen display floor. Pass the case's absolute Z and Y tolerances to
`DataFrame` as `zero_atol`. The table reports `missing` for those relative
entries, marks why they are unavailable, and leaves the raw comparison
untouched.

````@example gauntlet
errors = DataFrame(
    comparison;
    zero_atol = (Z = 1.0e-6, Y = 1.0e-9)
)
sort!(errors, [:quantity, :rms_absolute], rev = [false, true])
errors
````

The deliberately noisy Y mutual term still has its original relative value.
Only its table representation is suppressed.

````@example gauntlet
noise_row = only(eachrow(errors[(errors.quantity .== :Y) .& (errors.row .== 1) .& (errors.column .== 2), :]))
(
    raw_relative = comparison.Y.relative[1, 2],
    displayed_relative = noise_row.rms_relative,
    relative_status = noise_row.relative_status
)
````

## Execution measurements

The snapshot keeps the accepted local benchmark and PSCAD wall time. Timing
gates apply only when the current and recorded environment identities match.

````@example gauntlet
benchmark = snapshot["julia_benchmark"]
DataFrame(
    metric = [
        "Julia minimum", "Julia median", "allocated bytes",
        "allocations", "samples", "PSCAD wall time"
    ],
    value = [
        1.0e3 * benchmark.minimum_seconds,
        1.0e3 * benchmark.median_seconds,
        benchmark.bytes,
        benchmark.allocations,
        benchmark.samples,
        snapshot["pscad_elapsed_seconds"]
    ],
    unit = ["ms", "ms", "bytes", "count", "count", "s"]
)
````

The full recorded metadata can be inspected directly when a formulation or
execution environment needs closer review.

````@example gauntlet
benchmark.environment
````
