"""
	LineCableModels.Utils

The [`Utils`](@ref) module provides utility functions and constants for the  [`LineCableModels.jl`](index.md) package. This module includes functions for handling measurements, numerical comparisons, and other common tasks.

# Overview

- Provides general constants used throughout the package.
- Includes utility functions for numerical comparisons and handling measurements.
- Contains functions to compute uncertainties and bounds for measurements.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Utils

# Export public API
export resolve_T, coerce_to_T, _CLEANMETHODLIST, resolve_backend, is_headless

export to_nominal,
    to_certain,
    percent_to_uncertain,
    bias_to_uncertain,
    to_upper,
    to_lower,
    percent_error

# Load common dependencies
using ..LineCableModels
include("utils/commondeps.jl")

# Module-specific dependencies
using Measurements: Measurement, value, uncertainty, measurement
using Statistics
using Plots

"""
$(TYPEDSIGNATURES)

Extracts the nominal value from a measurement or returns the original value.

# Arguments

- `x`: Input value which can be a `Measurement` type or any other type.

# Returns

- The nominal value if `x` is a `Measurement`, otherwise returns `x` unchanged.

# Examples

```julia
using Measurements

$(FUNCTIONNAME)(1.0)  # Output: 1.0
$(FUNCTIONNAME)(5.2 ± 0.3)  # Output: 5.2
```
"""
function to_nominal(x)
    return x isa Measurement ? Measurements.value(x) : x
end

"""
$(TYPEDSIGNATURES)

Converts a measurement to a value with zero uncertainty, retaining the numeric type `Measurement`.

# Arguments

- `value`: Input value that may be a `Measurement` type or another type.

# Returns

- If input is a `Measurement`, returns the same value with zero uncertainty; otherwise returns the original value unchanged.

# Examples

```julia
x = 5.0 ± 0.1
result = $(FUNCTIONNAME)(x)  # Output: 5.0 ± 0.0

y = 10.0
result = $(FUNCTIONNAME)(y)  # Output: 10.0
```
"""
function to_certain(value)
    return value isa Measurement ? (Measurements.value(value) ± 0.0) : value
end

"""
$(TYPEDSIGNATURES)

Converts a value to a measurement with uncertainty based on percentage.

# Arguments

- `val`: The nominal value.
- `perc`: The percentage uncertainty (0 to 100).

# Returns

- A `Measurement` type with the given value and calculated uncertainty.

# Examples

```julia
using Measurements

$(FUNCTIONNAME)(100.0, 5)  # Output: 100.0 ± 5.0
$(FUNCTIONNAME)(10.0, 10)  # Output: 10.0 ± 1.0
```
"""
function percent_to_uncertain(val, perc) #perc from 0 to 100
    measurement(val, (perc * val) / 100)
end

"""
$(TYPEDSIGNATURES)

Computes the uncertainty of a measurement by incorporating systematic bias.

# Arguments

- `nominal`: The deterministic nominal value (Float64).
- `measurements`: A vector of `Measurement` values from the `Measurements.jl` package.

# Returns

- A new `Measurement` object representing the mean measurement value with an uncertainty that accounts for both statistical variation and systematic bias.

# Notes

- Computes the mean value and its associated uncertainty from the given `measurements`.
- Determines the **bias** as the absolute difference between the deterministic `nominal` value and the mean measurement.
- The final uncertainty is the sum of the standard uncertainty (`sigma_mean`) and the systematic bias.

# Examples
```julia
using Measurements

nominal = 10.0
measurements = [10.2 ± 0.1, 9.8 ± 0.2, 10.1 ± 0.15]
result = $(FUNCTIONNAME)(nominal, measurements)
println(result)  # Output: Measurement with adjusted uncertainty
```
"""
function bias_to_uncertain(nominal::Float64, measurements::Vector{<:Measurement})
    # Compute the mean value and uncertainty from the measurements
    mean_measurement = mean(measurements)
    mean_value = Measurements.value(mean_measurement)  # Central value
    sigma_mean = Measurements.uncertainty(mean_measurement)  # Uncertainty of the mean
    # Compute the bias (deterministic nominal value minus mean measurement)
    bias = abs(nominal - mean_value)
    return mean_value ± (sigma_mean + bias)
end

"""
$(TYPEDSIGNATURES)

Computes the upper bound of a measurement value.

# Arguments

- `m`: A numerical value, expected to be of type `Measurement` from the `Measurements.jl` package.

# Returns

- The upper bound of `m`, computed as `value(m) + uncertainty(m)` if `m` is a `Measurement`.
- `NaN` if `m` is not a `Measurement`.

# Examples

```julia
using Measurements

m = 10.0 ± 2.0
upper = $(FUNCTIONNAME)(m)  # Output: 12.0

not_a_measurement = 5.0
upper_invalid = $(FUNCTIONNAME)(not_a_measurement)  # Output: NaN
```
"""
function to_upper(m::Number)
    if m isa Measurement
        return Measurements.value(m) + Measurements.uncertainty(m)
    else
        return NaN
    end
end

"""
$(TYPEDSIGNATURES)

Computes the lower bound of a measurement value.

# Arguments

- `m`: A numerical value, expected to be of type `Measurement` from the `Measurements.jl` package.

# Returns

- The lower bound, computed as `value(m) - uncertainty(m)` if `m` is a `Measurement`.
- `NaN` if `m` is not a `Measurement`.

# Examples

```julia
using Measurements

m = 10.0 ± 2.0
lower = $(FUNCTIONNAME)(m)  # Output: 8.0

not_a_measurement = 5.0
lower_invalid = $(FUNCTIONNAME)(not_a_measurement)  # Output: NaN
```
"""
function to_lower(m::Number)
    if m isa Measurement
        return Measurements.value(m) - Measurements.uncertainty(m)
    else
        return NaN
    end
end

"""
$(TYPEDSIGNATURES)

Computes the percentage uncertainty of a measurement.

# Arguments

- `m`: A numerical value, expected to be of type `Measurement` from the `Measurements.jl` package.

# Returns

- The percentage uncertainty, computed as `100 * uncertainty(m) / value(m)`, if `m` is a `Measurement`.
- `NaN` if `m` is not a `Measurement`.

# Examples

```julia
using Measurements

m = 10.0 ± 2.0
percent_err = $(FUNCTIONNAME)(m)  # Output: 20.0

not_a_measurement = 5.0
percent_err_invalid = $(FUNCTIONNAME)(not_a_measurement)  # Output: NaN
```
"""
function percent_error(m::Number)
    if m isa Measurement
        return 100 * Measurements.uncertainty(m) / Measurements.value(m)
    else
        return NaN
    end
end

@inline _nudge_float(x::AbstractFloat) = isfinite(x) && x == trunc(x) ? nextfloat(x) : x #redundant and I dont care

_coerce_args_to_T(args...) =
    any(x -> x isa Measurement, args) ? Measurement{BASE_FLOAT} : BASE_FLOAT

# Promote scalar to T if T is Measurement; otherwise take nominal if x is Measurement.
function _coerce_scalar_to_T(x, ::Type{T}) where {T}
    if T <: Measurement
        return x isa Measurement ? x : (zero(T) + x)
    else
        return x isa Measurement ? T(value(x)) : convert(T, x)
    end
end

# Arrays: promote/demote elementwise, preserving shape. Arrays NEVER decide T.
function _coerce_array_to_T(A::AbstractArray, ::Type{T}) where {T}
    if T <: Measurement
        return (eltype(A) === T) ? A : (A .+ zero(T))             # Real → Measurement(σ=0)
    elseif eltype(A) <: Measurement
        B = value.(A)                                             # Measurement → Real (nominal)
        return (eltype(B) === T) ? B : convert.(T, B)
    else
        return (eltype(A) === T) ? A : convert.(T, A)
    end
end

"""
$(TYPEDSIGNATURES)

Determines if the current execution environment is headless (without display capability).

# Returns

- `true` if running in a continuous integration environment or without display access.
- `false` otherwise when a display is available.

# Examples

```julia
if $(FUNCTIONNAME)()
	# Use non-graphical backend
	gr()
else
	# Use interactive backend
	plotlyjs()
end
```
"""
function is_headless()::Bool
    # 1. Check for common CI environment variables
    if get(ENV, "CI", "false") == "true"
        return true
    end

    # 2. Check if a display is available (primarily for Linux)
    if !haskey(ENV, "DISPLAY") && Sys.islinux()
        return true
    end

    # 3. Check for GR backend's specific headless setting
    if get(ENV, "GKSwstype", "") in ("100", "nul", "nil")
        return true
    end

    return false
end

function display_path(file_name)
    return is_headless() ? basename(file_name) : relpath(file_name)
end

"""
$(TYPEDSIGNATURES)

Checks if the code is running inside a `@testset` by checking if `Test` is loaded
in the current session and then calling `get_testset_depth()`.
"""
function is_in_testset()
    # Start with the current module
    current_module = @__MODULE__

    # Walk up the module tree (e.g., from the sandbox to Main)
    while true
        if isdefined(current_module, :Test) &&
           isdefined(current_module.Test, :get_testset_depth)
            # Found the Test module, check the test set depth
            return current_module.Test.get_testset_depth() > 0
        end

        # Move to the parent module
        parent = parentmodule(current_module)
        if parent === current_module # Reached the top (Main)
            break
        end
        current_module = parent
    end

    return false
end

"""
$(TYPEDSIGNATURES)

Selects the appropriate plotting backend based on the environment.

# Arguments

- `backend`: Optional explicit backend to use. If provided, this backend will be activated.

# Returns

Nothing. The function activates the chosen backend.

# Notes

Automatically selects GR for headless environments (CI or no DISPLAY) and PlotlyJS
for interactive use when no backend is explicitly specified. This is particularly needed when running within CI environments.

# Examples

```julia
resolve_backend()           # Auto-selects based on environment
resolve_backend(pyplot)     # Explicitly use PyPlot backend
```
"""
function resolve_backend(backend=nothing)
    if isnothing(backend) # Check if running in a headless environment 
        if is_headless() # Use GR for CI/headless environments
            ENV["GKSwstype"] = "100"
            gr()
        else # Use PlotlyJS for interactive use 
            plotlyjs()
        end
    else # Use the specified backend if provided 
        backend()
    end
end

include("utils/docstringextension.jl")
include("utils/typecoercion.jl")
include("utils/macros.jl")

end # module Utils
