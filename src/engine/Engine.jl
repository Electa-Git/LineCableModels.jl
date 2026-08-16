"""
    LineCableModels.Engine

The [`Engine`](@ref) module provides the main numerical functionality of
[`LineCableModels.jl`](index.md). It implements the materialized problems,
formulations, and primitive results used to calculate frequency-dependent
electrical parameters (Z/Y matrices) of line and cable systems.

# Overview

- Calculation of frequency-dependent series impedance (Z) and shunt admittance (Y) matrices.
- Direct uncertainty propagation through the optional Measurements extension.
- Internal impedance computation for solid, tubular and multi-layered coaxial conductors.
- Earth return impedances/admittances for overhead lines and underground cables (valid up to 10 MHz).
- Support for frequency-dependent soil properties.
- Handling of arbitrary polyphase systems with multiple conductors per phase.
- Phase and sequence domain calculations for ordinary and uncertain values.
- Novel N-layer concentric cable formulation with semiconductor modeling.

# Dependencies

$(IMPORTS)

"""
module Engine

# Export public API
export CableConstantsProblem, LineParametersProblem,
       LineParameters, SeriesImpedance, ShuntAdmittance,
       Z, Y, R, X, L, G, B, C,
       series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance,
       frequencies, nconductors, nfrequencies, basis,
       kronify
export AbstractFormulation, EMTFormulation, Formulation, EMTOptions, ComputeOptions

export compute!, plot

# Module-specific dependencies
using Reexport, ForceImport
using LinearAlgebra
using ..Commons
import ..Commons: get_description, LineParamsDomain, PhaseDomain, ModalDomain, domain
import ..Commons: retired_fem_sector
import ..Commons: basis, Z, Y, R, X, L, G, B, C,
                  series_impedance, shunt_admittance,
                  resistance, reactance, inductance,
                  conductance, susceptance, capacitance,
                  frequencies, nconductors, nfrequencies

using ..Utils
using ..UnitHandler
using ..PlotBuilder
using ..Materials
using ..EarthProps: EarthModel
using ..DataModel: CableDesign, CableConstants, LineCableSystem
import ..DataModel
using ..Utils: levelfrom, TimestampLogger
using Logging, LoggingExtras

include("types.jl")

# Problem definitions
include("lineparamopts.jl")
include("problemdefs.jl")
include("lineparams.jl")

# Submodule `InternalImpedance`
include("internalimpedance/InternalImpedance.jl")
using .InternalImpedance: InternalImpedance

# Submodule `InsulationImpedance`
include("insulationimpedance/InsulationImpedance.jl")
using .InsulationImpedance: InsulationImpedance

# Submodule `EarthImpedance`
include("earthimpedance/EarthImpedance.jl")
using .EarthImpedance: EarthImpedance

# Submodule `InsulationAdmittance`
include("insulationadmittance/InsulationAdmittance.jl")
using .InsulationAdmittance: InsulationAdmittance

# Submodule `EarthAdmittance`
include("earthadmittance/EarthAdmittance.jl")
using .EarthAdmittance: EarthAdmittance

# Submodule `Transforms`
include("transforms/Transforms.jl")
using .Transforms

# Submodule `EHEM`
include("ehem/EHEM.jl")
using .EHEM

# Helpers
include("helpers.jl")

# Workspace definition
include("workspace.jl")

# Computation methods
include("solver.jl")
include("reduction.jl")

# Override I/O methods
include("base.jl")
include("plotspecs.jl")
include("dataframe.jl")

"""
    plot(parameters[, quantities]; kwargs...)
    plot(impedance, frequencies[, quantities]; kwargs...)
    plot(admittance, frequencies[, quantities]; kwargs...)

Plot computed line parameters with a loaded Makie backend. `quantities` is a
tuple of accessors such as `(R, L, G, C)` or `(abs, angle)`.

Without an explicit selection, [`LineParameters`](@ref) produces separate
series-impedance and shunt-admittance figures. Each figure places the real part
on the left and the imaginary part on the right. Every selected matrix element
is represented by one data series.

Load `CairoMakie`, `GLMakie`, or `WGLMakie` before calling this function.

# Returns

- A `Vector{UIPlot}` containing one figure for each selected matrix family.
"""
function plot end

function plot(args...; kwargs...)
    throw(
        ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling plot.",
    ),
    )
end

# Removed FEM entry points
include("retired.jl")

@reexport using .InternalImpedance: InternalImpedance
@reexport using .InsulationImpedance: InsulationImpedance
@reexport using .EarthImpedance: EarthImpedance
@reexport using .InsulationAdmittance: InsulationAdmittance
@reexport using .EarthAdmittance: EarthAdmittance
@reexport using .EHEM, .Transforms

end # module Engine
