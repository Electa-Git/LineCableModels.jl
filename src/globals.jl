# Export public API
export FormulationSet, DataFrame, add!, load!, export_data, save, preview
export f₀, μ₀, ε₀, ρ₀, T₀, TOL, ΔTmax
export setup_logging!
export BASE_FLOAT, REALSCALAR, COMPLEXSCALAR

# General constants
"Base power system frequency, f₀ = 50.0 [Hz]."
const f₀ = 50.0
"Magnetic constant (vacuum permeability), μ₀ = 4π * 1e-7 [H/m]."
const μ₀ = 4π * 1e-7
"Electric constant (vacuum permittivity), ε₀ = 8.8541878128e-12 [F/m]."
const ε₀ = 8.8541878128e-12
"Annealed copper reference resistivity, ρ₀ = 1.724e-08 [Ω·m]."
const ρ₀ = 1.724e-08
"Base temperature for conductor properties, T₀ = 20.0 [°C]."
const T₀ = 20.0
"Maximum tolerance for temperature variations, ΔTmax = 150 [°C]."
const ΔTmax = 150.0
"Default tolerance for floating-point comparisons, TOL = 1e-6."
const TOL = 1e-6

using Reexport, ForceImport

# Define aliases for the type constraints
using Measurements: Measurement
const BASE_FLOAT = Float64
const REALSCALAR = Union{BASE_FLOAT,Measurement{BASE_FLOAT}}
const COMPLEXSCALAR = Union{Complex{BASE_FLOAT},Complex{Measurement{BASE_FLOAT}}}

"""
    FormulationSet(...)

Constructs a specific formulation object based on the provided keyword arguments.
The system will infer the correct formulation type.
"""
FormulationSet(engine::Symbol; kwargs...) = FormulationSet(Val(engine); kwargs...)

function get_description end

function add! end

function load! end

function export_data end

function save end

function preview end


