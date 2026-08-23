"""
$(TYPEDEF)

Abstract type for numerical problem definitions in
[`LineCableModels.jl`](index.md).
"""
abstract type AbstractImpedanceFormulation <: AbstractFormulation end
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

abstract type AbstractTransformFormulation <: AbstractFormulation end
abstract type AbstractEarthPropertiesFormulation <: AbstractFormulation end

abstract type _FormulationFamily end
struct _LineParameterFormulationFamily <: _FormulationFamily end

const _ACTIVE_FORMULATION_BACKENDS = Dict{DataType, Symbol}(
    _LineParameterFormulationFamily => :analytical,
)

function _active_formulation_backend(::Type{F}) where {F <: _FormulationFamily}
    return get(_ACTIVE_FORMULATION_BACKENDS, F) do
        throw(ArgumentError("no active backend is registered for formulation family $F"))
    end
end

function _activate_formulation_backend!(
        ::Type{F}, backend::Symbol
) where {F <: _FormulationFamily}
    _ACTIVE_FORMULATION_BACKENDS[F] = backend
    return backend
end

"""
$(TYPEDSIGNATURES)

Construct a formulation selected by `backend`, using the analytical formulation by
default.
"""
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)
function Formulation(; kwargs...)
    backend = _active_formulation_backend(_LineParameterFormulationFamily)
    return Formulation(Val(backend); kwargs...)
end
Formulation(::Val{:FEM}; kwargs...) = retired_fem_sector("FEM formulation")

"""
$(TYPEDEF)

Abstract type representing equivalent homogeneous earth models used by the
multi-dispatch earth-property evaluation.

# Currently available formulations

- [`EnforceLayer`](@ref): Effective parameters defined according to a specific earth layer.
"""
abstract type AbstractEHEMFormulation <: AbstractFormulation end
