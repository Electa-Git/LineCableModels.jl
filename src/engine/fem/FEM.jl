"""
    LineCableModels.Engine.FEM

Provide the deprecated compatibility facade for the optional finite-element
integration. Load `Gmsh` before constructing or running an FEM formulation.
"""
module FEM

export Darwin, Electrodynamics, FormulationSet, MeshTransition, calc_domain_size,
       compute!, preview_results

import ...Engine: FormulationSet, compute!

const _DEPRECATION_MESSAGE = "The current FEM interface is deprecated and will be replaced by a simplified API in a future release."

_parent_package() = parentmodule(parentmodule(@__MODULE__))

function _extension()
    ext = Base.get_extension(_parent_package(), :LineCableModelsGmshExt)
    ext === nothing && throw(
        ArgumentError(
        "FEM is optional. Load Gmsh with `using Gmsh` before using LineCableModels.Engine.FEM.",
    ),
    )
    return ext
end

function _warn_fem(symbol::Symbol)
    Base.depwarn(_DEPRECATION_MESSAGE, symbol)
    return nothing
end

function Darwin(args...; kwargs...)
    _warn_fem(:Darwin)
    return _extension().Darwin(args...; kwargs...)
end

function Electrodynamics(args...; kwargs...)
    _warn_fem(:Electrodynamics)
    return _extension().Electrodynamics(args...; kwargs...)
end

function MeshTransition(args...; kwargs...)
    _warn_fem(:MeshTransition)
    return _extension().MeshTransition(args...; kwargs...)
end

"""
    calc_domain_size(earth_model, frequencies; min_radius=5.0, max_radius=5000.0)

Estimate the radial FEM domain size from the most resistive finite earth layer.

# Arguments

- `earth_model`: Earth model containing frequency-dependent resistivity and
  permeability.
- `frequencies`: FEM frequencies in hertz.

# Returns

- Domain radius in metres, clamped to `min_radius:max_radius`.

# Notes

The implementation evaluates the magnitude of the complex skin depth at the
first frequency,

```math
\\delta = \\left|\\sqrt{\\frac{\\rho}{\\mathrm{j}\\,2\\pi f\\mu}}\\right|,
```

using the finite earth layer with the largest resistivity.
"""
function calc_domain_size(args...; kwargs...)
    _warn_fem(:calc_domain_size)
    return _extension().calc_domain_size(args...; kwargs...)
end

function preview_results(args...; kwargs...)
    _warn_fem(:preview_results)
    return _extension().preview_results(args...; kwargs...)
end

function FormulationSet(::Val{:FEM}; kwargs...)
    _warn_fem(:FormulationSet)
    return _extension().formulation_set(; kwargs...)
end

end # module FEM
