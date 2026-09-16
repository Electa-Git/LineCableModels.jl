"""
    LineCableModels.Engine.ShuntModel

Select the cable-local shunt geometry model independently of dielectric
constitutive laws. Coaxial annuli are the default; the boundary approximation
resolves eligible open wires and tapes inside a closed circular shield.

# Dependencies

$(IMPORTS)
"""
module ShuntModel
import ...Grammar: FormulationOptions
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: off
# Expanded in the module docstring, outside this module's analyzed expressions.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
import ..Engine: ShuntModelFormulation
import ...LineCableModels: FormulaDefinition, formula_id, description
import ...Grammar: formulation_options

export Formula, formulas

include("interface.jl")
end
