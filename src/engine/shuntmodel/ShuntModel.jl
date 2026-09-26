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
import ...Materials: Material
import ...DataModel: CableDesign
import ...DataModel
import ...TextDisplay
import ..Engine: CableBlueprint, InternalShuntBlock, internal_shunt_response,
                 blueprint_dependencies, same_physical_state, numerical_magnitude
import ..Engine: InsulationAdmittance, SemiconAdmittance
import ...LineCableModels: nominal
using LinearAlgebra: eigvals!, norm, I, lu, mul!, qr!, ColumnNorm, eigen,
                     SymTridiagonal, UpperTriangular
using QuadGK: alloc_segbuf, quadgk!
import SpecialFunctions

export Formula, formulas

include("interface.jl")
include("geometry.jl")
include("charge_collocation.jl")
include("blueprint.jl")

public BoundarySolveError
end
