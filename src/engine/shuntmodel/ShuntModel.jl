"""
    LineCableModels.Engine.ShuntModel

Select the cable-local shunt geometry model independently of dielectric
constitutive laws. The equivalent annular layer is the default. The geometric boundary
approximation resolves exposed wires and tapes inside a layered circular reference.

# Dependencies

$(IMPORTS)
"""
module ShuntModel
import ...Commons: FormulationOptions, formulas
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: off
# Expanded in the module docstring, outside this module's analyzed expressions.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
import ..Engine: ShuntModelFormulation
import ...LineCableModels: FormulaDefinition, Expression, formula_id, description
import ...Commons: formulation_options, initialize_buffers
using ...Commons: vacuum_permittivity
import ...Materials: Material
import ...DataModel: CableDesign, CableBlueprint
import ...DataModel
import ...TextDisplay
import ..Engine: ShuntShield, ShuntJacket, tolerance_norm
import ..Engine: InsulationAdmittance, SemiconAdmittance
import ...LineCableModels: nominal, validate
using ...LineCableModels: constitutive
using LinearAlgebra: norm, opnorm, cond, lmul!, mul!, ldiv!, rdiv!, qr!, cholesky!,
                     Hermitian, eigen, SymTridiagonal, UpperTriangular
using QuadGK: alloc_segbuf, quadgk!
import SpecialFunctions

export Formula, formulas

include("interface.jl")
include("geometry.jl")
include("charge_collocation.jl")
include("blueprint.jl")

public BoundarySolveError, InternalShuntDiagnostic, DielectricGroup, dielectric_groups,
       internal_capacitance
end
