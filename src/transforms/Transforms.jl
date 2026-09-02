"""
    LineCableModels.Transforms

Transform fully coupled line-parameter matrices between phase and modal
coordinate domains independently of the backend that calculated them.

# Dependencies

$(IMPORTS)
"""
module Transforms

export ModalTransformationProblem, ModalTransformationFormulation
export LineCableModelsModal, ModalOperators, Formula
export operators, formula_id, assumptions, formulas, gamma, modal_quantities

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: FormulaMethod, nominal, FormulaSpec, _construction, validate
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  ComputationDetails, compute, computation_details,
                  computation_owner, details
import ..Engine: LineParameters, PhaseDomain, ModalDomain,
                 SeriesImpedance, ShuntAdmittance,
                 description, formula_id, selectdomain,
                 offdiagonal_ratio, reciprocity!
using LinearAlgebra: Diagonal, I, checksquare, cond, diag, dot, eigen,
                     issuccess, ldiv!, lu!, mul!, norm, rdiv!, svd
#! explicit-imports: on

include("interfaces.jl")
include("problems.jl")

"Registered modal-transformation formula selected by `:default`."
const DEFAULT = :Chrysochos2014

include("formulations.jl")
include("eigensystems.jl")

#! explicit-imports: off
const FORMULAS = let
    directory = joinpath(@__DIR__, "formulas")
    Base.include_dependency(directory)
    files = sort!(filter(
        path -> endswith(path, ".jl"),
        readdir(directory; join = true)
    ))
    identifiers = Symbol[]
    for file in files
        identifier = include(file)
        identifier isa Symbol || error(
            "modal-transformation formula file $(basename(file)) must return its Symbol identifier"
        )
        identifier in identifiers && error(
            "duplicate modal-transformation formula identifier :$identifier"
        )
        push!(identifiers, identifier)
    end
    Tuple(identifiers)
end
#! explicit-imports: on

"Return the built-in modal-transformation formula identifiers."
formulas() = FORMULAS

include("compute.jl")
include("quantities.jl")

end # module Transforms
