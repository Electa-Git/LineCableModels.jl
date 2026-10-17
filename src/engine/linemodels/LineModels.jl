"""
    LineCableModels.Engine.LineModels

Define frequency-dependent line-model fitting problems, formulations and results.

# Dependencies

$(IMPORTS)
"""
module LineModels

export LineModelProblem, LineModelFormulation, LineModel, WidebandModel, ModalModel
export poles, residues, delays, fitting_error, passivity, Formula, formula_id, formulas

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...LineCableModels: Expression, FormulaDefinition, formula, parameterize, validate,
                          line_length, Grid, Gridspace, verbosity, VerbosityLogger
import ...Commons: AbstractProblemDefinition, AbstractFormulation, AbstractCoreResult,
                  FormulationOptions, ComputationOptions, ComputationDetails,
                  compute, computation_options, computation_details, formulation_options,
                  details, initialize_buffers, formulas, observe, observables, basis
import ...Commons
import ...Units
import ...TextDisplay
import ..Engine: description, formula_id, frequencies, sigma, LineParametersFormulation, PhaseDomain
import ..ModalAnalysis: PropagationParameters, Ti, transform
using LinearAlgebra: norm, lu, mul!, ldiv!
using Logging: with_logger, current_logger
#! explicit-imports: on

include("problems.jl")
include("formulations.jl")
include("interfaces.jl")
public fit
include("results.jl")
include("compute.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/ulm1999.jl"),
    include("formulas/fdcm2016.jl"),
    include("formulas/marti1982.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return the built-in line-model formula identifiers.
"""
formulas(::Type{<:Formula}) = FORMULAS

include("observations.jl")
include("textdisplay.jl")

end # module LineModels
