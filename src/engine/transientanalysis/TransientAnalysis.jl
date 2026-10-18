"""
    LineCableModels.Engine.TransientAnalysis

Define transient line problems, terminal assignments and voltage and current responses.

# Dependencies

$(IMPORTS)
"""
module TransientAnalysis

export TerminalImpedance, Waveform, Step, Cosine, DoubleExponential, Matched
export TransientAnalysisProblem, TransientFormulation, TransientResponse
export voltage, current, Formula, formula_id, formulas
export peak, arrival_time, damping, dc_offset

#! explicit-imports: off
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ...LineCableModels: Expression, FormulaDefinition, formula, parameterize, validate
import ...LineCableModels: verbosity, VerbosityLogger, line_length
import ...Commons: AbstractProblemDefinition, AbstractFormulation, AbstractCoreResult,
                  FormulationOptions, ComputationOptions, ComputationDetails,
                  compute, computation_options, computation_details, formulation_options,
                  details, initialize_buffers, formulas, observe, basis
import ...Commons
import ...Units
import ..Engine: LineParameters, LineParametersProblem, LineParametersFormulation, Formulation,
                 description, formula_id, compare, RMSError, completed_inputs, retain_gridpoint
import ..LineModels: LineModel, ModalModel
import ...TextDisplay
using Logging: with_logger, current_logger
using LinearAlgebra: Diagonal, eigen!, lu, lu!, ldiv!, mul!
using Statistics: mean
#! explicit-imports: on

include("terminals.jl")
include("problems.jl")
include("formulations.jl")
include("interfaces.jl")
public respond!
include("results.jl")
include("kpis.jl")
include("compute.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/nlt.jl"),
    include("formulas/td.jl"),
    include("formulas/default.jl"),
)
#! explicit-imports: on

"""
Return the built-in transient formula identifiers.
"""
formulas(::Type{<:Formula}) = FORMULAS

include("observations.jl")
include("textdisplay.jl")

end # module TransientAnalysis
