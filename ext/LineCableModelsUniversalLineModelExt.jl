"""
    LineCableModelsUniversalLineModelExt

Fit wideband line models with UniversalLineModel.jl.
"""
module LineCableModelsUniversalLineModelExt

import LineCableModels
import UniversalLineModel
using LineCableModels: PropagationParameters, PhaseDomain, frequencies, line_length,
                       details, Yc, H, gamma, Ti, LineModel
import LineCableModels.Engine.LineModels: Formula, fit, WidebandModel
import LineCableModels: compute
import LineCableModels.Engine.TransientAnalysis
import LineCableModels.Engine.TransientAnalysis: TransientAnalysisProblem,
    TerminalImpedance, Waveform, Step, Cosine, DoubleExponential, Matched
using UniversalLineModel: LineFrequencyData, ULMModel, fit_ulm, fit_fdcm,
                          eval_yc, eval_h, passivity_check, enforce_passivity,
                          Network, add_line!, add_resistor!, add_inductor!, add_vsource!,
                          add_rational_admittance!,
                          StepSource, CosineSource, DoubleExpSource, simulate
using Logging: @warn, @info

include("LineCableModelsUniversalLineModelExt/adapter.jl")
include("LineCableModelsUniversalLineModelExt/result.jl")
include("LineCableModelsUniversalLineModelExt/passivity.jl")
include("LineCableModelsUniversalLineModelExt/fit.jl")
include("LineCableModelsUniversalLineModelExt/transients.jl")

end
