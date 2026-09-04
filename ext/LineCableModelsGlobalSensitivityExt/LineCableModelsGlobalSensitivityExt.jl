"""
    LineCableModelsGlobalSensitivityExt

Evaluate variance-based Sobol sensitivity products through GlobalSensitivity.jl.
"""
module LineCableModelsGlobalSensitivityExt

import Distributions
import GlobalSensitivity
import QuasiMonteCarlo
import Statistics

import LineCableModels
const Grammar = LineCableModels.Grammar
const ParametricBuilder = LineCableModels.ParametricBuilder
const UQ = LineCableModels.UQ
import LineCableModels.Grammar: compute

include("sobol.jl")

end
