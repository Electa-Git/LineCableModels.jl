"""
    LineCableModelsPolyChaosExt

Fit and independently validate non-intrusive polynomial-chaos expansions for
LineCableModels parametric cable calculations.
"""
module LineCableModelsPolyChaosExt

import PolyChaos
using LinearAlgebra: norm
import Random

import LineCableModels
import LineCableModels: basis, compute, nominal, observe, uncertainty
import LineCableModels: points, realize, realize_arguments, uncertainties
import LineCableModels.Grammar: computation_details
import LineCableModels.Engine
using LineCableModels: CableConstants, G, L, LineParameters, R, C
using LineCableModels: ParametricProblem, PolynomialChaos, PolynomialChaosResult

include("basis.jl")
include("projection.jl")
include("compute.jl")

end
