module GetDPFrontend

const GetDP = @__MODULE__

# Include helper functions
include("helpers.jl")

# Export helper functions
export add_raw_code, comment, make_args, get_getdp_executable

# Include base classes
include("getdp_object.jl")

# Export base classes
export AbstractGetDPObject, GetDPObject, SimpleItem, Base_, ObjectItem, CaseItem_, Case_
export code, add_raw_code!, add_comment!, add!

# Include main components
include("group.jl")
include("function.jl")
include("constraint.jl")
include("functionspace.jl")
include("jacobian.jl")
include("integration.jl")
include("formulation.jl")
include("resolution.jl")
include("postprocessing.jl")
include("postoperation.jl")
include("problem_definition.jl")

# Export public API
export Constraint, Formulation, Function, FunctionSpace, Group, Integration
export Jacobian, PostOperation, PostProcessing, Problem, Resolution
export define!, add!, add_list!, add_file!, add_analytic!, add_piecewise!, add_akima!
export add_quantity!, add_equation!, content, add_subspace!
export Region, Global, NodesOf, EdgesOf, FacetsOf, VolumesOf, ElementsOf
export GroupsOfNodesOf, GroupsOfEdgesOf, GroupsOfEdgesOnNodesOf, GroupOfRegionsOf
export EdgesOfTreeIn, FacetsOfTreeIn, DualNodesOf, DualEdgesOf, DualFacetsOf, DualVolumesOf
export get_code, make_problem!, write_file, include!, SystemItem#, write_multiple_problems

export add_constraint!, add_global_quantity!, error, add_operation!, add_basis_function!
export add_case!, VolSphShell, add_nested_case!, add_space!, add_constant!, assign!, case!
export for_loop!

# Math Functions
export Exp, Log, Log10, Sqrt
export Sin, Asin, Cos, Acos, Tan, Atan, Atan2
export Sinh, Cosh, Tanh, TanhC2
export Fabs, Abs, Floor, Ceil, Fmod, Min, Max, Sign
export Jn, dJn, Yn, dYn
export Cross, Hypot, Norm, SquNorm, Unit, Transpose, Inv, Det, Rotate, TTrace
export Cos_wt_p, Sin_wt_p
export Printf, Rand
export Normal, NormalSource, Tangent, TangentSource
export ElementVol, SurfaceArea, GetVolume, CompElementNum, GetNumElements, ElementNum
export QuadraturePointIndex, AtIndex
export InterpolationLinear, dInterpolationLinear
export InterpolationBilinear, dInterpolationBilinear
export InterpolationAkima, dInterpolationAkima
export Order, Field
export ScalarField, VectorField, TensorField
export ComplexScalarField, ComplexVectorField, ComplexTensorField
export GetCpuTime, GetWallClockTime, GetMemory
export SetNumberRunTime, GetNumberRunTime
export SetVariable, GetVariable
export ValueFromIndex, VectorFromIndex, ValueFromTable

end # module
