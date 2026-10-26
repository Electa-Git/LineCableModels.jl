"""
    LineCableModels.Engine

Compute the series impedance and the shunt admittance per unit length of the cables of a
line. A problem describes the cables, the earth, the operating temperature and the
frequencies.
A formulation selects the formula of each part of the computation, and a backend runs it.
Each formula family is a child module with one `Formula` type, whose identifier names the
formula.
"""
module Engine

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: compute, defaults, select, validate
using ..Commons: AbstractFormulation, AbstractBackend, FormulationOptions,
                 ComputationOptions
using ..Materials: Material
using ..Earth: EarthModel
using ..DataModel: Placement, CableSystem, support, flatten
using LinearAlgebra: I, checksquare, ldiv!, lu!, mul!
using QuadGK: alloc_segbuf, quadgk

export vacuum_permittivity, vacuum_permeability
export AbstractInternalImpedance, AbstractInsulationImpedance, AbstractInsulationAdmittance,
       AbstractSemiconAdmittance, AbstractEarthReturn, AbstractShunt, AbstractResistivity,
       AbstractEarthProperties
export admittivity, resistivity, properties, impedance, impedance!, potential!, prepare!,
       initialize_buffers
export LineParametersProblem, CoaxialFormulation, FEM, LCMBackend, GmshBackend,
       backend, LineParameters
export bundle!, eliminate!, transposition!, formulas
export SpectralIntegral, integrate, nominal, ray!
public CoaxialInput
public InternalImpedance, InsulationImpedance, InsulationAdmittance, SemiconAdmittance,
       EarthReturn, Shunt, Resistivity, EarthProperties

"Supertype of the formulas of the surface impedances of a tubular conductor."
abstract type AbstractInternalImpedance <: AbstractFormulation end
"Supertype of the formulas of the magnetic term across an annular dielectric."
abstract type AbstractInsulationImpedance <: AbstractFormulation end
"Supertype of the admittivity laws of insulating layers."
abstract type AbstractInsulationAdmittance <: AbstractFormulation end
"Supertype of the admittivity laws of semiconducting layers."
abstract type AbstractSemiconAdmittance <: AbstractFormulation end
"Supertype of the formulas of the earth-return impedance and potential coefficients."
abstract type AbstractEarthReturn <: AbstractFormulation end
"Supertype of the formulas that form the potential coefficients inside a cable."
abstract type AbstractShunt <: AbstractFormulation end
"Supertype of the temperature laws of the resistivity of a conductor."
abstract type AbstractResistivity <: AbstractFormulation end
"Supertype of the frequency laws of the properties of the earth."
abstract type AbstractEarthProperties <: AbstractFormulation end

"""
    admittivity(formula, material, f)

Return the admittivity of a dielectric layer of `material` at the frequency `f` \\[Hz\\] by
`formula`, in \\[S/m\\].
"""
function admittivity end

"""
    resistivity(formula, material, temperature)

Return the resistivity of a conductor of `material` at `temperature` \\[°C\\] by
`formula`, in \\[Ω·m\\].
"""
function resistivity end

"""
    properties(formula, layer, f)

Return the resistivity \\[Ω·m\\], the relative permittivity and the relative permeability of
the earth `layer` at the frequency `f` \\[Hz\\] by `formula`, as
`(rho = …, eps_r = …, mu_r = …)`.
"""
function properties end

"""
    impedance(formula, …)

Return the series impedance per unit length that `formula` gives \\[Ω/m\\], or a named tuple
of the impedances that it gives together.
"""
function impedance end

"""
    impedance!(formula, Z, …)

Write into the matrix `Z` the series impedances per unit length that `formula` gives
\\[Ω/m\\]. A method of several formulas takes them all before `Z`.
"""
function impedance! end

"""
    potential!(formula, P, …)

Write into the matrix `P` the potential coefficients per unit length that `formula` gives
\\[m/F\\].
"""
function potential! end

"""
    prepare!(formula, buffers, input, media, s, k)

Fill in `buffers` the state of `formula` at the frequency of index `k`, which its quantity
methods read. `s = j 2π f` \\[1/s\\] is the complex frequency, and `media` gives the
properties of air and of the earth at that frequency. A formula without such state has
nothing to prepare.
"""
prepare!(::AbstractFormulation, buffers, input, media, s, k) = buffers

"""
    initialize_buffers(formula, problem, input, buffers)

Return `buffers` with the preallocated arrays that `formula` needs to compute `problem`,
built once, before the loop over the frequencies. A formula without arrays of its own adds
none.
"""
initialize_buffers(::AbstractFormulation, problem, input, buffers) = buffers

"""
    nominal(value)

Return the nominal value of a quantity: a real number itself, and a complex number with the
nominal values of its real and imaginary parts. A decision that must not depend on the
uncertainty of a value reads its nominal value. Such decisions are the branch of a kernel and
the breakpoints of an integral.
"""
nominal(value::Real) = value
nominal(value::Complex) = complex(nominal(real(value)), nominal(imag(value)))

include("constants/constants.jl")
include("problem.jl")
include("input.jl")
include("integral.jl")
include("internalimpedance/InternalImpedance.jl")
include("insulationimpedance/InsulationImpedance.jl")
include("insulationadmittance/InsulationAdmittance.jl")
include("semiconadmittance/SemiconAdmittance.jl")
include("earthreturn/EarthReturn.jl")
include("shunt/Shunt.jl")
include("resistivity/Resistivity.jl")
include("earthproperties/EarthProperties.jl")
include("formulation.jl")
include("backends.jl")
include("result.jl")
include("assembly.jl")
include("reductions/reductions.jl")
include("compute.jl")

end
