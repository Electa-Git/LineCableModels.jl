"""
    LineCableModels.Engine

Compute the series impedances and shunt admittances of cables per unit length. A problem
describes the cables, the earth, the operating temperature, and the frequencies. A
formulation selects the formulas, and a backend executes the computation. The trait
[`backend`](@ref) gives a formulation type's default backend. A backend passed first to
[`compute`](@ref) chooses another. Each formula family is a child module with one `Formula`
type whose identifier names the formula.
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
export bundle!, eliminate!, transposition!, formulas, reductions
export SpectralIntegral, integrate, nominal, ray!
public CoaxialInput, solve!
public InternalImpedance, InsulationImpedance, InsulationAdmittance, SemiconAdmittance,
       EarthReturn, Shunt, Resistivity, EarthProperties

"""
$(TYPEDEF)

Represent the common supertype of surface-impedance formulas for tubular conductors.
"""
abstract type AbstractInternalImpedance <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of magnetic series-impedance formulas across annular
dielectrics.
"""
abstract type AbstractInsulationImpedance <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of admittivity laws for insulator layers.
"""
abstract type AbstractInsulationAdmittance <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of admittivity laws for semicon layers.
"""
abstract type AbstractSemiconAdmittance <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of earth-return impedance and potential-coefficient formulas.
"""
abstract type AbstractEarthReturn <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of formulas for potential coefficients inside cables.
"""
abstract type AbstractShunt <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of temperature laws for conductor resistivity.
"""
abstract type AbstractResistivity <: AbstractFormulation end
"""
$(TYPEDEF)

Represent the common supertype of frequency laws for earth properties.
"""
abstract type AbstractEarthProperties <: AbstractFormulation end

"""
    admittivity(formula, material, f)

Return the admittivity of a dielectric layer at a frequency.

# Arguments

- `formula`: admittivity law.
- `material`: material of the dielectric layer.
- `f`: frequency \\[Hz\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity end

"""
    resistivity(formula, material, temperature)

Return the conductor resistivity at a temperature.

# Arguments

- `formula`: resistivity law.
- `material`: material of the conductor.
- `temperature`: temperature of the conductor \\[°C\\].

# Returns

- Resistivity \\[Ω·m\\].
"""
function resistivity end

"""
    properties(formula, layer, f)

Return the material properties of an earth or air layer at a frequency.

# Arguments

- `formula`: earth-properties law.
- `layer`: earth or air layer.
- `f`: frequency \\[Hz\\].

# Returns

- `rho`: resistivity \\[Ω·m\\].
- `eps_r`: relative permittivity \\[dimensionless\\].
- `mu_r`: relative permeability \\[dimensionless\\].
"""
function properties end

"""
    impedance(formula, …)

Return the series-impedance terms of a formula.

# Arguments

- `formula`: series-impedance formula.
- `…`: physical inputs of the selected method.

# Returns

- An impedance \\[Ω/m\\] or a named tuple of impedance terms in the same unit.
"""
function impedance end

"""
    impedance!(formula, Z, …)

Write the series-impedance matrix of a formula. A method that combines formulas takes
them all before `Z`. An earth-return method takes
`buffers, input, media, s, k` after `Z`.

# Arguments

- `formula`: series-impedance formula.
- `Z`: destination for the matrix in \\[Ω/m\\].
- `buffers`: preallocated computation arrays for an earth-return method.
- `input`: numerical cable geometry and connections.
- `media`: air and earth properties at the frequency.
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.
- `…`: inputs of another selected method.

# Returns

- `Z`, with the computed entries \\[Ω/m\\].
"""
function impedance! end

"""
    potential!(formula, P, …)

Write the potential-coefficient matrix of a formula. An earth-return method takes
`buffers, input, media, s, k` after `P`.

# Arguments

- `formula`: potential-coefficient formula.
- `P`: destination for the matrix in \\[m/F\\].
- `buffers`: preallocated computation arrays for an earth-return method.
- `input`: numerical cable geometry and connections.
- `media`: air and earth properties at the frequency.
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.
- `…`: inputs of another selected method.

# Returns

- `P`, with the computed entries \\[m/F\\].
"""
function potential! end

"""
$(TYPEDSIGNATURES)

Fill the shared state of a formula in `buffers` at the frequency of index `k`. Its quantity
methods read this state. A formula without such state returns `buffers` unchanged, as this
method does.

# Arguments

- `formula`: formula whose state needs preparation.
- `buffers`: preallocated computation arrays.
- `input`: numerical cable geometry and connections.
- `media`: air and earth properties at the frequency.
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.

# Returns

- `buffers`.
"""
prepare!(::AbstractFormulation, buffers, input, media, s, k) = buffers

"""
$(TYPEDSIGNATURES)

Add the preallocated arrays that a formula needs to compute `problem` to `buffers` once,
before the loop over frequencies. A formula without arrays of its own returns `buffers`
unchanged, as this method does.

# Arguments

- `formula`: formula whose arrays the computation needs.
- `problem`: line-parameters problem.
- `input`: numerical cable geometry and connections.
- `buffers`: arrays already allocated.

# Returns

- `buffers`.
"""
initialize_buffers(::AbstractFormulation, problem, input, buffers) = buffers

"""
$(TYPEDSIGNATURES)

Return a real number unchanged or a complex number with nominal real and imaginary parts.
Kernel branches and integration breakpoints use nominal values to exclude uncertainty
from their decisions.

# Arguments

- `value`: numeric quantity.

# Returns

- The nominal value, in the unit of `value`.
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
