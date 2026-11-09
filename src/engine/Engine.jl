"""
    LineCableModels.Engine

Compute the series impedances and shunt admittances of cables per unit length. The problem
describes the cables, the earth, the operating temperature, and the frequencies. A
formulation selects the formulas, and a backend executes the computation. The trait
[`backend`](@ref) gives a formulation type's default backend. A backend passed first to
[`compute`](@ref) chooses another. Each formula family is a child module with one `Formula`
type whose identifier names the formula.
"""
module Engine

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: compute, defaults, select, validate
using ..Commons: nominal
using ..Commons: AbstractFormulation, AbstractBackend, FormulationOptions,
                 ComputationOptions
using ..Materials: Material
using ..Earth: EarthLayer, EarthModel
using ..BaseGeometry: Placement, support
using ..DataModel: CableSystem, flatten
using IntervalSets: ClosedInterval
using LinearAlgebra: BlasFloat, BlasInt, I, LAPACK, LU, SingularException, checksquare,
                     ldiv!, lu!, mul!
using QuadGK: alloc_segbuf, quadgk

export ε₀, μ₀
export AbstractInternalImpedance, AbstractInsulationImpedance, AbstractInsulationAdmittance,
       AbstractSemiconAdmittance, AbstractEarthReturn, AbstractShunt, AbstractResistivity,
       AbstractEarthProperties
export admittivity, impedance!, potential!, to_buffer!, from_buffer,
       initialize_buffers
export LineParametersProblem, CoaxialFormulation, FEM, LCMBackend, GmshBackend,
       backend, LineParameters
export @formulation, @air, @earth, @mixed
export bundle!, eliminate!, factorize!, transposition!, reductions
export SpectralIntegral, integrate
public CoaxialInput, domain, solve!, Λ
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

Represent the common supertype of the laws of earth properties and of the models of the
equivalent earth.
"""
abstract type AbstractEarthProperties <: AbstractFormulation end

"""
    admittivity(formulation, material, f, θ)

Return the admittivity of a material at a frequency and a temperature, by the law of its class
in the formulation. A conductor takes its resistivity at the temperature by the resistivity
law, with its relative permittivity. An insulator or a semicon takes the admittivity law of its
class, which does not depend on the temperature. The medium's forms, `admittivity(layer, s)`
for an earth layer and `admittivity(rho, eps_r, s)` for a resistivity and a relative
permittivity, return the conductivity plus the complex frequency times the permittivity.

# Arguments

- `formulation`: formulation that holds the material laws.
- `material`: material of a conductor or a dielectric layer.
- `f`: frequency \\[Hz\\].
- `θ`: temperature \\[°C\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity end

"""
    domain(formula)
    domain(problem, options)
    domain(input)

Return the closed interval of each quantity by its key: where a formula holds, or what a
computation presents to the formulas of its formulation.

A formula that holds within bounds declares them in its own file. Any other formula holds for
every value of every quantity, `(;)`. For each quantity, a computation presents the interval
from its least to its greatest nominal value, empty when it has none.
`validate(presented, formulation)` warns once for each quantity outside the domain of a
formula of the formulation.

| Key | Quantity | Unit | Read from |
| --- | --- | --- | --- |
| `frequencies` | the frequencies of the scan | \\[Hz\\] | `problem` |
| `mu_r_earth` | the relative permeability of each earth layer below the air | \\[dimensionless\\] | `problem` |
| `Γ` | the modulus of the longitudinal wavenumber | \\[1/m\\] | `options` |
| `mD` | the product ``m D`` of the deepest layer | \\[dimensionless\\] | `problem` |
| `Δθ` | the distance of the temperature from `T0` of each conductor material | \\[K\\] | `problem` |
| `wall` | ``(r_\\mathrm{ex} - r_\\mathrm{in}) / (r_\\mathrm{ex} + r_\\mathrm{in})`` of each tube with a hole | \\[dimensionless\\] | `input` |

The key `mD` bounds the product ``m D``, where

```math
m = \\sqrt{2 \\pi f \\mu_0 \\mu_r / \\rho}
```

takes the static resistivity and permeability of the deepest layer at the lowest and the
highest frequency. The length ``D`` is the distance between two buried cables or the outer
radius of one \\[m\\].

# Arguments

- `formula`: formula of a family.
- `problem`: line-parameters problem.
- `options`: options of a coaxial or finite-element formulation, with `Γ` \\[1/m\\].
- `input`: numerical input of a coaxial computation.

# Returns

- A named tuple of closed intervals by quantity key.
"""
function domain end

"""
    impedance!(formula, Z, …)

Write the series-impedance matrix of a formula. A method that combines formulas takes
them all before `Z`. The method on a formulation is the
recipe of the earth step, and it takes `buffers, input, media, s, k` after `Z`.

# Arguments

- `formula`: series-impedance formula.
- `Z`: destination for the matrix in \\[Ω/m\\].
- `buffers`: preallocated computation arrays for the earth step.
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

Write the potential-coefficient matrix of a formula. The method on a formulation is the
recipe of the earth step, and it takes `buffers, input, media, s, k` after `P`.

# Arguments

- `formula`: potential-coefficient formula.
- `P`: destination for the matrix in \\[m/F\\].
- `buffers`: preallocated computation arrays for the earth step.
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

Move the per-frequency terms of a formula into its block of `buffers`, at the frequency of index
`k`. A formula without such terms returns `buffers` unchanged, as this method does.

# Arguments

- `formula`: formula whose terms the recipes read.
- `buffers`: preallocated computation arrays.
- `input`: numerical cable geometry and connections.
- `media`: the air and the equivalent earth at the frequency.
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.

# Returns

- `buffers`.
"""
to_buffer!(::AbstractFormulation, buffers, input, media, s, k) = buffers

"""
    from_buffer(formula, buffers)

Return the block of `buffers` that holds the terms of `formula`, keyed by its identifier.

# Arguments

- `formula`: formula whose block the method reads.
- `buffers`: preallocated computation arrays.

# Returns

- The block of the formula, a named tuple of arrays.
"""
function from_buffer end

"""
$(TYPEDSIGNATURES)

Add the preallocated arrays that a formula needs to `buffers` once, before the loop over
frequencies. A formula without arrays of its own returns `buffers` unchanged, as this method
does.

# Arguments

- `formula`: formula whose arrays the computation needs.
- `input`: numerical cable geometry and connections.
- `buffers`: arrays already allocated.

# Returns

- `buffers`.
"""
initialize_buffers(::AbstractFormulation, input, buffers) = buffers

"""
    Λ(formula, Val(:self), κ, r, h)
    Λ(formula, Val(:mutual), κ, r, D, D′)

Return the circumferential mean of the direct and image fields of a source on a receiver in
one medium:

```math
\\Lambda_{ii} = K_0(\\kappa r) - I_0(\\kappa r) K_0(2 \\kappa h), \\qquad
\\Lambda_{ij} = I_0(\\kappa r) \\left[K_0(\\kappa D) - K_0(\\kappa D')\\right].
```

# Arguments

- `formula`: earth-return formula.
- `κ`: transverse propagation parameter of the medium \\[1/m\\].
- `r`: exterior radius of the receiver \\[m\\].
- `h`: distance from the receiver to the interface \\[m\\].
- `D`: distance from the center of the receiver to the source \\[m\\].
- `D′`: distance from the center of the receiver to the image of the source \\[m\\].

# Returns

- The mean \\[dimensionless\\].

# Notes

The earth-return family defines the method that each of its formulas gets. A formula with
another direct and image term defines its own method.
"""
function Λ end

include("constants/constants.jl")
include("problem.jl")
include("input.jl")
include("integral.jl")
include("reductions/reductions.jl")
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
include("compute.jl")

end
