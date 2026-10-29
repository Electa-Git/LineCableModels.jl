# Matrix formulation of line and cable parameters

LineCableModels calculates line and cable parameters by assembling the per-unit-length series-impedance matrix ``\mathbf Z`` and shunt-admittance matrix ``\mathbf Y``. These matrices contain the conductor self and mutual coefficients in the multiconductor transmission-line equations. They determine the propagation of voltages and currents along the line.

The matrix entries combine contributions from conducting materials, insulation and the external medium. Each contribution is evaluated with a selected formulation. Any formulation expressible in this framework can be incorporated with compatible conductor definitions, reference conventions and physical assumptions. The formulation sections linked below document the available equations and their domains of validity.

## Transmission-line equations

At angular frequency ``\omega``, with time dependence ``e^{j\omega t}``, the multiconductor transmission-line equations are

```math
\frac{\mathrm d\mathbf V}{\mathrm dx}=-\mathbf Z\mathbf I,
\qquad
\frac{\mathrm d\mathbf I}{\mathrm dx}=-\mathbf Y\mathbf V.
```

The coordinate ``x`` measures distance along the line. The matrices ``\mathbf Z`` and ``\mathbf Y`` have units ``\Omega/\mathrm m`` and ``\mathrm S/\mathrm m``, respectively. Their order equals the number of retained electrical conductors [Ametani1980, Ametani2015b](@cite).

In the conductor representation, each component of ``\mathbf V`` and ``\mathbf I`` corresponds to one retained conductor. The same external reference is used for all voltages, and all currents are positive in the same longitudinal direction. Matrix indices enumerate conductors independently of their physical cross-sectional coordinates.

Sheaths and armor remain electrical conductors when retained in this representation. Bonding conditions, return-current constraints, and modal transformations are applied separately from parameter assembly.

### Potential coefficients and shunt admittance

Shunt formulations often provide Maxwell potential coefficients rather than admittance entries. Under electrostatic, lossless-dielectric assumptions, the potential-coefficient matrix ``\mathbf P`` relates conductor voltages to per-unit-length charges ``\mathbf q``:

```math
\mathbf V=\mathbf P\mathbf q,
\qquad
\mathbf C=\mathbf P^{-1},
\qquad
\mathbf Y=j\omega\mathbf P^{-1}.
```

The units of ``\mathbf P`` and ``\mathbf C`` are ``\mathrm m/\mathrm F`` and ``\mathrm F/\mathrm m``, respectively. The admittance calculation requires the inverse of the complete potential-coefficient matrix [Ametani1980, Ametani2015b](@cite). Matrix inversion couples the coefficients, so ``Y_{jk}`` cannot be calculated as ``j\omega/P_{jk}``.

The normalization of a shunt formulation determines the conversion to ``\mathbf Y=\mathbf G+j\omega\mathbf C`` [Gustavsen2005](@cite). Under the charge-based, lossless definition above, ``\mathbf P=j\omega\mathbf Y^{-1}`` for invertible matrices and ``\omega\ne0``. The potential-coefficient matrix is distinct from an inverse-admittance matrix.

## Physical contributions and formulation families

The formulation families correspond to the physical contributions required by the impedance and admittance calculations.

| Physical contribution | Matrix contribution | Formulation section |
| --- | --- | --- |
| Electromagnetic response inside conducting materials | Local series response, including metallic surface and through-wall transfer impedances | [Conductors](conductors.md) |
| Fields in insulation, jackets, and semiconducting regions | Series-field contributions and local potential-coefficient or shunt contributions | [Dielectrics](dielectrics.md) |
| Series fields outside the local cable assemblies | External self and mutual impedances, including earth return | [Earth return](earth.md) |
| Electric fields in the external medium | External potential coefficients used in the shunt-admittance calculation | [Earth return](earth.md) |
| Earth conductivity, permittivity, and their frequency dependence | Material inputs to the external-field formulations | [Earth return](earth.md) |

Each contribution has a selectable formulation. Formulations used together require compatible geometry and material properties. Their voltage references and current conventions must agree, as must their propagation assumptions.

A formulation may provide a scalar coefficient, a coupled block or a complete matrix. Its output must correspond to the physical contribution and conductor representation required by the assembly. An earth-return correction requires the associated geometric term to form a complete external impedance. Radial branch admittances and loop impedances require conversion before use in a conductor matrix. Contributions already included in one term must not be repeated in another.

## Block matrix assembly

In the cylindrical construction considered here, each local unit consists of one round conductor or a concentric assembly of metallic conductors and dielectric regions. The shared regions comprise the air and earth outside the individual units. The following arbitrary-dimension notation reformulates the block patterns in [Ametani1980, Ametani2015b](@cite). The cylindrical assumptions apply to this decomposition, rather than to the transmission-line matrix formulation itself.

Consider ``m`` local units. Unit ``j`` contains ``N_j`` longitudinal metallic conductors, ordered from the center outward when the unit is coaxial. The conductor count is

```math
n=\sum_{j=1}^{m}N_j.
```

A coaxial unit may contain a solid central conductor followed by any number of concentric annular conductors, separated by dielectric regions. The values of ``N_j`` need not be identical across the installation.

Let ``\mathbf B_j`` and ``\mathbf D_j`` denote the local series-impedance and potential-coefficient blocks of unit ``j``, respectively. Both have order ``N_j``. Define

```math
\mathbf B=\operatorname{blockdiag}(\mathbf B_1,\ldots,\mathbf B_m),
\qquad
\mathbf D=\operatorname{blockdiag}(\mathbf D_1,\ldots,\mathbf D_m).
```

Each local block contains the self and mutual coefficients of the conductors within its unit. The block-diagonal structure excludes interactions between different units from ``\mathbf B`` and ``\mathbf D``.

### Current and charge aggregation

Define the conductor-membership matrix

```math
\mathbf S=\operatorname{blockdiag}
\left(\mathbf 1_{N_1},\ldots,\mathbf 1_{N_m}\right)
\in\mathbb R^{n\times m},
```

where ``\mathbf 1_{N_j}`` is a column of ones. Column ``j`` identifies the conductors belonging to unit ``j``. The total longitudinal current of each unit is

```math
\mathbf J=\mathbf S^{\mathsf T}\mathbf I,
\qquad
J_j=\sum_{k=1}^{N_j}I_{j,k}.
```

The corresponding total per-unit-length charge is

```math
\mathbf Q=\mathbf S^{\mathsf T}\mathbf q,
\qquad
Q_j=\sum_{k=1}^{N_j}q_{j,k}.
```

Under the cylindrical approximation, the shared-region formulations use ``J_j`` and ``Q_j`` as source amplitudes. Multiplication by ``\mathbf S`` assigns the resulting voltage gradient or potential contribution to each conductor in the corresponding unit.

For a shared-region impedance matrix ``\mathbf H_Z`` and potential-coefficient matrix ``\mathbf H_P``, the conductor-matrix contributions are

```math
\Delta\mathbf Z=\mathbf S\mathbf H_Z\mathbf S^{\mathsf T},
\qquad
\Delta\mathbf P=\mathbf S\mathbf H_P\mathbf S^{\mathsf T}.
```

The same membership matrix applies to both quantities. For impedance, the block between units ``j`` and ``k`` is

```math
\left(\mathbf S\mathbf H_Z\mathbf S^{\mathsf T}\right)_{jk}
=(H_Z)_{jk}\mathbf 1_{N_j}\mathbf 1_{N_k}^{\mathsf T}.
```

The repeated coefficients result from representing each unit by its aggregate source amplitude [Ametani1980, Ametani2015b](@cite).

Angular current or charge distributions omitted by the cylindrical approximation require a more detailed coupled formulation. When that formulation couples several units, the coupled response must remain in the combined matrix rather than in independent local blocks.

### Assembly

Let ``\mathbf Z_{\mathrm{env}}`` and ``\mathbf P_{\mathrm{env}}`` describe the shared medium outside the individual cable boundaries. The complete matrices are

```math
\boxed{
\begin{aligned}
\mathbf Z&=\mathbf B+\mathbf S\mathbf Z_{\mathrm{env}}\mathbf S^{\mathsf T},\\
\mathbf P&=\mathbf D+\mathbf S\mathbf P_{\mathrm{env}}\mathbf S^{\mathsf T}.
\end{aligned}
}
```

The local blocks contain the contributions within each unit. The external matrices contain the self and mutual contributions of the surrounding medium.

Under the lossless potential-coefficient convention,

```math
\boxed{
\mathbf Y=j\omega
\left(\mathbf D+\mathbf S\mathbf P_{\mathrm{env}}\mathbf S^{\mathsf T}\right)^{-1}.
}
```

Inversion applies to the full sum. In general, independently inverting the local and external potential blocks and adding their admittances does not recover ``\mathbf Y``.

## The assembly in this release

The coaxial formulation fills the local blocks ``\mathbf B_j`` and ``\mathbf D_j`` of each
cable with two methods. The loop recursion assembles the surface impedances of the conductors
and the magnetic terms of the dielectric intervals into ``\mathbf B_j``. The formula
`Shunt.Formula{:series}` fills ``\mathbf D_j`` by the series chain of the potential
coefficients of the layers. The earth-return formula gives ``\mathbf Z_{\mathrm{env}}`` and
``\mathbf P_{\mathrm{env}}``, the matrices `Zg` and `Pg` by cable, as
[Earth return](earth.md) describes.

```@meta
CurrentModule = LineCableModels.Engine
```

```@docs
impedance!(internal::InternalImpedance.Formula, insulation::InsulationImpedance.Formula, Zin::AbstractMatrix{Complex{T}}, input::CoaxialInput{T}, rho::AbstractVector{T}, s::Complex{T}) where {T <: Real}
```

```@autodocs
Modules = [Shunt]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

```@docs
potential!(::Shunt.Formula{:series}, Pin::AbstractMatrix{Complex{T}}, input::CoaxialInput{T}, kappa::AbstractVector{Complex{T}}, s::Complex{T}) where {T <: Real}
```

## Reductions

The reductions act on the assembled matrices in order. The bundle first merges the conductors
that share a connection index, and the Kron elimination then removes the conductors at zero
voltage. When selected, the ideal transposition acts on `Z` and `Y` after the admittance is
computed.

```@docs
reductions
bundle!
eliminate!
transposition!
```

## The scan

The computation of [`LCMBackend`](@ref) evaluates the formulas at each frequency and writes
`Z` and `Y`:

```@docs
solve!(Z::AbstractArray{Complex{T}, 3}, Y::AbstractArray{Complex{T}, 3}, trace, problem::LineParametersProblem{T}, input::CoaxialInput{T}, formulation::CoaxialFormulation, buffers, options::ComputationOptions) where {T <: Real}
```

The formula families extend these functions of the engine. A formula method takes the formula
first.

```@docs
impedance
impedance!
potential!
prepare!(::AbstractFormulation, buffers, input, media, s, k)
initialize_buffers(::AbstractFormulation, problem, input, buffers)
initialize_buffers(formulation::CoaxialFormulation, problem::LineParametersProblem, input::CoaxialInput{T}) where {T <: Real}
admittivity
admittivity(formulation::Union{CoaxialFormulation, FEM}, material::Material{T, :insulator}, f) where {T <: Real}
resistivity
properties
```

[Back to Contents](contents.md)
