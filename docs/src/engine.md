# Computational engine

LineCableModels separates the physical problem, selected equations, and numerical
execution. Source equations and bibliography belong in the implementing formula
file. Backend/formula methods select an implementation through Julia dispatch.

Every family owns a concrete `:default` identifier. There is no author alias or
forwarding map. Author registrations preserve comparison formulas and the former earth defaults:

| Family | Registered choices |
|---|---|
| Internal impedance, insulation impedance, pipe impedance | `:default` |
| Insulation admittance, semicon admittance | `:default`, `:Ametani2004` |
| Earth impedance | `:default`, `:Carson1926`, `:Pollaczek1926`, `:Gary1976`, `:WedepohlWilcox1973`, `:Saad1996`, `:Ametani2009`, `:Lucca1994`, `:Wise1934`, `:Xue2018` |
| Earth admittance | `:default`, `:Pollaczek1926`, `:Wise1948`, `:Xue2018` |
| Frequency-dependent soil properties, equivalent earth, modal transformation | `:default` |

The internal default retains Schelkunoff's tubular conductor expressions;
insulation impedance retains the annular magnetic term documented by Ametani.
Both earth defaults implement the supplied circumferentially averaged framework
with complete enclosed-current normalization. The former air defaults remain
available as `:Wise1934` for impedance and `:Wise1948` for potential coefficients;
both former buried defaults are `:Xue2018`. Dielectric defaults are lossless; explicit `:Ametani2004`
retains conductivity and the material's supplied polarization losses. The FrequencyDependent
default preserves static properties, EquivalentHomogeneous selects the basement when explicitly
requested, and the modal default performs Levenberg–Marquardt tracking.
Bibliography remains attached to the equations despite the package-owned names.

```julia
selected = Formulation(
    earth_impedance=formula(:default;
        options=(integration=(method=:quad, options=(rtol=1e-8,)),)),
    earth_admittance=formula(:default),
    earth_properties=formula(:default),
)
result = compute(problem, selected; options=(trace=true,))
```

`FormulaDefinition` carries the identifier, explicit physical `parameters`,
`hooks`, numerical `options`, and optional formula-local `equivalent_earth`.
The receiving family resolves the selection. A bound `Functor` separates physical
`state`, concrete callables, interaction `binding`, and normalized numerical
`options`. Mutable integration storage belongs to the calculation workspace.

Numerical requirements are declared by `computation_options(binding::FormulaMethod)`.
The binding includes the owning equation function, identifier and semantic selectors.
An empty declaration admits no numerical controls. A missing declaration is an error.
`Engine.hooks(binding)` declares the external case's physical hook defaults and
admitted overrides; it does not declare equation availability.
Special functions, algebraic approximations, spectral integrals and iterative solvers
therefore share a grammar without a binary “closed/integral” classification.

An override follows one route from `formula(...; hooks=(... ,))` through
`Formulation`, indexed validation and functor construction to execution. For example:

```julia
my_Γ(jω, materials, layers) = zero(jω)
selected = Formulation(earth_impedance=formula(:default; hooks=(Γ=my_Γ,)))
result = compute(problem, selected)
```

`Γ(jω, materials, (s,t))` returns one finite scalar [1/m]. Its square is derived
from that scalar. An explicit problem `Γ` and explicit Γ hook conflict, including
when both return zero. The new earth defaults accept a common prescribed Γ; all retained author earth
equations require Γ=0. No propagation-constant root solver is implied. Medium propagation laws have signature
`air/earth(jω, μ, σ, ε)`, the permeability hook has signature `permeability(μ)`, and
a complete earth contribution has signature `contribution(functor, pair, workspace)`.
A hook unused by the selected indexed equation is rejected. Hook arities are not
guessed or retried. Modified selections are recorded in result details.

Scalar families expose `hooks=(contribution=my_law,)` with the signature documented
by their `Formula` constructor. Internal impedance exposes its inner, outer and
mutual surface callables. Unknown parameters or hooks fail at construction or indexed preflight. External
hook names and numerical sections are admitted together by the required cases.

`EarthPair` carries conductor row/column indices, integer source/target layer indices,
heights, horizontal separation and an explicit self radius. Self means the same
conductor; distinct conductors in one layer remain mutuals. A self pair has zero
horizontal separation, with its radius supplied separately. The geometry substitution
needed by a published self expression occurs at equation evaluation.

The external equation signatures are:

```julia
earth_impedance(::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, functor, pair, workspace)
earth_potential_coefficient(::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, functor, pair, workspace)
```

`Kind` is `:self` or `:mutual`; source `S` is the matrix column, and target `T`
is the row. Layer 1 is air; soils occupy layers 2 through N. The shared `validate`
uses native method selection on this canonical signature and excludes the throwing
fallback. Domain-defining methods accept the three runtime payloads without extra
subtype constraints. Numerical specializations can optimize an admitted case.
Adding a canonical equation changes admission without a second capability table.

Each earth slot also accepts a NamedTuple for a physical air/soil two-half-space model:

```julia
selected = Formulation(earth_impedance = (
    air = formula(:Carson1926),
    earth = formula(:Pollaczek1926),
    mixed = formula(:Lucca1994),
))
```

The same syntax applies independently to `earth_admittance`. The labels resolve
`(1,1)`, `(2,2)` and the two cross-layer directions before the existing indexed
method dispatch. Each case retains its own formula, hooks, numerical options and
material preparation. There is no combined formula identity. An unused leaf does
not supply missing cases or execute a kernel. True layered inputs require a scalar
selection; scalar multilayer and explicit EHEM behavior are unchanged.

The default potential formulation supplies both mixed directions. Selecting a
default leaf for only part of a matrix still assembles its complete auxiliary
system and then selects the requested final entries. Retained author entries
never become inputs to its K/H kernels.

The workspace binds every required ordered pair before frequency evaluation.
Geometry and layer indices follow the same order. Both directions are evaluated;
assembly, reduction and modal transformation preserve the returned ordered entries.
No implicit reciprocity operation supplies a missing equation or averages its result.

The medium inventory is a separate physical restriction. A homogeneous formula
consumes exactly air and one soil half-space. A finite-layer model consumes its
whole declared inventory and interfaces. Explicit `Val(S), Val(T)` methods describe
its cases; arbitrary-layer Green-function generation remains deferred. Buried
placement in a vertical multilayer earth is rejected because its physical layer
indexing has no defined origin in the present geometry contract.

Carson admits only `(1,1)`. Both Pollaczek families admit only `(2,2)` and reject
air or mixed pairs. Ametani2009 and Lucca1994 remain mixed-only equations and cannot
assemble a complete native matrix by themselves. Their missing self terms are
never filled by another source. The defaults supply air/air, earth/earth and both mixed directions with
independent medium permeabilities. Every circumference must lie wholly in its
half-space and exterior circles must not overlap. The default requires homogeneous
earth or an explicitly globally consistent equivalent earth. Pair-dependent
reductions remain available through author selections.

The complete earth source kernels satisfy

```math
K=\mathcal Z-\Gamma^2\mathcal P_\phi/s,\qquad
L=A_r^{-1}-F_rK,\qquad
P_eL=H,\quad Z_eL=K+\Gamma^2H/s,\quad Y_eH=sL.
```

Here `s=jω`, Ze is in Ω/m, Pe in m/F and Ye in S/m. Matrices are solved on the
right, with source columns and receiver rows. Pe and Ye retain the direction of
the voltage path; they are not symmetrized. The default reference is common deep
earth. Select `parameters=(reference=:interface,)`, a positive finite reference
depth below every circumference, or `reference=:scalar` for the distinct scalar
potential diagnostic. Use the same physical reference for both owners when a
consistent Ze/Pe pair is required.

Public `compute` prepares the complete system. Calling a default pair callback
without that context is an error; explicit author formulas retain their pairwise
contract. Contribution overrides that declare integration resources receive the
prepared context at the final-entry stage. Internal conductor and insulation
contributions retain the existing cable composition and terminal reductions.
With `options=(trace=true,)`, `trace.Zg` and `trace.Pg` expose the exterior
matrices; the returned total line admittance also includes insulation effects.

Select an equivalent homogeneous earth on each consuming formula:

```julia
selected = Formulation(
    earth_impedance=formula(:default;
        equivalent_earth=formula(:default; order=:before)),
    earth_admittance=formula(:default;
        equivalent_earth=formula(:default; order=:after)),
    earth_properties=formula(:default),
)
```

The reduction receives the physical `(kind,s,t)` selectors, all physical layer
properties, model, pair and frequency. Its runtime suffix is
`(rho, eps_r, mu_r, model, pair, frequency, parameters, options, workspace)` and
its result is one `EarthMaterial`. It owns its numerical sections independently
of the external equation. The consuming source explicitly admits compatible
reductions. A full multilayer consumer rejects reductions.

The current reduction default explicitly selects the bottommost soil; it is a
layer-selection settings, not a derived general recursion. `:after` applies the one
backend-selected frequency law to physical layers first. `:before` reduces static
properties and applies that same law to the resulting material. Physical and
effective pairs remain distinct in the binding, and reductions run for each
ordered interaction on which they depend. Layerwise evaluated properties are reused
between consumers when needed; the air material remains static.

A complete contribution override must also declare its numerical defaults:

```julia
using LineCableModels: FormulaMethod, computation_options
const II = LineCableModels.Engine.InternalImpedance
my_outer(functor, workspace) = zero(functor.state.jω)
computation_options(
    ::FormulaMethod{:default,typeof(II.internal_impedance),Tuple{Val{:outer}}},
    ::typeof(my_outer),
) = (;)
```

This replaces only an admitted case and receives the canonical runtime suffix.
An algebraic replacement of an integral rejects unused integration controls;
an integral replacement declares its own integration section. Small physical
hooks retain the operation they customize. No callback inherits an unrelated
provider's numerical contract or expands its physical domain.

`InternalImpedance.surface_impedances(resolved_formula, r_in, r_ex, rho, mu_r, jω)`
returns `(inner,outer,mutual)` coefficients in Ω/m, with hooks and per-kind
numerical options applied. Internal kinds have no earth-layer selectors. Assemblers
own the current-basis transformation and matrix placement. The deferred pipe
contribution concerns one contained metal and its enclosing pipe; recursive
assembly and pipe equations are outside this implementation.

`ComputationOptions` remains an alias for `NamedTuple`. The existing
`computation_options` constructor validates and normalizes execution controls once:

| `integration.method` | Numerical operation |
|---|---|
| `:quad` | Adaptive QuadGK on an admissible contour with explicit physical feature breakpoints and error estimates. |
| `:trapz` | Double-exponential transformed trapezoids from DoubleExponentialFormulas, with independently verified local subdomains. |
| `:cim` | Matrix-pencil/GPOF exponential fit along the spectral coordinate at fixed physical frequency; analytic image integration. |

`SpectralIntegral` exposes the kernel, analytic weight, spectral scale and admissible
contour. Cosine and radial Sommerfeld weights use different analytic image identities.
The overhead potential kernel additionally declares an exact simple-pole contribution;
CIM fits only its rationalized remainder. Cases without a declared integral receive
no integration section or integration scratch.

CIM fits spectral λ (or the explicitly declared radial spectral coordinate), never
physical frequency. Built-in earth kernels retain up to 64 image representations
in the computation workspace. Reuse requires an exact copied material/frequency/Γ
identity, the same kernel and transformation, and an error certificate meeting
the requested tolerance. Material-only radial responses can share fits across
depths, separations and scalar source normalizations. Geometry-dependent kernels
retain those parameters in their identity. Arbitrary callbacks are not cached.
Uniform local pencil regions generate a global exponential representation;
construction starts with a compact region set and expands it when verification
requires more work. Sample and Hankel buffers are reused.
Amplitude fitting uses the integration contour and analytic tail penalties;
an independent integral of the absolute weighted residual checks its complete
continuation. Quadrature verifies the image sum and never supplies a value
labelled `:cim`. Geometry certificates integrate an absolute residual envelope
over the complete contour, including the tail. They apply to greater weight
heights and smaller separation-plus-radius on that contour. New geometries
outside the certified range require another certificate or fit; valid cache hits
evaluate the images without quadrature or refitting. One kernel prototype
evaluation still checks the physical scalar contract. Reported certificates
remain numerical estimates.
An unresolved value-only integral raises an error. The complete
earth assembler instead consumes the numerical estimate, propagates prefactors
and right-solve sensitivity, and tightens the influential interactions until
individual Ze/Pe/Ye entries meet their budgets. A small correction can therefore
receive a larger integral tolerance without relaxing final-matrix accuracy.

CIM accepts Float32/Float64 physical inputs. Its fits and residual verification
use Float64 arithmetic, with output rounding included in the error estimate.
BigFloat and Measurements inputs use quadrature or trapz; CIM rejects them
explicitly. Trapz retains correlated physical uncertainty while its numerical
norm and sampling coordinates are nominal. A zero nominal value with nonzero
uncertainty remains visible to refinement.

For trapz, `samples` controls the initial transformed-rule spacing
(`h0=min(1,128/samples)`), `max_refinements` limits DE levels, and
`max_tail_refinements` limits independently verified subdomain refinements over
the mapped semi-infinite interval. Local phase and envelope variation guide
initial subdivision. Their finite seeding horizon does not truncate the
integral: the final interval still extends to infinity. Independent panel
quadrature references supply normalization and verification in one pass; their
absolute errors cannot cancel across panels. Successful panels are retained;
the DE package reuses nested samples within each panel. Rule tables are reused for identical
spacing, level limits and precision. Mathematical feature metadata is required to guard
against known narrow structures; finite black-box samples cannot certify the
absence of an undeclared feature. Reported errors remain estimates unless the
contributing regularity and tail information supplies bounds.

Formula discovery includes sorted `formulas/*.jl` files, each returning one unique
identifier. `FormulaMethod` binds that identifier to the family-owned equation generic.
The hot loop has no lookup registry. Later unchanged reproductions of an equation
receive no entry; distinct contributions require their own verified equations.

Analytical result details retain requested and effective identities for every formulation
slot, explicit modification flags, independent equivalent-earth selections and
orders, and normalized numerical options for internal surfaces, scalar material laws,
external cases and each required reduction case. Absence of a selected reduction remains `nothing` in this calculation records.

PSCAD extends the same equation generics with a `Val(:pscad)` execution payload:

```julia
earth_impedance(::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad})
earth_potential_coefficient(::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad})
internal_impedance(::Val{ID}, ::Val{Kind}, ::Val{:pscad})
```

These methods compile native settings. Every actual ordered pair is validated
using the Engine's physical geometry. A complete native settings record is used
for project export, execution, readback and numerical-input fingerprinting.
Unsupported potential selections fail before export; native `:default` potential
behavior belongs to PSCAD. Native conductor approximations likewise retain their
backend-owned `:default` rather than an alias to LCM's exact default.

PSCAD dispatch maps retained equations to native settings. Gary1976 maps to PSCAD's
`DERISEMLYEN` spelling; this creates no second mathematical registration. Carson1926
(overhead) and Pollaczek1926 (underground) map to native direct numerical integration.
These names follow PSCAD's [documented earth-return selections](https://www.pscad.com/webhelp-pscad-v5.1.0-ol/EMTDC/Transmission_Lines/Mutual_Impedance_with_Earth_Return.htm).
They identify the requested native controls, not a guarantee that native equations,
material assumptions or results equal LCM's implementations. The exported ground
permittivity remains the supplied material value.
PSCAD's `:default` selects that native setting, or native Lucca for a mixed arrangement.
Fixed backend calculations are recorded as such. PSCAD rejects analytical
hook overrides it cannot execute. FEM accepts only its four constitutive
selections and rejects analytical kernel keywords at construction. It executes
resolved material contributions without a second author registration. Constitutive
overrides passed to PSCAD remain subject to its documented export limits.

PSCAD export applies the selected temperature law to the same resolved conductor
materials used by the analytical engine. It evaluates each physical dielectric
layer before radial homogenization, explicitly enables native loss-tangent
handling, and retains the equivalent dielectric at its 50 Hz reference frequency.
The native loss-tangent cap is 10;
the aerial shunt setting uses the component's minimum, `1e-38 S/m`. These native
limits and the complete exported project accompany the results. Frequency-dependent
soil laws are rejected until their native parameter convention is verified; a
Julia constitutive callback is never converted into guessed Portela coefficients.

FEM batches reuse a field solve only when effective material, mesh and execution
inputs agree. Every request retains its metadata and independent result arrays.
Saved-run checks include inputs, implementation sources, executable identity and
artifact checksums. Incomplete compatible runs resume missing jobs; UI and explicit
remeshing requests execute separately. See [FEM](fem.md) for execution details.

`pipe_impedance=formula(:default)` uses the same selection grammar. Concentric
assemblies require no additional pipe term. Eccentric or multicore conducting
enclosures fail explicitly on the coaxial backend; FEM retains its supported
physical enclosure geometry.

## Cable-material temperature dependence

Select the constitutive law in the formulation and prescribe the operating
condition in the problem:

```julia
selected = Formulation(temperature_dependence=formula(:default))
problem = LineParametersProblem(system; temperature=80.0, frequencies=[50.0,1000.0],
    earth_props=homogeneous(rho=100.0))
result = compute(problem, selected)
```

The Materials-owned `TemperatureDependent` family evaluates electrical
resistivity. Its own `:default` implements
``\rho(T)=\rho_0[1+\alpha(T-T_0)]`` using each material's reference calibration.
`temperature_dependence=nothing` retains reference resistivity. The same slot
is available in `CableConstantsFormulation` and `LineCableModelsFEM`.
Temperature is prescribed in this electromagnetic calculation; no thermal
rating or temperature-field equation is implied.

Conductors consume the evaluated resistivity. Insulation/semicon constitutive
relations receive an ephemeral material with evaluated resistivity before their
electromagnetic equation; the stored reference material remains unchanged.
Original radial dielectric constituents are evaluated before aggregation.

A custom temperature contribution has the signature
`f(material, temperature, parameters, options, workspace) -> rho` in Ω·m.
Register its numerical defaults using the existing `FormulaMethod` grammar:

```julia
const TD = LineCableModels.Materials.TemperatureDependent
my_rho(m, t, parameters, options, workspace) = m.rho * exp((t-m.T0)/1000)
LineCableModels.computation_options(
    ::LineCableModels.FormulaMethod{:default,typeof(TD.temperature_resistivity)},
    ::typeof(my_rho)) = (;)
selected = Formulation(temperature_dependence=formula(:default;
    hooks=(contribution=my_rho,)))
```

A replacement owns its validity domain; all responses require positive real
resistivity, finite for conductors. The default approximation also enforces
``|T-T_0|<150`` K and a positive finite linear factor. The problem itself validates
finite temperature without imposing an unselected law. Numerical options remain
owned by the actual equation; the built-in linear law needs no integration.

## Finite formulation selection

The final formulation constructors participate in the same `Gridspace`
grammar as physical problem construction. Any line-parameter method slot or
the complete options tuple may be an explicit finite source:

```julia
formulations = Formulation(
    insulation_admittance = Grid((
        :Ametani2004,
        :default,
    )),
    earth_impedance = Grid((
        :Pollaczek1926,
        :Saad1996,
    )),
    combine = :product,
)
```

The result is `Gridspace{LineParametersFormulation}`. Every point is a
completely resolved scalar formulation: no `Grid`, symbol selector, or
`FormulaDefinition` reaches `compute`. `combine=:zip` aligns formulation fields and
broadcasts singleton fields. This local composition is separate from
`Combinatorial`, which always evaluates the Cartesian product of problem and
formulation points.

`CableConstantsFormulation`, `ModalTransformationFormulation`, and backend
formulation constructors follow the same rule. A deterministic `Grid` of
already completed, potentially external formulations is also accepted by
`Combinatorial`.

For each selected problem, the Coaxial collection dispatch validates and
lowers the physical declaration once. LineParameters flattens each design and
constructs `LocalCableData` plus geometry/index input once before creating a
separate workspace for every formulation. CableConstants performs its own
independent one-flatten orchestration. Formula-dependent mutable matrices,
earth/EquivalentHomogeneous values, reduction maps, and trace buffers remain workspace-local.
The generic collection dispatch simply invokes established scalar `compute`
methods and therefore supports external problem/formulation pairs without a
new registration layer.

## Earth-free cable constants

`CableConstantsProblem`, `CableConstantsFormulation`, and `CableConstants`
belong to Engine. They reuse the registered internal-impedance,
insulation-impedance, insulation-admittance, and semicon-admittance formulas,
and the same earth-free local primitive assemblers used by LineParameters.
Their public solve orchestration and reduction remain separate.
The default bundle is:

```julia
CableConstantsFormulation(
    internal_impedance = formula(:default),
    insulation_impedance = formula(:default),
    insulation_admittance = formula(:default),
    semicon_admittance = formula(:default),
)
```

`Engine.flatten(LineCableModelsCoaxial(), design)` supplies a
frequency-independent, unreduced `CableBlueprint`. Contiguous components
sharing one radial centre form one concentric assembly. Constitutive relations
are evaluated only after the workspace has been allocated. The Engine retains
each assembly's innermost terminal, grounds every additional outward terminal,
assembles and reduces the local N-terminal series-impedance matrix, and combines
the physical dielectric layers in radial series. A one-terminal assembly uses
the declared outer dielectric boundary directly; it does not require a metallic
sheath. Earth impedance, earth admittance, EquivalentHomogeneous, Γ, position, transposition,
and bundle reduction never enter this workflow.

`CableConstants(design; temperature=20, frequency=50)` is the convenience
entry point. CableConstants admits only the 50 Hz and 60 Hz datasheet base
frequencies. The result owns `cores`, aligned `R/L/C/G` vectors, and the
evaluation frequency. A conventional coaxial cable has one row and supports
`only(constants)`.

## Completed-result read side

`AbstractCoreResult` marks direct LineCableModels-owned computation results;
`CableConstants` and `LineParameters` are the current core result types.
`AbstractResultSpace{T}` marks completed finite collections of stored core
results. Its element type remains open so an external solver's concrete result
can be stored without subtyping a LineCableModels type. Result-space
constructors reject abstract element types and nested result-space envelopes.

Core results own their scientific extraction methods.
[`observe`](@ref) reads native numerical values through function-object
selectors:

```julia
observe(parameters, Z)
observe(parameters, L, 1, 1, Colon())
observe(parameters, Y, angle, 1, 1, Colon())
```

The public `Z`, `L`, and other laconic accessors delegate to these methods.
Consumers do not inspect `LineParameters` storage or repeat the R/X/L/G/B/C
formulae.

Direct numerical access remains available:

```julia
parameters.Z[1, 1, :]
@view parameters.Y[1, 1, :]
Z(parameters, 1, 2)
@observe parameters L[1, 2, :]
```

`@observe` expands the indexed expression to `observe(parameters, L,
1, 2, Colon())`. With no source argument it constructs the same plain request
tuple without reading a result:

```julia
request = @observe R[:, :, :]
magnitude_request = @observe (Z, abs)[:, :, :]
magnitude = @observe parameters (Z, abs)[1, 2, :]
```

The request tuple is an implementation representation, not a second selector
type. `quantity(Z, abs)` and `quantity(Z, angle)` remain the transformed
scientific identities.

[`observables`](@ref) publishes only explicitly requested values:

```julia
published = observables(
    parameters,
    (
        (frequencies, Colon()),
        @observe(R[1, 1, :]),
    );
    units = (
        LineCableModels.Units.units(:base, :hertz),
        LineCableModels.Units.units(
            :base,
            :ohm;
            per = (:kilo, :meter),
        ),
    ),
)
```

Every positional payload contains only `values`, `quantity`, and `unit`.
Publication converts and detaches `values`; it does not attach labels, result
objects, execution options, Gridspace points, or Monte Carlo context.

Line plotting accepts explicit observable requests. Its public convenience
forms expand selectors such as `Z`, `real`, and `angle` once, at the optional
Makie surface, then enter the same request path. The plotting extension groups completed
observations with the qualified `Units.family(::Quantity)` metadata. Series and
shunt identities return `Val(:series)` and `Val(:shunt)` respectively. Neither
the plotting extension nor ReportBuilder owns another quantity or family map.

The qualified `Grammar.validate_observables` method is the single request
validation used by direct publication and generic reports. It validates the source
declaration, request identities, and positional unit alignment.
`Grammar.unit_targets` resolves a tuple of requests to aligned `UnitExpr`
values. A unit override may be `nothing`, a metric-prefix `Symbol`, an explicit
`UnitExpr`, or a quantity-keyed collection used by an entry-point normalizer.
Line plots, Monte Carlo plots, and reports all use this path.

`LineCableModels.Units` owns `Unit`, `UnitExpr`, `Quantity`, `units`,
`quantity`, `native_unit`, `display_unit`, `scale_factor`, `label`, and
`symbol`. The plotting extension derives scientific axes from publication payloads.
ReportBuilder derives human-facing tables through `report`. Neither consumer
owns a quantity map, physical transform, or ordinary scientific unit string.

`Quantity{Q}` is a fieldless typed identity for extension methods and internal
publication payloads. Ordinary calls use scientific selector functions:

```julia
label(R)
symbol(Z, angle)
native_unit(R, :pul)
display_unit(Z, abs, :total)
```

The selector methods delegate through `quantity`; they do not contain a second
label or unit map. An external selector adds its identity and metadata at the
Units boundary:

```julia
function profile_response end

LineCableModels.Units.quantity(::typeof(profile_response)) =
    LineCableModels.Units.Quantity{:profile_response}()

LineCableModels.Units.native_unit(
    ::LineCableModels.Units.Quantity{:profile_response},
) = LineCableModels.Units.units(:base, :ohm)

LineCableModels.Units.display_unit(
    ::LineCableModels.Units.Quantity{:profile_response},
) = LineCableModels.Units.units(:milli, :ohm)

LineCableModels.Units.label(
    ::LineCableModels.Units.Quantity{:profile_response},
) = "Profile response"

LineCableModels.Units.symbol(
    ::LineCableModels.Units.Quantity{:profile_response},
) = "u"

label(profile_response)
display_unit(profile_response)
```

`Quantity`, `Unit`, and `UnitExpr` remain qualified extension vocabulary. The
package root exports only the six metadata functions used with scientific
selectors.

Higher-order results remain containers of owned products. `result`,
`statistics`, `samples`, and `histograms` select those products; they are not
zero-argument aliases for publication. UQ reads trial results with
`observe`, while its retained statistics, samples, and histograms implement
the same selector grammar for later publication. Monte Carlo run settings and
resolved point values are read through `root_seed`, `point_seed`,
`trial_count`, `confidence`, `cdf_tolerance`, and `sampling_distribution`.
Consumers do not inspect the result or its formulation fields.

## Coaxial workspace and supplemental output

`LineCableModelsCoaxial` solves concentric coaxial assemblies. A sector,
stranded, or otherwise nonconcentric part must be represented by the equivalent
round/concentric properties owned by DataModel before it reaches this backend.
The backend then owns frequency scans, self and mutual line parameters, earth
effects, and reduction; it does not redefine cable-design equivalence or modal
coordinates.

`LineParametersWorkspace` is the coaxial backend's per-computation working
state. Its constructor adapts a completed physical system once, evaluates
earth data, constructs cable and reduction indices, and allocates the matrices
and adaptive-quadrature segment storage used by the frequency loop. The
workspace is not a result, public data model, or alternate cable
representation. Its shunt solvers consume DataModel's ordered physical
dielectric layers directly, so the analysis is independent of the frequency
used when a lossy homogeneous export representation is requested.

The workspace separates four owned concerns:

- `input`: immutable numerical input derived from the problem;
- `invariants`: reusable physical values and index maps;
- `buffers`: mutable storage reused while solving every frequency;
- `capture`: optional diagnostic matrices allocated before the loop.

The ordinary result is always `LineParameters`. Requesting
`options=(trace=true,)` retains completed diagnostic arrays under
`details(parameters).trace`; it does not select another result type.

## Modal transformations

Modal decomposition is independent of the backend that produced fully coupled
phase-domain matrices. `LineCableModels.Transforms` owns its own problem,
formulation, registered formula files, and default backend:

```julia
phase = compute(line_problem, line_formulation)
modal = compute(
    ModalTransformationProblem(phase),
    ModalTransformationFormulation(
        formula(:default; options=(iteration=(convergence=1e-8,),)),
    ),
)
rebuilt = compute(ModalTransformationProblem(modal))
```

`LineCableModelsModal` is the default backend for this workflow. Each formula
has one route and returns a `ModalOperators` value containing the complete
frequency-dependent phase-to-modal voltage and current tensors. The shared
backend applies those operators to both `Z` and `Y`.

The modal `LineParameters` result carries `ModalDomain(operators, formula)` as
its domain value. The stored operators preserve the resolved mode
order, scaling, and complex phase convention and make the transformation
bidirectional without rerunning the decomposition. An operator-less modal
result cannot be constructed through the admitted `LineParameters` interface.
The formula value is retained through one formula-family storage parameter; its
specific author identity does not parameterize the domain. Modal results from
different registered routes therefore remain one concrete result-space element
type. Numerical inverse dispatch uses the concrete operator tensor and does not
inspect formula calculation records.

The retained modal formula is selected by `ModalTransformationFormulation()`.
Explicit controls use `formula(:default; options=(iteration=(convergence=1e-8,),))`.
A custom decomposition uses `hooks=(contribution=my_route,)` and returns
`ModalOperators` through the same application and inverse-transformation code.
The default tracks eigenpairs with Levenberg–Marquardt iteration, retaining a
matched conventional eigensolution when iteration fails. Its bibliography stays
in `src/transforms/formulas/default.jl`.

[`ComputationDetails`](@ref) is an alias for `NamedTuple`.
[`computation_details`](@ref) reads the fixed-key details tuple owned by a
registered formulation type. There is no general method: an unregistered
formulation raises `MethodError`. Higher-order calculations dispatch directly
on `typeof(formulation)` while collecting retained records; no owner registry
or wrapper token intervenes.

[`ParametricResult`](@ref), [`LinearErrorResult`](@ref), and
[`MonteCarloResult`](@ref) store the concrete details tuple type. Retention is
disabled by default, so `details(result) === (;)`. The higher-order formulation
owns the retention option:

```julia
Combinatorial(formulation; options=(retain_details=true,))
LinearError(formulation; options=(retain_details=true,))
MonteCarlo(formulation; trials=100, options=(retain_details=true,))
```

Parametric and linear calculations retain `(points=records,)`, with one record
per core result. Monte Carlo retains `trials`, `failures`, and
`failure_summary`, each aligned by Gridspace point. `trials` contains one inner
computation record per accepted trial. Each failure record contains the
attempt, target trial, failure stage, realised argument tuple, error type and
message, and a bounded stack summary. Statistics, samples, histograms, seeds,
and accepted-trial counts remain dedicated result fields.

## Formulation options

[`formulation_options`](@ref) validates values that alter the mathematical
calculation represented by a formulation. Dispatch uses the formulation owner
type rather than the public construction selector:

```julia
formulation_options(LineParametersFormulation, options)
```

The default line-parameter formulation owns:

- bundle and Kron reduction.
- ideal transposition.

The normalised named tuple is stored in `LineParametersFormulation.options`.
`PSCADFormulation` uses the shared physical options and currently requires
unreduced, untransposed matrices.

`Formulation()` constructs the default method bundle without a backend
tag. `LineParametersFormulation` owns the formulation options;
`LineCableModelsCoaxial` separately owns execution. Symbol and `Val` selectors
remain available for external backends, but there is no
`:line_cable_models` or legacy `:analytical` selector.

Modal formulas carry an `iteration` section containing convergence, iteration
count, damping and the `:matched` or `:error` fallback settings. The computation action
accepts `offdiagonal_tolerance` separately. Frequency continuation belongs to one
run; result details record the frequency indices where matched eigensolutions were
used. The stored voltage/current operators are retained for inverse transformation.

## Computation options

[`computation_options`](@ref) validates values belonging to one execution.
Formula numerical options select how the owning equation is evaluated. Backend
execution options govern output, tracing, logging and callbacks.

The coaxial backend accepts:

```julia
(
    verbosity = (default = 0,),
    output_basis = :pul,
    trace = false,
    on_result = nothing,
)
```

`trace=true` preallocates diagnostic capture with the workspace and attaches
the retained matrices to `details(result).trace` after computation.

Coaxial, FEM, and PSCAD computations accept an optional callable
`on_result(problem, index, result)`. It runs synchronously after each completed
formulation, including a reused result, before computing the next selection.
The index refers to the submitted formulation collection (`1` for a scalar
call). This allows manual campaigns to save completed results without waiting
for the whole batch. The callback must not mutate its arguments; its return
value is ignored and any exception stops execution. It is an execution option,
not a formula or a frequency-loop operation. Ordinary calls leave it as `nothing`.

The PSCAD backend accepts:

```julia
(
    output_stem = "case_name",
    remote = remote_config,
    verbosity = (default = 0, PSCAD = 2),
    output_basis = :pul,
    on_result = nothing,
)
```

`remote` must be a `PSCAD.RemoteConfig`. `output_stem` names files
created by that execution. Neither value belongs to `PSCADFormulation`.

Both option sets are ordinary `NamedTuple`s, aliased as
[`FormulationOptions`](@ref) and [`ComputationOptions`](@ref). Callers can
compose them with `merge`. Each owner rejects unknown keys and returns a
fixed-key normalised tuple. There is no general fallback and no conversion
from dictionaries, pairs, or `nothing`.

`MonteCarlo` owns a separate outer computation-option tuple. Its normalised
keys are `retain_details`, `on_error`, and `max_failures`. `on_error=:fail` is
the default and rethrows every exception. `on_error=:retry` requires
`retain_details=true` and rejects only `DomainError` realisations until the
requested accepted-trial count is reached or `max_failures` is exhausted.
Other exception types always propagate immediately.

## Gauntlet routing

`BenchmarkDefinition` coordinates two computations but does not own either backend's
keys. Its computation options form an outer tuple:

```julia
(
    output_basis = :pul,
    reference = (
        output_stem = "case_name",
        remote = remote_config,
        verbosity = (default = 0, PSCAD = 2),
    ),
    candidate = (
        verbosity = (default = 0,),
    ),
    benchmark = (
        samples = 10,
        seconds = 10.0,
    ),
)
```

The runner validates only this outer shape. The runner passes `reference` and
`candidate` to the corresponding `computation_options` methods and forces the same
`output_basis` into both so their results are comparable. Live and record runs
load a configured remote endpoint only when `reference.remote` is absent.
Snapshot runs never load remote configuration. Run mode, snapshot writing,
comparison tolerances, expected dimensions, and port ordering are not
computation options.

## Extending the engine

An external package may own a backend identity and a separate formulation
type. The backend's `compute` method normalises execution options before doing
work:

```julia
import LineCableModels:
    AbstractFormulation,
    AbstractCoreResult,
    ComputationOptions,
    FormulationOptions,
    ComputationDetails,
    computation_details,
    computation_options,
    compute,
    formulation_options

struct ExternalEngine end

struct ExternalFormulation{O <: NamedTuple} <: AbstractFormulation
    options::O
end

function formulation_options(
    ::Type{ExternalFormulation},
    options::NamedTuple,
)::FormulationOptions
    isempty(options) || throw(ArgumentError("unsupported formulation option"))
    return (;)
end

function computation_options(
    ::Type{ExternalEngine},
    options::NamedTuple,
)::ComputationOptions
    unknown = filter(key -> key != :tolerance, keys(options))
    isempty(unknown) || throw(ArgumentError("unsupported computation option"))
    normalized = merge((tolerance = 1.0e-8,), options)
    normalized.tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    return (tolerance = Float64(normalized.tolerance),)
end

function compute(
    ::ExternalEngine,
    problem,
    formulation::ExternalFormulation;
    options::NamedTuple = (;),
)
    execution = computation_options(ExternalEngine, options)
    # Use `problem`, `formulation`, and `execution` here.
end
```

An external implementation does not need a dedicated options struct or private
wrapper around the two normalisation functions. If it omits either Grammar
method, Julia raises `MethodError`.

The same backend may expose supplemental output without changing the generic
higher-order result types:

```julia
import LineCableModels: ComputationDetails, computation_details

struct ExternalResult <: AbstractCoreResult
    parameters
    diagnostics::NamedTuple
    raw::Dict{String,Any}
end

function computation_details(
    ::Type{<:ExternalFormulation},
    output::ExternalResult,
)::ComputationDetails
    return (
        diagnostics=output.diagnostics,
        raw=output.raw,
    )
end
```

An external backend that cannot modify or wrap its solver's concrete return
type may store that type directly in a result space. `AbstractCoreResult` marks
owned direct results; it is not an admission requirement for external result
payloads.

The outer keys and their types are fixed for `ExternalFormulation`. Dynamic
vendor channels remain inside the explicit `raw` leaf. ParametricBuilder and
UQ collect these records only when `retain_details=true`; they do not inspect
the fields.

## Reports and XLSX output

[`report`](@ref) executes `select`, `tabulate`, `illustrate`, `encode`, and
`write`, then constructs a `ReportArtifact`. `select` and `tabulate` are
required. Optional stages inherit the abstract-root no-op and in-memory reports
return `ReportArtifact.output === nothing`.

For formulation comparisons, ReportBuilder retains unformatted data in
`artifact.published` and exposes `summary`, `maxima`, `formulations`, `calculations`,
`comparisons` and `terms` through `artifact.table`:

```julia
using LineCableModels.ReportBuilder: BenchmarkTableDefinition

candidates = compute(problem, Formulation(earth_impedance=Grid((:default, :Pollaczek1926))))
reference = compute(problem, Formulation())
artifact = report(BenchmarkTableDefinition(), (; reference, candidate=candidates))
artifact.table.summary
artifact.table.terms
# After loading a Makie backend:
plot(artifact, (Z, Y))
```

The default request compares all Z/Y/R/L/G/C matrix terms in five bands. It creates
no figure. A retained publication is selected without recalculating RMS; new
numerical settings require explicit reanalysis. Native output terminal identities
must agree. Multiple problems require an explicit plot selection; a single
reference is overlaid once alongside all selected formulations. The
[Gauntlet guide](gauntlet.md) explains saved results, summaries and publication.

[`XLSXReportDefinition`](@ref) owns the human-facing line-parameter workbook:

```julia
using XLSX

artifact = report(
    XLSXReportDefinition(file_name="line_parameters.xlsx"),
    parameters,
)
artifact.output
```

ReportBuilder selects values through `observables`, builds one wide table with
coordinate columns followed by one column per observed quantity, and encodes a complete
[`LineCableModels.ReportBuilder.XLSXWorkbook`](@ref) containing the destination,
ordered sheet names, and final cell strings. Loading XLSX activates the package
extension that writes only this encoded description and records its path in
[`ReportArtifact`](@ref). Relative and default paths resolve from the caller's
current working directory; the package source tree is never the implicit
destination. `export_data(:xlsx, parameters; ...)` remains a thin ImportExport
convenience call that returns the same path. ImportExport owns no second
workbook implementation.


Scalar calculation selections are retained in `details(result).formulations`.
Its `requested` and `methods` fields hold complete requested and resolved records,
including physical parameters and explicit hooks. Formula identifiers are available
as `record.requested.earth_admittance.identifier` (or through the corresponding
`air`, `earth`, `mixed` leaf). The existing `effective`, `modified`, `numerical` and
`equivalent_earth` records retain the applied analytical interaction information.
Reports use these records directly; they do not infer a selection from numerical
agreement or collapse different voltage references into the same formula label.
