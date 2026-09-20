# Computational engine

LineCableModels separates the physical problem, selected equations, and numerical
execution. Source equations and bibliography belong in the implementing formula
file. Backend/formula methods select an implementation through Julia dispatch.

Every family owns a `:default` routing identifier and resolves it to an explicit
implementation before evaluation. Literature references remain attached to
the equations without determining their software names.

| Family | Registered choices |
|---|---|
| Internal impedance | `:default`, `:schelkunoff1934` |
| Insulation impedance | `:default`, `:ametani1980` |
| Pipe impedance | `:default`, `:none` |
| Insulation admittance, semicon admittance | `:default`, `:lossless`, `:lossy` |
| Local shunt geometry | `:default` (coaxial), `:coaxial`, `:boundary` |
| Earth impedance | `:default`, `:unified`, `:carson1926`, `:pollaczek1926`, `:gary1976`, `:wedepohl1973`, `:saad1996`, `:ametani2009`, `:lucca1994`, `:wise1934`, `:xue2018` |
| Earth admittance | `:default`, `:unified`, `:pollaczek1926`, `:wise1948`, `:xue2018` |
| Frequency-dependent soil properties | `:default`, `:constant`, `:alipio2014`, `:cigre2019`, `:datsios2019`, `:longmire1975`, `:messier1985`, `:portela1999`, `:scott1967`, `:visacro1987`, `:visacro2012` |
| Equivalent earth | `:default`, `:bottommost` |
| Modal transformation | `:default`, `:chrysochos2014` |
| Temperature-dependent resistivity | `:default`, `:linear` |

The internal default and `:schelkunoff1934` retain Schelkunoff's tubular
conductor expressions; insulation impedance's default and `:ametani1980`
retain the annular magnetic term documented by Ametani.
Both earth defaults route to `:unified`, which implements the supplied circumferentially averaged framework
with complete enclosed-current normalization. Only `:default` and `:unified`
execute in the two built-in coaxial earth families. Other registered identities
throw "not yet implemented" when requested; they never substitute Unified.
Independent external-backend and consumer-defined implementations are unaffected.
Dielectric `:default` selections
route to explicit `:lossless` equations; `:lossy` retains conductivity and the
material's supplied polarization losses. The FrequencyDependent `:default`
routes to `:constant`, which preserves static properties. Its explicit
literature relations model measured soil dispersion; their references remain
attached to the equations without determining their software names.
EquivalentHomogeneous selects the basement when explicitly requested, and the
modal default performs Levenberg–Marquardt tracking. The EquivalentHomogeneous
default is a package-owned policy rather than an author equation.

## Coaxial computation

The workflow is design, system, problem, formulation, then `compute`. Initialization
flattens the selected designs into blueprint tables, binds equations to their
required inputs and selected output entries, and allocates one
`LineParametersWorkspace`.

The calculation evaluates fixed-temperature conductor properties and required
frequency-dependent earth tables, then follows this order at each frequency:

1. Complete equivalent-earth evaluation and dielectric material admittivities.
2. Calculate cable-local impedance and potential coefficients.
3. Calculate the selected exterior impedance and potential coefficients together
   through `earth!`.
4. Assemble primitive matrices, apply the selected reductions and store Z/Y.

`BeforeFD` and `AfterFD` retain their distinct ordering. Material-law calls receive
valid material objects. Field formulas receive completed material properties;
their wave numbers and mathematical approximations do not replace those values.
Boundary-shunt coefficients remain blueprint-time quantities.

The workspace input/invariants retain geometric tables and equation/input/output
bindings. Its buffers retain evaluated material tables, local coefficients,
separate exterior Z/P destinations, primitive/reduction matrices, QuadGK storage
and any numerical scratch arrays required by the selected equations. Trace
storage is optional. Independent computations share no mutable calculation
buffers.

Ordinary earth equations evaluate their indexed cases. A coupled equation may
require full-system inputs even when only some of its entries are selected.
Those inputs and compatible Z/P consumers are bound during initialization.
Compatible consumers calculate one response per frequency; different
configurations calculate separately and publish their selected entries before
reusing scratch. The binding records the fixed correspondence, not the validity
of a previously calculated response. Every call calculates anew: there is no
subordinate workspace or readiness/reset lifecycle.

## Internal shunt geometry

`shunt_model` selects the cable-local geometry approximation independently of
the dielectric material laws. Both `Formulation()` and
`CableConstantsFormulation()` default to `shunt_model=:default`: ordinary
coaxial annuli, with no boundary solve. `:coaxial` selects this explicitly.

`insulation_admittance` and `semicon_admittance` still select material
admittivity κ [S/m]. They do not select geometry. For one homogeneous annulus,
the shunt branch is `y = 2πκ/log(ro/ri)` [S/m]; successive dielectric layers
combine in series. This is the same reduced geometry used by the coaxial
backend, including its equivalent wire-screen geometry.

Select `shunt_model=:boundary` to resolve eligible open wire and finite-tape
groups inside a closed circular shield before the frequency sweep. This is an
explicit local refinement of the coaxial backend, not another earth formula.
It supplies a coupled terminal operator, not a fitted scalar permittivity:
an annular chain cannot retain direct coupling across an open intermediate
screen. The same blueprint construction supplies `CableConstants`.

The resolved local calculation preserves physical filler permittivity, finite
tape thickness, conductor terminal membership and concentric dielectric layers.
It couples all intermediate terminals in a qualifying domain together. Round
wires use auxiliary sources; tapes use integrated corner-weighted charges on
their two circular faces and two ends. The inner core retains the existing
equivalent-core assumption, including bounded circular, rectangular, compacted
and sector stranded constructions. Strand packing and clearance are unchanged.
Series impedance, earth return, terminal ordering and matrix reductions are
also unchanged. Only the applicable internal shunt contribution changes.
Existing coaxial lowering requirements, including radial terminal ordering,
remain in force; this does not add support for declarations that lowering
already rejects.

Qualification currently requires an explicitly filled annular host, continuous
concentric dielectric paths, exposed whole wire/tape boundaries, and a closed
circular reference shield. Separate shielded assemblies are handled in their
own local frames. Overlapping same-terminal faces, interacting courses in
different hosts, nonconcentric dielectric interfaces and noncircular reference
shields retain the existing equivalent-coaxial treatment. They are not fed to
an inapplicable circular Green function. The resolved path currently uses the
built-in `:lossless` insulation/semicon laws or their `:default` aliases.
Lossy or custom constitutive selections are unsupported by `:boundary` and raise `BoundarySolveError` for
eligible domains. Select `:coaxial` to retain their radial calculation without
losing conductivity or frequency dependence.

Inspect `details(result).data.shunt_model` for requested/effective model, domain
terminal ranges, solve counts, and diagnostics. Boundary-resolved domains
coexist with ordinary radial intervals outside their coverage. `effective`
describes the qualifying domains; `:mixed` indicates explicit fallback in some
of them. Numerical grid residuals and small-coupling indicators are not
certified error bounds. Whole-matrix convergence can conceal substantial
relative changes in weak individual couplings.

Formulation-aware blueprint construction uses bounded dense storage and
in-place pivoted QR. The completed `CableBlueprint` owns lossless terminal
capacitance and potential coefficients, local terminal coverage and numerical
outcomes. Workspaces consume these coefficients; they perform no boundary solve.
Repeated equivalent domains share coefficient matrices within the same
construction call. Formulations with identical local selections share blueprints
even when earth-return choices differ. Independent uncertainty sources prevent
sharing, and dense boundary matrices are released after construction.

Select the formulation and call `compute` directly:

```julia
formulation = Formulation(shunt_model=:boundary)
result = @time compute(problem, formulation)
```

There is no separate preparation object or execution option. Default coaxial
blueprints contain no boundary blocks, extract no boundary domains, and evaluate
no boundary material law. Their reports describe the concentric assembly ranges;
no boundary audit is implied.
Boundary coefficients are independent of the frequency grid and external earth
model. They are reused across the frequency sweep, but a fresh `compute` call
constructs fresh blueprints: there is no process-global or cross-call cache.
Changed geometry and Monte Carlo realizations therefore receive new coefficients.
Do not mutate coefficient arrays shared by completed blueprints.

Production reports unmet quadrature targets, estimated rank, reciprocity and
positive-capacitance concerns as warnings. The computed finite result is retained
without retry, repair or substitution. Dense-storage limits, degenerate columns,
actual solve failure and nonfinite physical results remain errors. The independent
validation grid and derivative step-refinement are opt-in:

```julia
boundary = formula(:boundary;
    parameters=(fallback=:error,),
    options=(
        resolution=(wire=64, order=32, quadrature=256, modes=1024),
        integration=(rtol=1e-8, atol=1e-10, maxevals=100_000),
        audit=false,
    ))
formulation = Formulation(shunt_model=boundary)
```

Set `audit=true` for the independent checks. Unaudited residual fields are
`nothing`, not zero. The audit warns on unmet criteria and does not alter the
computed terminal matrix.
The integration tolerances control dimensionless logarithmic moments, not
the error in an individual terminal coupling.

See [`BoundarySolveError`](@ref) and [`LineCableModels.Engine.ShuntModel.Formula`](@ref)
for the failure and selection contracts.

Actual failures propagate by default. An explicit `parameters=(fallback=:coaxial,)`
permits annular replacement only after recognized numerical or unsupported-law
failures, with a warning and recorded reason. Finite-result quality warnings do
not trigger that fallback. Invalid inputs and unexpected
exceptions propagate. UQ wrappers reject automatic fallback; choose strict
`:boundary` or `:coaxial` for the whole study. Monte Carlo also checks that
shunt-model coverage stays fixed across realizations and never retries a
`BoundarySolveError` as a geometry rejection.

With Measurements loaded, local sensitivities preserve the original correlated
inputs. The nominal QR is reused in an implicit least-squares derivative that
includes its residual term. Centered kernel derivatives use one step in
production; the audit additionally checks step halving. Derivative matrices
are streamed in bounded blocks rather than stored
as dense Measurement arrays. This is fixed-topology linear propagation, not a
claim about differentiability across a contact or strand-count transition.
Each Monte Carlo realization prepares its own physical geometry/material
operator; identical cables/formulations within that realization share it.
Concentricity checks include coordinate dependencies: two centres with equal
nominal positions but independent uncertainty do not qualify as concentric.
The local boundary calculation uses Float64 workspaces; it does not promise
arbitrary-precision boundary accuracy when surrounding scalar inputs use a
wider type.

### Numerical method and attribution

The annular Green function is assembled from classical cylindrical Laplace
harmonics: `r^m` and `r^-m`, with `1` and `log(r)` for the zero mode.
[Schelkunoff1934](@cite), Eqs. (122)–(124), p. 573, gives this radial basis
in its treatment of cylindrical fields in the small-radius-to-wavelength
limit. Here the same basis solves the electrostatic potential problem.
[Sunde1968](@cite), Section 1.5, Eq. (1.20), p. 11, supports the coaxial
capacitance normalization; Section 1.6, Eqs. (1.39)–(1.42), p. 15, gives
the load/reflection transformation applied here in log-radius coordinates.
Matching the dielectric interfaces and grounded inner/outer boundaries
assembles these ingredients into the implemented layered annular kernel.

The finite-thickness, layered-Green-function, whole-face charge approach has
precedent in [Bernal1997](@cite), especially Section IV, Eq. (10). Their
Maxwell-weighted Chebyshev/Galerkin formulation is not reproduced verbatim:
this implementation uses a concentric circular Green function, normalized
Jacobi face charges and oversampled boundary collocation with pivoted QR.

[Campione2018](@cite), Section 2, Eq. (3), provides a cable-screen application
of dielectric image factors `(εhost − εadjacent)/(εhost + εadjacent)`, also
used in the extracted nearest-interface terms here. Their cylindrical braid
is locally approximated by a plane; that work is not a derivation of this
implementation's concentric annular kernel or radial layer recursion.

The corner weights encode the finite-energy edge behaviour of
[Meixner1972](@cite). Corners on dielectric interfaces use the actual adjacent
permittivities, following the metal–dielectric wedge principle of
[VanBladel1985](@cite), rather than assigning every corner the homogeneous
exponent. [Classen2011](@cite), Sections 2 and 4, supports incorporating known
singularities into an approximation space; its FIT/DG methods are not this
tape element. The implemented corner equation and normalized charge measure
are documented below.

The need for special evaluation near source boundaries is discussed by
[HelsingOjala2008](@cite), Section 1. Here logarithmic direct/image moments use
adaptive Gauss–Kronrod quadrature after a cosine change of variable, while
the smooth remainder uses Gauss–Jacobi quadrature. This is a different
algorithm from their rational-quadrature scheme.

The following internal numerical methods document the equations, charge
normalization and quadrature needed to reproduce the calculation. They are
implementation details, not additional public modeling APIs. Source placement,
resolution controls, validation thresholds and uncertainty differentiation
remain package-specific choices, not accuracy guarantees supplied by these
references.

```@docs
LineCableModels.Engine.ShuntModel._shunt_load
LineCableModels.Engine.ShuntModel._shunt_kernel_coefficients
LineCableModels.Engine.ShuntModel._shunt_junction_exponent
LineCableModels.Engine.ShuntModel._shunt_tape_faces
LineCableModels.Engine.ShuntModel._shunt_log_moments
LineCableModels.Engine.ShuntModel._shunt_capacitance
```

## Formula selection

```julia
selected = Formulation(
    earth_impedance=formula(:default;
        options=(integration=(method=:quad, options=(rtol=1e-8,)),)),
    earth_admittance=formula(:default),
    earth_properties=formula(:default),
)
result = compute(problem, selected; options=(trace=true,))
```

`FormulaDefinition` carries an identifier, physical `parameters`, numerical
`options`, and an optional formula-local `equivalent_earth`. It is a passive
request. A family constructor resolves a symbol or declaration to a concrete
selection; a completed user-owned selection passes through unchanged. Built-in
catalogues are explicit inventories, not admission registries for user code.
There is no `hooks` field or callable-override channel.

A `FormulaMethod(selection, operation, Val(...), ...)` binds the actual selected
type to the family operation. The first argument of every equation is that
selection—not its identity symbol. `formula_id` is descriptive metadata; it
cannot redirect execution to a different implementation. `:default` only
routes, never implements an equation and never catches a failed calculation.

For the analytical families, the default routes are:

| Family | Explicit implementation |
|---|---|
| Internal impedance | `:schelkunoff1934` |
| Insulation impedance | `:ametani1980` |
| Earth impedance and potential coefficients | `:unified` |
| Insulation and semicon admittivity | `:lossless` |
| Frequency-dependent earth properties | `:constant` |
| Equivalent homogeneous earth | `:bottommost` |
| Temperature-dependent resistivity | `:linear` |
| Modal transformation | `:chrysochos2014` |
| Local shunt geometry | `:coaxial` |
| Pipe contribution | `:none` |

These routes do not expand applicability: `:none` does not implement a pipe
solver, and `:unified` still requires its supported earth geometry. PSCAD owns
separate native selections; FEM owns its field equations and accepts supported
material selections, not analytical impedance equations.

Numerical defaults belong to `formulation_options(::FormulaMethod{<:MySelection,
typeof(operation), ...})`. Empty defaults admit no controls. Unknown or unused
numerical sections are errors. The selected formula provisions QuadGK storage
through `initialize_buffers` when its required indexed consumers need integration;
option-key presence is not an allocation policy. Full-system earth physics also
belongs to the selected physical formulation.

Custom types implement the existing family operation and expose `parameters`
and `options` records. Internal selections expose per-surface options;
earth selections expose `assumptions` and
`equivalent_earth`. Constructors own parameter checking and option normalization.
Implement `formula_id`, `description`, `NamedTuple`, and `formulation_options`
for scientific inspection and persistence. A saved declaration is not executable
code; unknown saved leaf identities remain passive identities.

Unified accepts a prescribed longitudinal argument [1/m] through
`formula(:unified; options=(Γ=value,))`. A finite scalar applies to every
sample; a finite vector aligns one-to-one with the problem's frequency vector,
without sorting or interpolation. Zero is the default. This formulation option is not a
problem field or a common earth-formula requirement. No propagation-constant
root solver is implied. Material constitutive
laws evaluate valid material objects before field calculations. Selected earth
equations consume those properties and calculate their own wave numbers and
field approximations without redefining the material values.

The implemented equations use `s=jω` and the positive-time `exp(jωt)` phasor
convention. For the same real field expressed as `F̂₋ exp(Γ₋ x-jωt)`, the
positive-time representation has `F̂₊=conj(F̂₋)` and `Γ₊=conj(Γ₋)`; complex
material coefficients and response phasors must use that convention consistently.
The API accepts Γ verbatim and does not silently perform this conversion.
Γ is a prescribed longitudinal wavenumber, not an independent UQ input or a
request for modal iterations. Ordinary uncertainty in physical materials retains
its existing treatment.

`EarthPair` carries conductor row/column indices, integer source/target layer indices,
heights, horizontal separation and an explicit self radius. Self means the same
conductor; distinct conductors in one layer remain mutuals. A self pair has zero
horizontal separation, with its radius supplied separately. The geometry substitution
needed by a published self expression occurs at equation evaluation.

The external equation signatures are:

```julia
earth_impedance(selection::MyEarthImpedance, ::Val{Kind}, ::Val{S}, ::Val{T}, functor, pair, workspace)
earth_potential_coefficient(selection::MyEarthPotential, ::Val{Kind}, ::Val{S}, ::Val{T}, functor, pair, workspace)
```

`Kind` is `:self` or `:mutual`; source `S` is the matrix column, and target `T`
is the row. Layer 1 is air; soils occupy layers 2 through N. Required methods
are invoked through indexed dispatch; unsupported cases reach the throwing
fallback. There is no reflected equation-coverage preflight. Domain-defining
methods accept the three runtime payloads without extra
subtype constraints. Numerical specializations can optimize an admitted case.
Julia method availability determines which equations can be selected.

Each earth slot also accepts an explicit NamedTuple recipe:

```julia
selected = Formulation(earth_impedance = (
    earth = formula(:unified),
))
```

The same syntax applies independently to `earth_admittance`. The labels resolve
`(1,1)`, `(2,2)` and the two cross-layer directions before the existing indexed
method dispatch. Each case retains its own formula, numerical options and
material inputs. There is no combined formula identity. The example covers an
all-buried problem, not an overhead or mixed problem. Whole-family omission or
`nothing` routes to `:default`; an explicit recipe receives no implicit completion.
An unused leaf does not allocate numerical buffers, supply missing cases or
execute a kernel. True layered inputs require a scalar
selection; scalar multilayer and explicit EHEM behavior are unchanged.

The default potential formulation supplies both mixed directions. Selecting a
default leaf for only part of a matrix still assembles its complete auxiliary
system and then selects the requested final entries. Other selected output
equations never become inputs to its K/H kernels.

The workspace binds every required ordered pair before frequency evaluation.
Geometry and layer indices follow the same order. Both directions are evaluated;
assembly, reduction and modal transformation preserve the returned ordered entries.
No implicit reciprocity operation supplies a missing equation or averages its result.

The medium inventory is a separate physical restriction. A homogeneous formula
consumes exactly air and one soil half-space. A finite-layer model consumes its
whole declared inventory and interfaces. Explicit `Val(S), Val(T)` methods describe
its cases; arbitrary-layer Green-function generation remains deferred. Buried
placement in a vertical multilayer earth is rejected because its physical layer
indexing has no defined origin in the present geometry definition.

Only `:unified` and explicit `:default → :unified` routing are executable built-in
coaxial earth implementations, in both impedance and potential families. Other
author identifiers remain registered scientific identities with explicit
"not yet implemented" methods for their indexed cases. They never fall back to
Unified. Independent external-backend and consumer-defined methods are unaffected.
The defaults supply air/air, earth/earth and both mixed directions with
independent medium permeabilities. Every circumference must lie wholly in its
half-space and exterior circles must not overlap. The default requires homogeneous
earth or an explicitly globally consistent equivalent earth. Indexed dispatch
for arbitrary layer numbers is an extension contract, not an implementation of
new multilayer Unified equations.

The complete earth source kernels satisfy

```math
K=\mathcal Z-\Gamma^2\mathcal P_\phi/s,\qquad
L=A_r^{-1}-F_rK,\qquad
P_eL=H,\quad Z_eL=K+\Gamma^2H/s,\quad Y_eH=sL.
```

Here `s=jω`, Ze is in Ω/m, Pe in m/F and Ye in S/m. Matrices are solved on the
right, with source columns and receiver rows. Pe and Ye retain the direction of
the voltage path; they are not symmetrized. Air receivers use the interface
voltage; earth receivers use deep-earth voltage. These references are fixed by
the receiving layer, not selectable parameters.

Each formula owns its equations. Unified assembles its
complete current system once per frequency and shares compatible impedance and
potential calculations. Its indexed methods return physical self/mutual entries
from that coupled system; an isolated-pair call cannot supply it. Engine uses
the same explicit material and earth-calculation stages for every formula and does not know these
physical equations. Internal conductor and insulation
contributions retain the existing cable composition and terminal reductions.
With `options=(trace=true,)`, `trace.Zg` and `trace.Pg` expose the exterior
matrices; the returned total line admittance also includes insulation effects.
`trace.integrals` retains the native integral values and QuadGK error estimates,
with formula, frequency, receiver/source and term identifiers. These estimates
are not final-matrix error bounds. Trace-off computations retain no such history.

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

The default selects the deepest soil layer. Martins-Britto et al.
[Martins-BrittoLopes2020](@cite) found that deep-layer conductivity predominated
in magnetic ground-return impedance for the multilayer soil cases they studied.
This provides a qualified rationale for the default resistivity. Accuracy depends
on layer contrasts and frequency; selecting one layer does not implement the
paper's equivalent-conductivity formula or establish the accuracy of the selected
permittivity and permeability. Choose another built-in rule or a user-owned subtype of
`EquivalentHomogeneous.AbstractRule` implementing `equivalent_material`.

`:after` applies the selected frequency law to physical layers first. `:before` reduces static
properties and applies that same law to the resulting material. Physical and
effective pairs remain distinct in the binding, and reductions run for each
ordered interaction on which they depend. Layerwise evaluated properties are reused
between consumers when needed; the air material remains static.

For example, the numerical declaration for a user-owned outer surface equation is:

```julia
using LineCableModels: FormulaMethod, formulation_options, FormulationOptions
const II = LineCableModels.Engine.InternalImpedance
# MyConductor is a concrete InternalImpedanceFormulation owned by the user.
formulation_options(
    ::FormulaMethod{<:MyConductor,typeof(II.internal_impedance),Tuple{Val{:outer}}},
) = FormulationOptions()
```

The equation is `II.internal_impedance(selected::MyConductor, ::Val{:outer},
functor, workspace)`. Its state constructor prepares one `II.Functor` per
conductor and frequency. The same state is shared by all surfaces using that
selection. An integral equation declares its own integration section; an
algebraic equation does not inherit another formula's numerical options.

`InternalImpedance.surface_impedances(resolved_formula, r_in, r_ex, rho, mu_r, jω)`
returns `(inner,outer,transfer)` coefficients in Ω/m, with per-kind
numerical options applied. Internal kinds have no earth-layer selectors. Assemblers
own the current-basis transformation and matrix placement. The deferred pipe
contribution concerns one contained metal and its enclosing pipe; recursive
assembly and pipe equations are outside this implementation.

Like the earth selections, internal impedance accepts either one formula or a
complete named selection:

```julia
selected = Formulation(internal_impedance=(
    inner=formula(:default),
    outer=formula(:default),
    transfer=formula(:default),
))
```

Each surface can select a built-in or user-owned formulation that implements that
surface, with its own parameters and numerical options. A custom transfer
selection belongs directly in that leaf, as `transfer=my_transfer_model`.
Identical complete
selections share their prepared conductor state. Scalar shorthand retains its
existing numerical behavior. The internal term is called `transfer`; earth
`self`/`mutual` interaction names are unchanged.

Equation numerical sections belong to `FormulationOptions`. The existing
`formulation_options` methods validate and normalize them at the owning stage:

Spectral integration uses adaptive Gauss–Kronrod quadrature (`method=:quad`)
over the full spectral interval. Its controls are:

```julia
integration = (method = :quad, options = (rtol = 1e-8, atol = 0.0, maxevals = 10^7))
```

`SpectralIntegral(f)` stores only the complete callable integrated over
`[0, Inf)`. A formula includes its weights, any admissible contour and the
coordinate Jacobian in `f`, and supplies numerical subdivision points to
`integrate`. Engine knows none of those physical choices. Its half-line map
has no finite cutoff. Algebraic cases make no integration call and declare no
integration controls; they still receive the same computation workspace.

QuadGK retains physical scalar types, including BigFloat and correlated
Measurements inputs. Numerical error norms and sampling coordinates are nominal;
physical values and derivatives are not replaced by nominal values.

Quadrature returns its native `(value, estimated_error)`. An unmet requested
target produces a warning and returns the finite value without an outer retry,
tolerance tightening or matrix-error rejection. Invalid inputs, nonfinite values
and singular physical solves remain errors. The human judges accuracy; meaningful
matrix-accuracy assertions belong in tests.

Selected numerical formulas within one calculation share reusable segment storage.
The main workspace holds formula-owned subdivision and coupled-response arrays.
`initialize_buffers` is the single allocation extension point for both earth and
local formulas; numerical options do not determine allocation by key inspection.
Each formula supplies the complete callable and its own numeric subdivision hints.
Engine contains no physical feature finder or spectral sampler. Compatible
unified consumers share one current calculation per frequency and publish the
completed entries directly; independent calculations never share mutable buffers.

Each family explicitly includes its supported formula files, each returning one
identifier. `FormulaMethod` binds a selected formulation and its indexed case to
the family-owned equation generic.
The hot loop has no lookup registry. Later unchanged reproductions of an equation
receive no entry; distinct contributions require their own verified equations.

Result details retain one requested/resolved formulation record. Its `requested`
field preserves the declaration; `methods` records the selected identities,
physical parameters and owner-held numerical controls, including each explicit
equivalent-earth choice and order. Unspecified numerical defaults remain owned by
the selected equation rather than a second retained inventory of indexed calls.
Absence of an equivalent-earth selection remains `nothing`.

PSCAD extends the same equation generics with a `Val(:pscad)` execution payload:

```julia
earth_impedance(selection::PSCAD.NativeFormula, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad})
earth_potential_coefficient(selection::PSCAD.NativeFormula, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad})
internal_impedance(selection::PSCAD.NativeFormula, ::Val{Kind}, ::Val{:pscad})
```

These methods compile native settings. Every actual ordered pair is validated
using the Engine's physical geometry. A complete native settings record is used
for project export, execution, readback and numerical-input fingerprinting.
Unsupported potential selections fail before export. PSCAD defaults resolve to
`:direct_lucca` for earth impedance, `:coupled` for potential coefficients, and
`:cable_coax` for conductor and insulation impedance. These backend-owned
selections do not pretend to execute the analytical equations.

PSCAD dispatch maps retained equations to native settings. Gary1976 maps to PSCAD's
`DERISEMLYEN` spelling; this creates no second mathematical registration. Carson1926
(overhead) and Pollaczek1926 (underground) map to native direct numerical integration.
These names follow PSCAD's [documented earth-return selections](https://www.pscad.com/webhelp-pscad-v5.1.0-ol/EMTDC/Transmission_Lines/Mutual_Impedance_with_Earth_Return.htm).
They identify the requested native controls, not a guarantee that native equations,
material assumptions or results equal LCM's implementations. The exported ground
permittivity remains the supplied material value.
PSCAD's `:default` selects that native setting, or native Lucca for a mixed arrangement.
Fixed backend calculations are recorded as such. PSCAD rejects analytical
kernel selections or numerical controls it cannot execute. FEM accepts only its four constitutive
selections and rejects analytical kernel keywords at construction. It executes
resolved material contributions without a second author registration. Constitutive
selections passed to PSCAD remain subject to its documented export limits.

PSCAD export applies the selected temperature law to the same resolved conductor
materials used by the analytical engine. It evaluates each physical dielectric
layer before radial homogenization, explicitly enables native loss-tangent
handling, and retains the equivalent dielectric at its 50 Hz reference frequency.
The native loss-tangent cap is 10;
the aerial shunt setting uses the component's minimum, `1e-38 S/m`. These native
limits and the complete exported project accompany the results. Frequency-dependent
soil laws are rejected until their native parameter convention is verified; a
Julia material law is never converted into guessed Portela coefficients.

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
resistivity. Its `:default` routes to `:linear`, which implements
``\rho(T)=\rho_0[1+\alpha(T-T_0)]`` using each material's reference calibration.
`temperature_dependence=nothing` retains reference resistivity. The same slot
is available in `CableConstantsFormulation` and `LineCableModelsFEM`.
Temperature is prescribed in this electromagnetic calculation; no thermal
rating or temperature-field equation is implied.

Conductors consume the evaluated resistivity. Insulation/semicon constitutive
relations receive an ephemeral material with evaluated resistivity before their
electromagnetic equation; the stored reference material remains unchanged.
Original radial dielectric constituents are evaluated before aggregation.

A custom temperature law is a concrete selection, with no function-valued field:

```julia
const TD = LineCableModels.Materials.TemperatureDependent
struct ExponentialResistivity{P,O<:FormulationOptions} <: TD.TemperatureDependentFormulation
    parameters::P
    options::O
end
function ExponentialResistivity(; scale=1000.0)
    isfinite(scale) && scale > 0 || throw(ArgumentError("scale must be positive [K]"))
    ExponentialResistivity((scale=scale,), FormulationOptions())
end
TD.temperature_resistivity(::ExponentialResistivity, m, t, p, o, workspace) =
    m.rho * exp((t-m.T0)/p.scale)
LineCableModels.formulation_options(
    ::LineCableModels.FormulaMethod{<:ExponentialResistivity,typeof(TD.temperature_resistivity)}) = FormulationOptions()
LineCableModels.formula_id(::ExponentialResistivity) = :exponential_resistivity
LineCableModels.description(::ExponentialResistivity; compact=false) = "Exponential resistivity"
Base.NamedTuple(law::ExponentialResistivity) =
    (identifier=formula_id(law), parameters=law.parameters, options=law.options.data)
LineCableModels.formulation_options(law::ExponentialResistivity) =
    law.options
selected = Formulation(temperature_dependence=ExponentialResistivity())
```

The law owns its validity domain; all responses require positive real
resistivity [Ω·m], finite for conductors. The built-in linear approximation also
enforces `|T-T₀| < 150` K and a positive finite linear factor. The problem itself
validates finite temperature without imposing an unselected law.

Scalar equation suffixes are uniform: physical inputs, `parameters`,
`options`, then `workspace`. Insulation and semicon laws return finite
admittivity [S/m], normalized to `Complex{T}` for the input scalar type `T`.
A response requiring a wider scalar type is rejected; precision and measurement
uncertainty are never silently discarded. Soil laws return `EarthMaterial`;
their fitted coefficients are checked at construction, before a frequency sweep.

## Finite formulation selection

The final formulation constructors participate in the same `Gridspace`
grammar as physical problem construction. Any line-parameter method slot or
the complete options tuple may be an explicit finite source:

```julia
formulations = Formulation(
    insulation_admittance = Grid((
        :lossy,
        :default,
    )),
    earth_properties = Grid((
        :constant,
        :longmire1975,
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
calculation sequence, also flattening each design once. Formula-dependent mutable matrices,
earth/EquivalentHomogeneous values, reduction maps, and trace buffers remain workspace-local.
The generic collection dispatch simply invokes established scalar `compute`
methods and therefore supports external problem/formulation pairs without a
new registration layer.

## Earth-free cable constants

`CableConstantsProblem`, `CableConstantsFormulation`, and `CableConstants`
belong to Engine. They reuse the registered internal-impedance,
insulation-impedance, insulation-admittance, and semicon-admittance formulas,
and the same earth-free local primitive assemblers used by LineParameters.
Each calculation has its own solve and reduction methods.
The default bundle is:

```julia
CableConstantsFormulation(
    internal_impedance = formula(:default),
    insulation_impedance = formula(:default),
    insulation_admittance = formula(:default),
    semicon_admittance = formula(:default),
)
```

`Engine.flatten(LineCableModelsCoaxial(), design, formulation)` supplies a
frequency-independent, unreduced `CableBlueprint`; omitting the formulation
uses the default annular model. Contiguous components sharing one radial center
form one concentric assembly. Explicit boundary shunt coefficients are completed
during flattening. Frequency-dependent constitutive evaluation and conductor
temperature corrections remain in the calculation. The Engine retains
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
Makie API, then use the same observation requests. The plotting extension groups completed
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
state. Its constructor adapts a completed physical system once, binds equations,
constructs cable and reduction indices, and allocates numerical buffers. Material
evaluation is an explicit calculation stage, not a constructor side effect.
QuadGK arrays remain empty when no selected equation requires integration. The
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
`details(parameters).data.trace`; it does not select another result type.

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
A custom decomposition is a completed formulation implementing
`Transforms.modal_operators(selected::MyModalModel, parameters, model_parameters,
options, workspace)`. It returns `ModalOperators` through the same application
and inverse-transformation code. `ModalTransformationFormulation` retains its
requested declaration and resolved selection for the common inspection protocol.
The default tracks eigenpairs with Levenberg–Marquardt iteration, retaining a
matched conventional eigensolution when iteration fails. Its bibliography stays
in `src/transforms/formulas/chrysochos2014.jl`.

[`ComputationDetails`](@ref) is an immutable, nominal record parameterized by its
named-tuple payload. Access its fields through `.data`; it is not a tuple and
cannot be used as either kind of options.
[`computation_details`](@ref) returns the fixed-key details record owned by a
registered formulation type. There is no general method: an unregistered
formulation raises `MethodError`. Higher-order calculations dispatch directly
on `typeof(formulation)` while collecting retained records; no owner registry
or wrapper token intervenes.

[`ParametricResult`](@ref), [`LinearErrorResult`](@ref), and
[`MonteCarloResult`](@ref) store the concrete details record type. Retention is
disabled by default, so `details(result) == ComputationDetails()`. The higher-order formulation
owns the retention option:

```julia
Combinatorial(formulation; options=(retain_details=true,))
LinearError(formulation; options=(retain_details=true,))
MonteCarlo(formulation; trials=100, options=(retain_details=true,))
```

Parametric and linear calculations retain `ComputationDetails(points=records)`, with one record
per core result. Monte Carlo retains `trials`, `failures`, and
`failure_summary`, each aligned by Gridspace point. `trials` contains one inner
computation record per accepted trial. Each failure record contains the
attempt, target trial, failure stage, realized argument tuple, error type and
message, and a bounded stack summary. Statistics, samples, histograms, seeds,
and accepted-trial counts remain dedicated result fields.

## Formulation options

[`formulation_options`](@ref) validates values that alter the mathematical
calculation represented by a formulation. Dispatch uses the formulation owner
type rather than the public construction selector:

```julia
formulation_options(LineParametersFormulation, FormulationOptions(reduce_bundle=false))
```

The default line-parameter formulation owns:

- bundle and Kron reduction.
- ideal transposition.

The normalized `FormulationOptions` record is stored in
`LineParametersFormulation.options`; read its payload through `.data`.
`PSCADFormulation` uses the shared physical options and currently requires
unreduced, untransposed matrices.

Connection assignments use one-based active phase IDs. A zero assignment marks
a grounded/eliminated conductor, while repeated active IDs identify conductors
that belong to the same bundle.

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
Selected equation controls belong to `FormulationOptions`, even when they set
quadrature or iteration tolerances. Backend execution options govern output,
tracing, logging and callbacks.

Public `options=(...)` keywords accept ordinary named tuples or the appropriate
owned record. They wrap tuples once before passing them to the owner. Direct
`computation_options(Owner, options)` calls require `ComputationOptions`; direct
`formulation_options` calls require `FormulationOptions`. The owner checks keys,
fills defaults and validates values. Constructing a wrapper alone does none of
that. Normalized controls are forwarded without reapplying normalization.

The three record types preserve the exact payload type, including callback
and sampler types. Access payload fields through `.data`. They provide no
implicit conversions or tuple forwarding. Immutability is shallow: contained
arrays are not copied or frozen. Nested numerical groups, physical parameters,
scientific products, axes and plotting attributes remain ordinary named tuples.

`Combinatorial`, `LinearError` and `MonteCarlo` normalize their own execution
options inside their constructors, including positional construction.
All Monte Carlo controls live in `MonteCarlo.options` and pass through
`computation_options(MonteCarlo, options)`. Keyword shorthand remains available:

```julia
MonteCarlo(formulation; options=(trials=100, seed=42, retain_details=true))
MonteCarlo(formulation; trials=100, seed=42, options=(retain_details=true,))
```

These two calls are equivalent. Supplying the same key both as a keyword and
inside `options` raises `ArgumentError`. Unknown keys and invalid values also
raise `ArgumentError`.

`ParametricProblem(space, ComputationOptions(...))` stores options for the **inner computation**.
Grid, batch, combinatorial and uncertainty traversal forward that record to the
selected core solver, whose `computation_options` method validates it. The
problem cannot normalize it at construction because the solver has not yet
been selected. Traversal retention and Monte Carlo sampling controls belong to
the higher-order formulation's own `ComputationOptions` record.

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
the retained matrices to `details(result).data.trace` after computation.

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
fixed-key normalized tuple. There is no general fallback and no conversion
from dictionaries, pairs, or `nothing`.

`MonteCarlo` owns a separate outer computation-option tuple. Its normalized
keys are `retain_details`, `on_error`, and `max_failures`. `on_error=:fail` is
the default and rethrows every exception. `on_error=:retry` requires
`retain_details=true` and rejects only `DomainError` realizations until the
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
type. The backend's `compute` method normalizes execution options before doing
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

struct ExternalFormulation{O <: FormulationOptions} <: AbstractFormulation
    options::O
end

function formulation_options(
    ::Type{ExternalFormulation},
    options::FormulationOptions,
)::FormulationOptions
    isempty(options.data) || throw(ArgumentError("unsupported formulation option"))
    return FormulationOptions()
end

function ExternalFormulation(;
    options::Union{NamedTuple,FormulationOptions} = FormulationOptions(),
)
    options = options isa NamedTuple ? FormulationOptions(options) : options
    return ExternalFormulation(formulation_options(ExternalFormulation, options))
end

function computation_options(
    ::Type{ExternalEngine},
    options::ComputationOptions,
)::ComputationOptions
    unknown = filter(key -> key != :tolerance, keys(options.data))
    isempty(unknown) || throw(ArgumentError("unsupported computation option"))
    normalized = merge((tolerance = 1.0e-8,), options.data)
    normalized.tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    return ComputationOptions(tolerance = Float64(normalized.tolerance))
end

function compute(
    ::ExternalEngine,
    problem,
    formulation::ExternalFormulation;
    options::Union{NamedTuple,ComputationOptions} = ComputationOptions(),
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    execution = computation_options(ExternalEngine, options)
    # Use `problem`, `formulation`, and `execution` here.
end
```

An external implementation does not need a dedicated options struct or private
wrapper around the two normalization functions. If it omits either Grammar
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
    return ComputationDetails(;
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

candidates = compute(problem, Formulation(earth_properties=Grid((:constant, :longmire1975))))
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


Scalar calculation selections are retained in `details(result).data.formulations`.
Its `requested` and `methods` fields hold complete requested and resolved records,
including physical parameters and numerical controls. Formula identifiers are available
as `record.requested.earth_admittance.identifier` (or through the corresponding
`air`, `earth`, `mixed` leaf). Resolved identities are under `methods`; a leaf's
`equivalent_earth` field contains its reduction choice and order. Realized local
shunt outcomes, including explicit fallback, remain in `details(result).data.shunt_model`.
There are no parallel `effective` or `numerical` inventories. Optional `trace`
contains detached evaluation data, not another formulation authority.
Reports use these records directly; they do not infer a selection from numerical
agreement or collapse different voltage references into the same formula label.
