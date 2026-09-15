# Computational engine

LineCableModels separates the physical problem, selected equations, and numerical
execution. Source equations and bibliography belong in the implementing formula
file. Backend/formula methods select an implementation through Julia dispatch.

Every family owns a `:default` routing identifier and may expose the concrete
equations it routes to explicitly. Literature references remain attached to
the equations without determining their software names.

| Family | Registered choices |
|---|---|
| Internal impedance | `:default`, `:schelkunoff1934` |
| Insulation impedance | `:default`, `:ametani1980` |
| Pipe impedance | `:default` |
| Insulation admittance, semicon admittance | `:default`, `:lossless`, `:lossy` |
| Local shunt geometry | `:default` (coaxial), `:coaxial`, `:boundary` |
| Earth impedance | `:default`, `:carson1926`, `:pollaczek1926`, `:gary1976`, `:wedepohl1973`, `:saad1996`, `:ametani2009`, `:lucca1994`, `:wise1934`, `:xue2018` |
| Earth admittance | `:default`, `:pollaczek1926`, `:wise1948`, `:xue2018` |
| Frequency-dependent soil properties | `:default` |
| Equivalent earth | `:default` |
| Modal transformation | `:default`, `:chrysochos2014` |

The internal default and `:schelkunoff1934` retain Schelkunoff's tubular
conductor expressions; insulation impedance's default and `:ametani1980`
retain the annular magnetic term documented by Ametani.
Both earth defaults implement the supplied circumferentially averaged framework
with complete enclosed-current normalization. The former air defaults remain
available as `:wise1934` for impedance and `:wise1948` for potential coefficients;
both former buried defaults are `:xue2018`. Dielectric `:default` selections
route to explicit `:lossless` equations; `:lossy` retains conductivity and the
material's supplied polarization losses. The FrequencyDependent default
preserves static properties.
EquivalentHomogeneous selects the basement when explicitly requested, and the
modal default performs Levenberg–Marquardt tracking. The EquivalentHomogeneous
default is a package-owned policy rather than an author equation.

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
unmodified `:lossless` insulation/semicon laws or their `:default` aliases.
Lossy or custom constitutive selections are unsupported by `:boundary` and raise `BoundarySolveError` for
eligible domains. Select `:coaxial` to retain their radial calculation without
losing conductivity or frequency dependence.

Inspect `details(result).shunt_model` for requested/effective model, domain
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
blueprints contain no boundary blocks and evaluate no boundary material law.
Boundary coefficients are independent of the frequency grid and external earth
model. They are reused across the frequency sweep, but a fresh `compute` call
constructs fresh blueprints: there is no process-global or cross-call cache.
Changed geometry and Monte Carlo realizations therefore receive new coefficients.
Do not mutate coefficient arrays shared by completed blueprints.

Production retains quadrature convergence, dense-storage budget, rank,
finite-value, reciprocity and positive terminal-capacitance checks. The
independent validation grid and derivative step-refinement are opt-in:

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
`nothing`, not zero. The audit does not alter the accepted terminal matrix.
The integration tolerances control dimensionless logarithmic moments, not
the error in an individual terminal coupling.

```@docs
BoundarySolveError
LineCableModels.Engine.ShuntModel.Formula
```

Strict failure is the default. An explicit `parameters=(fallback=:coaxial,)`
permits annular replacement only after recognized numerical or unsupported-law
failures, with a warning and recorded reason. Invalid inputs and unexpected
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
LineCableModels.Engine._shunt_load
LineCableModels.Engine._shunt_kernel_coefficients
LineCableModels.Engine._shunt_junction_exponent
LineCableModels.Engine._shunt_tape_faces
LineCableModels.Engine._shunt_log_moments
LineCableModels.Engine._shunt_capacitance
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
transfer surface callables. Unknown parameters or hooks fail at construction or indexed preflight. External
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
uses native method selection on this required signature and excludes the throwing
fallback. Domain-defining methods accept the three runtime payloads without extra
subtype constraints. Numerical specializations can optimize an admitted case.
Julia method availability determines which equations can be selected.

Each earth slot also accepts a NamedTuple for a physical air/soil two-half-space model:

```julia
selected = Formulation(earth_impedance = (
    air = formula(:carson1926),
    earth = formula(:pollaczek1926),
    mixed = formula(:lucca1994),
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
indexing has no defined origin in the present geometry definition.

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
signatures. Contribution overrides that declare integration resources receive the
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

The default selects the deepest soil layer. Martins-Britto et al.
[Martins-BrittoLopes2020](@cite) found that deep-layer conductivity predominated
in magnetic ground-return impedance for the multilayer soil cases they studied.
This provides a qualified rationale for the default resistivity. Accuracy depends
on layer contrasts and frequency; selecting one layer does not implement the
paper's equivalent-conductivity formula or establish the accuracy of the selected
permittivity and permeability. Choose another formula or supply a
`contribution` hook to change the material-selection rule.

`:after` applies the selected frequency law to physical layers first. `:before` reduces static
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

This replaces only an admitted case and receives the trailing runtime arguments shown above.
An algebraic replacement of an integral rejects unused integration controls;
an integral replacement declares its own integration section. Small physical
hooks retain the operation they customize. No callback inherits an unrelated
formula's numerical options or expands its physical domain.

`InternalImpedance.surface_impedances(resolved_formula, r_in, r_ex, rho, mu_r, jω)`
returns `(inner,outer,transfer)` coefficients in Ω/m, with hooks and per-kind
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

Each surface can select a different registered formula that implements that
surface, with its own parameters, hooks and numerical options. A replacement
for the transfer coefficient belongs to the transfer leaf, as
`transfer=formula(:default; hooks=(transfer=my_transfer,))`. Identical complete
selections share their prepared conductor state. Scalar shorthand retains its
existing numerical behavior. The internal term is called `transfer`; earth
`self`/`mutual` interaction names are unchanged.

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
requires more work. Candidate image orders reuse each window's projected Hankel
factorization before requesting more samples. Rectangular pencils use all samples
while limiting their row count by the candidate image order. Sample and Hankel
buffers are reused.
Amplitude fitting uses the integration contour and analytic tail penalties;
an independent integral of the absolute weighted residual checks its complete
continuation. With an analytic tail bound, this integral uses a finite interval
and adds bounds on both the true-kernel and image-expansion remainders. Quadrature
verifies the image sum and never supplies a value
labeled `:cim`. Geometry certificates integrate an absolute residual envelope
over the complete contour, including the tail. They apply to greater weight
heights and smaller separation-plus-radius on that contour. New geometries
outside the certified range require another certificate or fit; valid cache hits
evaluate the images without quadrature or refitting. One kernel prototype
evaluation still checks the physical scalar requirements. Reported certificates
remain numerical estimates.
Before caching a fit, continuation bounds are tightened to the measured finite
fit-error scale. Reuse selects the strongest applicable certificate, so a loose
tail estimate does not repeatedly force matrix-level refinements of the same fit.
The assembler can retain a certified expansion of an insignificant correction
even when that correction misses its local relative tolerance. Subsequent reuse
still requires its certificate to meet the newly requested absolute error budget.
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

The generalized half-space kernels provide typed absolute remainder envelopes
for overhead, underground and ordered mixed interactions, including voltage
reference paths and radial transformations. The shared planner finds a cutoff
from the requested `rtol`/`atol` budget and rechecks it against the computed value.
Only features inside this justified range seed finite panels and image fits.
There is no universal maximum wavenumber. Kernels without a usable analytic
envelope retain verification over the infinite interval; uncertain physical
inputs also retain that path because nominal envelopes do not bound derivatives.

For both trapz and CIM, `samples=nothing` (the default) selects adaptive
construction. An integer `samples=N`, with `N≥16`, caps **construction kernel
evaluations per scalar integral**. Trapz counts evaluations performed by the DE
rule; CIM counts coverage probes, pencil samples and amplitude-fitting samples.
Repeated evaluations count again. Prototype, pilot and independent verification
evaluations are outside this budget. An accepted cached image fit requires no
construction samples. The budget never grows silently: insufficient samples
raise an explicit error. This intentionally replaces the former starting-density
meaning of `samples`; an old explicit value such as `128` may need increasing.
Numerical tolerance remains a separate requirement:

```julia
integration = (method = :trapz, options = (rtol = 1e-6, atol = 0.0, samples = nothing))
integration = (method = :cim, options = (rtol = 1e-6, samples = 8192))
```

`SpectralEstimate.samples` reports construction evaluations and
`evaluations-samples` reports verification and pilot evaluations. Earth workspace
diagnostics aggregate both counters across integral solves and matrix refinements.

For trapz, `max_refinements` limits DE levels and `max_tail_refinements` limits
cutoff searches and independently verified subdomain refinements. Local phase
and envelope variation guide initial subdivision. Independent panel
quadrature references supply normalization and verification in one pass; their
absolute errors cannot cancel across panels. Successful panels are retained;
the DE package reuses nested samples within each panel. Rule tables are reused for identical
level limits and precision. Mathematical feature metadata is required to guard
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
external cases and each required reduction case. Absence of a selected reduction remains `nothing` in these calculation records.

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
        :lossy,
        :default,
    )),
    earth_impedance = Grid((
        :pollaczek1926,
        :saad1996,
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
attempt, target trial, failure stage, realized argument tuple, error type and
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

The normalized named tuple is stored in `LineParametersFormulation.options`.
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
Formula numerical options select how the owning equation is evaluated. Backend
execution options govern output, tracing, logging and callbacks.

Callers supply ordinary named tuples; `ComputationOptions` is an alias for
`NamedTuple`, not a separate type to construct or cast to. The owning
`computation_options(Owner, options)` method checks supported keys, fills
defaults and validates values. Its concrete returned tuple retains the types
of callbacks and other supplied values.

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

`ParametricProblem(space, options)` stores options for the **inner computation**.
Grid, batch, combinatorial and uncertainty traversal forward that tuple to the
selected core solver, whose `computation_options` method validates it. The
problem cannot normalize it at construction because the solver has not yet
been selected. Traversal retention and Monte Carlo sampling controls belong to
the higher-order formulation's own `options` tuple.

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

candidates = compute(problem, Formulation(earth_impedance=Grid((:default, :pollaczek1926))))
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
