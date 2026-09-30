# Physical PML mesh qualification

Execution plan: [fem-pml-physical-mesh.md](../../../docs/plans/fem-pml-physical-mesh.md).
The optional production control is `pml_resolution=(interpolation_cells=72,
coefficient_change=0.12)`. Qualification preceded implementation. Evidence
is in `.linecablemodels/fem/pml-physical-mesh/`; `live.log` is the single log.
The completed selection retains a 24-skin-depth physical domain. It does not
qualify arbitrary finite longitudinal propagation constants. The authorized
two-skin-depth probes were subsequently completed and rejected below.

## Why this distribution

Let `u` be outward distance divided by layer thickness `L`. The current cubic
stretch is `s(u)=1+(1-i)b*u^3`, with integrated coordinate
`z(u)=u+(1-i)b*u^4/4`. Material propagation, tangential wavenumber and the
longitudinal model determine `Q=L*sqrt(k_t^2+gamma^2-Gamma^2)`.
This qualification used the Gamma-to-zero quasi-fw equations. Finite Gamma
was a theoretical stress test here. Separately implemented finite-Gamma support
now supplies transverse wavenumbers and the solver's stretch slope to the mesh
owner; the coupled results below do not qualify that additional physics.

For the scalar modal equation

```math
-\frac{d}{du}\left(\frac{1}{s}\frac{dv}{du}\right)+Q^2s v=0,
\qquad v(0)=1,\quad v(1)=0,
```

the exact input response is `T=Q*coth(Q*z(1))`, or `1/z(1)` at `Q=0`.
The outgoing response is `Q`. The harness separates finite-wall reflection
from the error of the discrete boundary response. A static/grazing mode cannot
have a uniform relative outgoing-reflection guarantee from a finite layer.

The leading P1 derivative-interpolation contribution is proportional to
`integral(abs(v'')^2*h^2/abs(s))/12`. Minimizing this estimate at fixed cell
count gives a node density proportional to its integrand coefficient raised
to the one-third power. The candidate takes an envelope over prescribed modes,
normalized by each finite-layer input response and weighted by
`exp(-2*real(q)*clearance)`. This is an interpolation estimate and mode-selection
assumption, not an error bound on coupled cable conductance.

The modal spectrum includes normal propagation, air grazing cutoffs, and
evanescent tangential wavenumbers from `0.01/clearance` to `12/clearance`.
It includes both media on the side, air at the top and earth at the bottom.
The static mode also measures transformed-coefficient error. The sampled
density resolves the transition `b*u^3≈1`. This design uses wavelength,
attenuation, clearance and tensor variation; no grading parameter is fitted
to analytical cable conductance.

Six geometric strips approximate this density, with twelve cells per strip.
They use native Gmsh transfinite curves/surfaces and share corner and reference
path nodes. There is no node relocation after meshing. A qualification-only
in-memory override changes the exterior CAD construction; it does not edit
production sources or conductor meshing.

## Initial modal results

At the retained physical extent and thickness, the worst weighted relative
side/top boundary-response errors are:

| Frequency, soil | Existing 144 | Proposed 72, six strips |
|---|---:|---:|
| 0.1 Hz, 0.1 ohm m | 0.1939% | 0.1685% |
| 21.54435 Hz, 0.1 ohm m | 0.1972% | 0.1218% |
| 1 MHz, 0.1 ohm m | 0.1955% | 0.03193% |
| 1 MHz, 100 ohm m | 0.1764% | 0.02954% |

These are scalar modal metrics, not cable errors. Unweighted errors, including
modes attenuated before reaching the PML, remain in `modal.csv`/`strips.csv`.
The finite-wall reflection near grazing is about 0.357 in both distributions;
redistribution cannot fix continuous-layer truncation. A 16-to-32-point
quadrature check changes the low-frequency fitted-grid modal response by at
most 5.82e-8 at the original domain. That is small compared with the reported
modal discretization error, not a bound on native GetDP quadrature.

The candidate has 41,472 corner triangles versus 138,240 in the retained
144/144/96 prescription: 70% fewer. Coupled solve cost and all raw signs still
require the native tests. No full sweep or production adoption is justified
by this scalar result alone.

Finite-Gamma stress probes at 0.99 times the earth propagation constant show
severe under-resolution at low frequency after domain reduction. They cannot
be advertised as validated by this design. They were not coupled-solver tests.

## Reproduction

`modal.jl` and `strips.jl` use Julia standard libraries only. `qualify.jl`
uses public `compute` with the implemented controls and skips completed cases.
`verify_prescription.jl` compares all 79 frozen layouts without solving.
`verify_execution.jl` checks managed compute, detached `.pro` and explicit resume.
`compare.jl` reads saved exploratory matrices and reports R/X/G/B, unclipped.

`prototype.jl` and `coupled.jl` record pre-feature exploration. Their temporary
source substitution requires the original geometry owner; do not use them to
run the implemented feature. All completed numerical data are retained.

```bash
tail -n 60 -F .linecablemodels/fem/pml-physical-mesh/live.log
```

## Coupled diagnosis

The initial interpolation-only prescription failed both probes: its G12 was
`+5.95269e-23` at 0.1 Hz and `+1.98139e-18` at 21.54435 Hz, versus analytical
values `-6.79782e-25` and `-7.03197e-20` S/m. All four G entries had wrong signs.
Native solves fell to 15.47/15.39 s and 221123/221459 DOFs, but those savings
are unusable with the observed error. No domain reduction followed this failure.

The same candidate mesh with 13 rather than 12 native triangle integration
points gave G12 `+5.952689e-23` S/m. This rules out that quadrature change as a
useful fix for this failure; it is not a universal quadrature-convergence proof.

Changing only side grading gave `+3.04713e-23`; changing only top grading gave
`+3.02630e-23` S/m at 0.1 Hz. Both used the original bottom mesh, and each
retained the original mesh in the other air-facing direction. This isolates
comparable side and top contributions rather than a bottom-only failure.

The largest cell-wise change `abs(log(s_right)-log(s_left))` increased from
0.1554 in the working side/top mesh to 1.0447 in the initial candidate.
Optimizing the scalar response norm had neglected this local coefficient
resolution in the coupled problem. The bottom change, in contrast, fell from
0.1740 to 0.0238. These are diagnostics, not a proof that a particular
coefficient tolerance guarantees a cable conductance error.

Adding the density `abs(s'/s)/0.15` to the existing interpolation envelope,
and splitting a native strip at `b*u^3=1`, prescribed 109/109/76 intervals at
0.1 Hz. Its realized maximum side change was 0.2012 after geometric fitting.
It reduced the G12 absolute error by 62.68 times, but G12 was still positive:
`+2.80794e-25` S/m. Self G signs recovered. Native time was 24.99 s with
334245 DOFs. This candidate also failed the mandatory mutual-sign gate.

The final two exploratory controls tighten that density prescription from
0.15 to 0.12, holding all physical extents fixed. This is a prescribed
qualification control, not a runtime error estimator or automatic retry.
Their result determines the bounded feasibility outcome; the total limit is
eight native frequency solves. There is no approval or production change
between trials.

`mesh-audit.csv` verifies identical conductor nodes/contours, PML entrance
nodes, voltage paths and reference nodes between the initial candidate and
the retained low-frequency native mesh. Physical air/earth meshes have the
same node/element counts but do not have bitwise-identical coordinates.
Accordingly this is a comparison with identical meshing targets, not a claim
of a completely frozen physical mesh. Source hashes remain unchanged.

## Surviving prescription and final selection

The 0.12 density control passed all four raw G signs at both probes:

| Frequency | G12, S/m | Analytical G12, S/m | Native seconds | DOFs |
|---|---:|---:|---:|---:|
| 0.1 Hz | -1.85293e-25 | -6.79782e-25 | 29.506 | 392215 |
| 21.54435 Hz | -2.89295e-20 | -7.03197e-20 | 24.342 | 335403 |

These timings are 19.2% and 32.4% below the retained 144/144/96 native
controls. They exclude Julia package loading, compilation and mesh preparation;
they are not a whole-study speed claim. G12 magnitude errors are 72.74% and
58.86%: the first improves on the retained control, the second is larger.
Signs and magnitude errors remain separate observations. All R/X/G/B values,
signed differences and changes from the retained results are in `components.csv`.

This one survivor proceeds to the existing 79-frequency final fixture selection
on the prototype, before production implementation. It keeps 24 skin depths,
the original PML thickness and continuous stretch. The eight exploratory solves
were consumed by diagnosis, so the smaller-domain steps have not been coupled
tested and no reduced-domain claim is made. The 30% time reduction remains a
study target, not a passed whole-study result.

The qualification override writes actual strips/counts under `resolved/<case>`.
Its legacy `pml_layers`/`pml_grading` fields are placeholders in the input record,
not descriptions of that mesh. This temporary mismatch stays in the manual
harness. A separate runtime/cache under the evidence root prevents prototype
meshes from entering production caches. The eventual feature must own resolved
prescriptions and cache identity explicitly for both managed and native export.

`qualify.jl` retains the original physical fixtures and raw comparison readers;
it now calls the implemented API directly. `detached.jl` and `postprocess.jl`
read earlier native parity and public plotting assessments. Keep the completed
`pml-conductance-cost` and `conductor-mesh-qualification` evidence directories:
they contain the independent references and physical inputs.

## Completed selection and limitations

All 79 frequencies / 167 source columns completed on the prototype before
production changes. All 363 nonzero reference conductance signs were retained.
The selection includes seven two-wire scans, buried and mixed three-wire
layouts, and screen, tube and sector fixtures. Sixteen R/X/G/B figures exist
as SVG and PNG under `plots/`; no data clipping is used.

Seven two-wire scans took 561.6 s versus 889.2 s for the retained fixed mesh,
a 36.8% elapsed reduction in these two batches. Each used four frequency
workers and one native thread. First-call compilation was 28.7 s versus
26.4 s respectively; no repeated-run confidence interval is implied. The
two serial low/mid-frequency probes independently measured 19.2%/32.4% native
solve reductions.

Summed native worker time across those seven scans decreased from 2822.2 s
to 1659.8 s (41.2%). This is native work summed over concurrent workers,
distinct from the 36.8% elapsed reduction. Both exceed the study's 30% target
at the seven-scan level; the two individual exploratory probes varied.

Preserved signs do not imply uniform magnitude convergence. The largest
two-wire G relative error remains 76.9%. Maximum two-wire self-R error is
0.466%, essentially the retained result. The mixed fixture's largest G error
increased from 56.0% to 65.3% in a tiny aerial self term, while buried self
terms remain within about 1.3–1.8% of the reference.

Screen G[3,1] is -4.9446341e-14 S/m versus the retained reference
-3.9614064e-14 S/m: 24.8% relative or 9.83e-15 S/m absolute error. The former
fixed mesh was much closer for that entry. This is a measured degradation,
not established roundoff. Tube and sector maximum G errors are 1.61% and
0.54% respectively. Thus this optional prescription addresses tested signs
and cost; it is not an accuracy upgrade for every fixture. Explicit fixed
controls remain supported; global defaults are unchanged.

The eight exploratory solves were consumed diagnosing the mesh-distribution
failure. Moving the PML entrance to two skin depths remains unqualified.
The user subsequently authorized four B/C probes beyond the initial eight-solve
budget. Their negative result is documented below. No smaller-domain preset
is adopted.

## Prepared smaller-domain meshes

`domain_meshes.jl` exports and meshes the four pending B/C cases without
calling GetDP. All use the original ten-frequency problem, the qualified
physical resolution controls and unchanged conductor targets. Each detached
bundle contains its own domain/thickness/stretch data. Native output and
per-case mesh hashes are retained under `domain-probes/`; the single live log
contains the preparation transcript.

| Physical halfwidth / PML thickness | Nodes at 0.1 Hz | Nodes at 21.54 Hz |
|---|---:|---:|
| Retained 24δ / 24δ | 99,199 | 86,799 |
| Prepared 2δ / 24δ | 67,022 | 56,924 |
| Prepared 2δ / 2δ | 86,686 | 69,554 |

Moving the PML entrance inward while retaining its thickness reduces actual
nodes by 32.4% and 34.4%. Also reducing thickness saves only 12.6% and 19.9%:
the stronger coordinate stretch requires more normal intervals. The low
frequency side/top counts grow from 126/126 to 150/150 for the thin layer.
The full partition and timings are in `domain-probes/mesh-costs.csv`.

These preparation measurements are mesh costs, separate from the solver
timing and accuracy results below. No coupled solves were consumed by preparation.

## Completed smaller-domain assessment: rejected

The user authorized four additional solves with "yes proceed". They ran
serially through each frozen native `.pro` entry, both source columns, with
Gamma=0 and all original material, excitation and conductor controls. No
production source changed for these probes. Native completion markers and
mesh hashes allow resumption without rerunning completed solves.

| Domain / PML | G12 at 0.1 Hz [S/m] | G12 at 21.54 Hz [S/m] | Native seconds, low / mid |
|---|---:|---:|---:|
| Analytical | -6.79782e-25 | -7.03197e-20 | — |
| Retained 24δ / 24δ | -1.85293e-25 | -2.89295e-20 | 29.51 / 24.34 |
| 2δ / 24δ | +3.48056e-22 | +1.85165e-17 | 20.57 / 16.44 |
| 2δ / 2δ | +3.40124e-23 | +3.24938e-18 | 27.01 / 22.12 |

Both reduced-domain choices give four wrong G signs at each frequency:
16 sign failures in 16 comparisons across the four solves. The thinner PML
reduces the error but does not restore the signs. Maximum self-R relative
errors remain small: 0.0421% / 0.00106% for the thick PML, and
0.0419% / 0.00487% for the thin PML. Thus this is not a recurrence of the
conductor self-resistance resolution failure.

DOFs are 263987 / 216389 for the thick PML and 342199 / 266591 for the thin
PML. Peak native RSS and separate constraint, assembly, solve and output times
are in each `domain-probes/<case>/solve.toml`. Each `components.csv` retains
every signed and absolute R/X/G/B error against the analytical reference and
the change from the retained 24δ result. These are native timings, excluding
Julia compilation and mesh generation; each control was measured once.

Decision: retain the 24δ physical domain and PML thickness with the new mesh
prescription. Wavelength/attenuation-aware distribution reduces its cost;
it does not make this tested 2δ reduction accurate for the very small real
admittance. Keeping the normal-wave absorption setting alone is insufficient.
The result does not prove smaller domains impossible under every other
discretization; it rejects these two fixed controls without another tuning
campaign, scientific acceptance logic in production, or sign correction.

Reproduce or resume only these controls with `domain_solves.jl b` then
`domain_solves.jl c`. Completed native results are read directly on repetition.
