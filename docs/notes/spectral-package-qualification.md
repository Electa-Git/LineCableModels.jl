# Spectral backend qualification — 2026-09-09

The subsequent fit-reuse and local-panel optimization is documented in
[`spectral-backend-optimization.md`](spectral-backend-optimization.md). The cost
table below records the original qualified implementation, before that follow-up.

The user's uncertainty requirement applies to the numerical dependency and to
its adapter, physical kernels and matrix solves. Keeping a Measurements element
type is insufficient: correlated derivatives must survive integration and the
complete K/H/L closure.

## Dependency decisions

- **QuadGK:** retain the existing dependency. Its public finite-breakpoint API
  and custom norm support complex Measurements values. The adapter uses nominal
  integration coordinates while keeping physical parameters uncertain.
- **DoubleExponentialFormulas 0.1.0:** added for transformed trapezoids. Its source
  is generic for complex Float32/Float64/BigFloat and Measurements callbacks.
  Acceptance requires normalization, mandatory feature subdomains and an
  independent error check; its raw stopping estimate is insufficient here.
- **ExpFit:** rejected; it is not a dependency. In the inspected `1.0.0-DEV`
  source, `hankel_matrix` constructs `Matrix{ComplexF64}`, while the amplitude
  solver requires `dt::Float64` and `ComplexF64` exponents. Executing the
  unmodified matrix-pencil and ESPRIT component sources demonstrated failures
  for Float32 spacing and Measurements inputs. The pencil zero-kernel case
  also failed. This was a component qualification, not installation or a claim
  to have tested every algorithm in the package.

The registered DE source tree is
`5a0a853bafeb3915e8a5b48a511d71838c1308c4` and its source files matched the
qualified source byte for byte. Downloaded source archives had SHA-256:

| Candidate | SHA-256 |
|---|---|
| DoubleExponentialFormulas | `957a2833c3e3b11143dba88d0a2be9f02befdd56a124289e628eba5b0c3d4a06` |
| ExpFit | `346c91fc967e4c701c6395f3275b7175a4d73b651c318ff11efb2ec962f668e0` |

Primary interfaces: [QuadGK](https://juliamath.github.io/QuadGK.jl/stable/api/),
[DoubleExponentialFormulas](https://machakann.github.io/DoubleExponentialFormulas.jl/stable/),
[ExpFit source](https://github.com/DOC-Package/ExpFit.jl).

## Demonstrated adapter requirements

DE's estimate contains a squared difference between successive integral values.
It is amplitude sensitive. A Gaussian centred at 0.371 with width 0.003 and
amplitude 1e-8 returned approximately 4.5188e-11 instead of 5.317e-11, despite a
reported error of 5.3e-19 and supplied feature breakpoints. A raw 1 MHz distant
voltage-kernel integral similarly differed by 1.26%, with reported error 8.69e-15.
Normalization and independent verification are necessary, including when the
backend has successfully converged by its own stopping condition.

A finite DE subinterval's affine map can round an interior node to the endpoint
1. The adapter's semi-infinite map must evaluate the nearest representable
interior coordinate instead of manufacturing an infinite Bessel argument. It
also avoids evaluating Bessel functions after the complete weighted value and
its uncertainty have underflowed to zero.

SpecialFunctions promotes some complex Float32 Bessel calls when the order is
an integer. Explicit Float32 wrappers preserve the physical scalar type. Its
missing complex BigFloat Bessel methods require the existing Engine-owned
extension, with scaled angular and exponentially decaying integral forms.

## Uncertainty evidence

An analytic complex exponential with uncertain amplitude and decay was checked
against its closed-form integral, including correlated derivatives. An
independent full two-wire 1 MHz closure with rho = 0.1 ± 0.001 Ω·m gave DE/quad
nominal relative differences around 2e-15 and relative uncertainty-component
differences below 2.3e-14.

The manuscript implementation was then exercised with independent uncertainty
in resistivity, radius, depth and separation. For the buried two-wire 1 MHz
case, DE/quad correlated differences were below 1.3e-16 relative. Finite-difference
derivative discrepancies were below 1.7e-5 relative. Mixed air/earth integration
also completed after the endpoint fix. Its resolved mutual derivatives agree;
some self derivatives are too small relative to their base values for a plain
Float64 subtraction test. Those checks must include subtraction roundoff or a
higher-precision reference; they are not evidence of a dependency conversion.

The reusable analytic regression is
`test/unit/engine/spectral_uncertainty.jl`, including a zero-nominal uncertain
integral and explicit CIM rejection. The complete public API passed 85 checks
of type inference, correlated material/geometry derivatives, DE/quad agreement,
shared Z/P preparation, central finite differences, and rejection of cache reuse
for independent uncertainties with identical nominal values and standard deviations.
CIM keeps its existing Float32/Float64
restriction; adding DE does not imply an uncertain matrix-pencil/SVD implementation.

The DE adapter supplies a numerical norm through a small scalar wrapper. This
avoids `abs(Complex{Measurement})` at zero nominal magnitude introducing NaN
derivatives into the package's error estimator. Sampling coordinates remain
nominal, while every physical value retains its derivative information.

## Production cost

Measured with Julia 1.12.7 and one BLAS thread, after compilation, at 1 MHz with
the accepted geometry. The reproducible driver is
`test/gauntlet/unified_earth_validation.jl`; its JSON records source/fixture
hashes, complete matrices, entrywise tolerance ratios and evaluation counts.
These are complete exterior Ze/Pe/Ye preparations in reused workspaces.

| Wires | Method | Minimum warmed time | Allocated bytes | Kernel evaluations |
|---:|---|---:|---:|---:|
| 2 | quad | 0.622 ms | 784 | 4,320 |
| 2 | trapz | 2.27 ms | 71,632 | 15,664 |
| 2 | cim | 772 ms | 528,285,776 | 94,476 |
| 3 | quad | 1.89 ms | 784 | 11,850 |
| 3 | trapz | 6.54 ms | 190,560 | 42,148 |
| 3 | cim | 1,778 ms | 1,198,708,552 | 214,788 |

CIM's verified image construction is much more expensive than quadrature; its
local matrix pencils and global least-squares decompositions dominate allocation.
Quadrature remains the default. The physical matrix and sensitivity arrays are
allocated once; the 784 bytes in the direct quadrature benchmark are primarily
small factorization metadata. There is no allocation in a Float64 spectral
callback. Trapz retains local DE-panel results while refining failed subdomains.

A separate three-cable, two-frequency public solver probe measured the old HEAD
at 3,616 bytes and 2.56 ms, the retained Xue selection at 3,872 bytes and 4.43 ms,
and the complete closure with Pollaczek series selection at 12,960 bytes and
7.81 ms. The original 4,096-byte retained-path regression ceiling is preserved.
The additional retained-author time buys explicit feature coverage and numerical
diagnostics. Matching default Z/P selections share one complete preparation.

## Resolved accuracy defects

For mixed air/earth voltage at 50 Hz, a cosine-residual CIM fit left Pe12 in
error by about 0.034 m/F and correctly failed the 0.00795 m/F final-entry budget.
The exact change `u²=λ²+κ_air²` removes the air-root singularity from the fitted
residual. The resulting image value differs from quadrature by roughly
`−1.42e−5+j3.29e−5 m/F` and passes the unchanged propagated budget. At exactly
κ_air=0 the already continuous cosine representation is retained.

Image candidates are drawn from every declared window before a weighted,
column-pivoted QR selection applies the image-count cap. Earlier first-come
truncation discarded later physical scales. The complete weighted residual and
true tail are still independently integrated; there is no quadrature fallback.
Candidates must also remain representable on the curved physical contour;
nonfinite columns are discarded before factorization and very large finite
columns are normalized during the fit. This prevents a LAPACK failure for a
mixed-layer case with prescribed Γ and unequal permeability.
If the initial pencil fit is unresolved, additional real exponential candidates
span the declared scales before the same weighted basis selection. This reduced
the difficult prescribed-Γ scalar integral's estimated error from approximately
1.3e-7 to 2.8e-9 without changing its kernel or acceptance tolerance.

The interface-reference voltage combines the receiver and endpoint expressions
before integration. The numerator contains
`I0(κ_receiver*r)*exp(-h_receiver*a_receiver)-J0(r*λ)`, evaluated with small-argument
series and `expm1`. This preserves the manuscript expression while avoiding
separate numerical errors in terms that nearly cancel.
Thirty direct checks validate this combination against the original endpoint
terms and the exact removable branch-point limits.

Radial validation and physical fitting retain the original λ coordinate. Forming
the small air root from `ag²-κg²+κair²` loses information when the two medium
roots differ substantially; this previously exhausted the quadrature verifier
in the scalar-potential diagnostic. The kernel can now evaluate `λ²+κair²`
directly, including through evaluation-counting wrappers. Thirty-three checks
compare the radial and original λ expressions with a ten-order root contrast.
For the scalar kernel, the smaller air root defines the radial coordinate even
for buried interactions. This removes its strongest small-scale feature from
the fitted residual; choosing the earth root compressed that feature beside a
much larger shift. The transformation is exact and the final value still comes
from the independently verified image sum.

The complete mixed-layer sweep in `test/gauntlet/unified_earth_references.jl`
passed 432 entrywise comparisons of Ze/Pe/Ye against quadrature. It uses unequal
radii, heights and permeability, prescribed Γ, and all four references: deep,
interface, finite depth and scalar potential. Both trapz and CIM satisfy the
unchanged component tolerances. The scalar-kernel regression also covers Γ=0.

Large electrical-size quadrature/trapz kernels combine their decay with the
analytic weight before exponentiation, preventing artificial overflow from a
split growing residual and vanishing weight.


## Repository quality check limitation

The broad engine/line-parameter suite passed 6,712 checks. The expanded public
uncertainty and scalar/air-root tests passed 267 checks; the spectral uncertainty
regression, including nominally equal uncertain media, passed 13 checks.
The subsequent engine/line-parameter regression passed all 4,565 checks after
the additional uncertainty, scalar, large-argument and allocation assertions.
After the stable radial evaluation change, all 3,888 frozen-matrix checks and
507 ownership/import checks passed again. The spectral, uncertainty and
retained-formula run passed 351 checks, including the new root-contrast limits
and scalar-coordinate comparisons.
The ownership/catalogue checks passed, and Aqua passed its ambiguity, type,
export, compatibility, stale-dependency and piracy checks.

Aqua's separate cold-precompile/persistent-task test fails in this local
Julia/depot environment. The original run fails to resolve an EarCut symbol
while precompiling GeometryBasics. Normal runtime triangulation succeeds.
Preloading the actual EarCut library removes that loader error, but the check
still times out, including with a 120-second budget. The untouched HEAD was
checked in an isolated archive and reproduces the same persistent-task failure.
The repository test is retained unchanged; this check is not reported as passing.
