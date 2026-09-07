# Høidalen solid-core proximity correction with low-frequency subtraction

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Round solid cores separated by ``d_{km}`` in the tested pipe-cable geometry. |
| Calculated quantities | Core-to-core proximity contribution to per-length series impedance; published differential/common-mode multipliers for the symmetrical three-core case |
| Earth structure | Not applicable; the correction is independent of earth layering or burial depth. |
| Model and approximation | The source subtracts the low-frequency inductive term from Kane's expression and adds the ``k=1`` Dwight term. The sum over ``n`` remains infinite, and no discarded-series error bound is supplied. |
| Main source | Hans Kr. Høidalen, “Analysis of Pipe-Type Cable Impedance Formulations at Low Frequencies,” 2013, DOI `10.1109/TPWRD.2013.2272343` |
| Citation key(s) | `:Hoidalen2013` |
| Evidence status | Equation (36) and its dependencies checked. |

**Description.** Series-impedance contribution from eddy-current proximity interaction of round solid cable cores. The formulation subtracts a low-frequency inductive contribution from an earlier core-to-core series and includes an additional finite-solid-conductor term. It is studied for three symmetrical screenless insulated cores inside a common conducting pipe, but the correction itself contains no pipe or earth parameter.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated as an imposed parameter; no ``Γ`` or axial exponential occurs in (25),(26),(36). The introduction says the study addresses longitudinal propagation effects and neglects transversal ones; that wording does not explicitly prescribe ``Γ=0``. | Stated/equation-implied — manuscript p. 1, final introduction paragraph; pp. 4,7 equations. |
| Air propagation constant ``γ_air`` | Not applicable to this conductor-proximity correction; no air bulk propagation term. | Equation-implied — (26),(36); the output is the conductor correction, not an earth-return kernel. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Equation-implied — (36); separate ground impedance is explicitly outside the paper's analysis, p. 3 after (18). |
| Earth permittivity and displacement current | Not applicable to this correction. | Equation-implied — (26),(36); no earth constitutive parameter enters. |
| Range of validity | Final (36) explicitly takes ``\mu_{rk}=1``. Author tests symmetrical three-core geometry; the correction does not resolve unequal-cross-section mutual asymmetry and over/under-compensates different mode inductances at high frequency. The 1 Hz–1 MHz FEM range is a test range, not a universal bound. No numerical small-``g`` error bound is stated. Section V.A uses 20 as the subsequent pipe-series upper limit after discussing (3),(9), without a separate convergence discussion for (36). | Stated — below (36), p. 7; discussion p. 8; FEM setup p. 2; section V.A, p. 4. |
| Earth permeability ``μ_earth`` | Not applicable; ``\mu_0`` is the reference permeability in the conductor formula, not a fitted earth permeability. | Equation-implied — (26),(36). |
| Arrangement | Internal core-to-core interaction, not an overhead/underground earth classification. Published mode factors concern symmetrical three-core geometry inside a common pipe. | Stated — Fig. 1, p. 1; (30),(33), p. 4; section VI, pp. 7–8. |
| Earth structure | Not applicable; the correction is independent of earth layering or burial depth. | Equation-implied — (36); separate ``Z_g`` discussion p. 3. |
| Conductor and insulation geometry | Round solid cores, radius ``r_{1k}``, center-to-center separation ``d_{km}``; screenless insulated cores in the tested pipe cable. The predecessor corresponds to solid/filamentary interaction; the extra term is intended to account for two solid conductors. | Stated — Fig. 1 and Table 1, p. 1; section VI before (36), p. 7. |
| Constitutive and field assumptions | Uniform scalar conductor conductivity ``\sigma_{ck}``; ``z_k=r_{1k}\sqrt{j\omega\mu_0\mu_{rk}\sigma_{ck}}`` retains conductive skin behavior and no conductor displacement-current term. Parent (25) permits ``\mu_{rk}``; final (36) restricts it to one. Proximity is approximate, not a complete coupled-field solution. | Equation-implied — (25),(26), p. 4; stated — below (36), p. 7 and conclusion p. 8. |
| Conventions | Positive ``j\omega`` factors; explicit time/axial exponent not stated. Parent text calls ``\Delta Z_{prox,km}`` the additional impedance in conductor ``m`` from eddy currents in conductor ``k``. Common mode has equal core currents and return in the pipe; differential excitation uses opposing currents. Impedances are per length, plotted in ``\Omega/\mathrm{km}``. | Stated/equation-implied — p. 4 preceding (25), sections IV.D–E; FEM setup p. 2; figures pp. 4,7–8. |

**Expression.** Høidalen's correction, manuscript equation (36), p. 7. The modulus square, infinite summation, distinct Bessel ratios, and logarithm of a square root are retained exactly in their published layout.

```math
\Delta Z_{prox,km}=
\frac{j\omega\mu_0}{2\pi}
\left(
\sum_{n=1}^{\infty}
\frac{g^n}{z_k}\cdot\frac{I_n(z_k)}{I_{n-1}(z_k)}
\cdot\left|1+\frac{n\cdot k_1}{1-k_1}\right|^2
-\ln\sqrt{\frac{1}{1-g}}
\right).
\qquad\text{(36)}
```

Definitions immediately below (36), with the required conductor parameter from (26):

```math
\begin{aligned}
g&=\left(\frac{r_{1k}}{d_{km}}\right)^2 \\
k_1&=g\cdot\frac{I_2(z_k)}{I_0(z_k)} \\
\mu_{rk}&=1,
\end{aligned}
```

```math
z_k=r_{1k}\cdot\sqrt{j\omega\cdot\mu_0\cdot\mu_{rk}\cdot\sigma_{ck}}.
\qquad\text{(26)}
```

``r_{1k}`` and ``d_{km}`` are lengths; ``\sigma_{ck}`` is conductor conductivity, ``\mu_{rk}`` relative conductor permeability, and ``\mu_0`` reference permeability. ``g,k_1,z_k`` are dimensionless. The manuscript does not explicitly define the Bessel family, but its cited Kane et al. (1995) source, printed p. 1647, defines ``I_n`` as the modified Bessel function of the first kind of order ``n``. The manuscript does not state the square-root branch. The formula provides a proximity increment, not the whole core impedance, pipe surface impedance, or cable terminal impedance.

For the source's symmetrical three-core arrangement, its displayed differential/common-mode proximity contributions are:

```math
\Delta Z_{prox,1}=\Delta Z_{prox,km},\qquad\text{(30)}
```

```math
\Delta Z_{prox,0}=4\cdot\Delta Z_{prox,km}.\qquad\text{(33)}
```

These are the author's mode multipliers, not a general matrix assembly for unequal cables. Equation (36) modifies the pair contribution; a new arbitrary-geometry diagonal formula is not printed.

**Approximation.** The author explicitly starts from Kane et al.'s formula reproduced as (25), subtracts the low-frequency inductive expression (28), then divides the impedance by two. An additional term, identified as ``k=1`` in Dwight's ``B_n,C_n`` summations [11], supplies the factor involving ``k_1`` in (36). The author calls the inductance modification first order and says further terms have higher powers of ``g``; no error bound or fully specified discarded series is printed. The displayed sum over ``n`` remains infinite. “Closed form” in the source's prose is not interpreted as exact finite-form evaluation.

**Limitations.** The author states that unequal-cross-section mutual asymmetry remains and that high-frequency differential/common-mode inductance is respectively over- and under-compensated for the tested case. The root branch and publisher-edition match remain unresolved; Bessel kind is now corroborated in the cited original; no attempt is made to symmetrize or improve the formula.

**Reference.** [Hoidalen2013](@cite).  Hans Kr. Høidalen, “Analysis of Pipe-Type Cable Impedance Formulations at Low Frequencies” (2013), DOI `10.1109/TPWRD.2013.2272343`; actual inspected manuscript (36), p. 7, (25)–(28),(30),(33), p. 4, and discussion pp. 7–8.

**Transcription source.** Original-author manuscript-format witness, not confirmed publisher-edition pagination. Formula (36) and dependencies were checked visually against PDF pages 1,4,7; subsequent restrictions were read on p. 8. The [author's NTNU publication list](https://www.ntnu.edu/employees/hans.hoidalen) independently lists this title in 2013 and links the same DOI, establishing publication identity but not equality of manuscript and publisher equation pages. The DOI retrieval attempt failed; no publisher page was seen. Earlier Kane/Dwight expressions are secondary witnesses only where reproduced by Høidalen.

## Source transcription

The formula section gives (36) in source notation; (26) is its unchanged dependency. The following complete parent expressions preserve the actual comparison and approximation source basis, not alternative formulas silently blended into (36).

Equation (25), manuscript p. 4, explicitly attributed to Kane, Ahmad and Auriol (1995), reference [3]:

```math
\Delta Z_{prox,km}=\frac{j\omega\mu_0}{2\pi}\cdot
\sum_{n=1}^{\infty}
\frac{2\mu_{rk}\cdot(r_{1k}/d_{km})^{2n}}
{n\cdot(\mu_{rk}-1)+z_k\cdot\dfrac{I_{n-1}(z_k)}{I_n(z_k)}}.
\qquad\text{(25)}
```

Equation (26) follows, as transcribed above. Equation (27), the **parent's** multi-core diagonal contribution, is:

```math
\Delta Z_{prox,kk}=\frac{j\omega\mu_0}{2\pi}\cdot
\sum_{\substack{m=1\\m\ne k}}^{P}\sum_{n=1}^{\infty}
\frac{2\mu_{rm}\cdot(r_{1m}/d_{km})^{2n}}
{n\cdot(\mu_{rm}-1)+z_m\cdot\dfrac{I_{n-1}(z_m)}{I_n(z_m)}}.
\qquad\text{(27)}
```

Here ``P`` is the number of core conductors; the source uses the corresponding material/radius labels ``m`` in the inner summand. This is not labelled a corrected diagonal companion to (36), because the source does not print such a replacement.

The next paragraph states that ``z_m I_{n-1}(z_m)/I_n(z_m)`` approaches ``2n`` at low frequency and then prints:

```math
\lim_{\omega\to0}\Delta Z_{prox,km}=
\frac{j\omega\mu_0}{2\pi}\cdot
\frac{2\mu_{rk}}{\mu_{rk}+1}\cdot
\ln\left(\frac{d_{km}^2}{d_{km}^2-r_{1k}^2}\right).
\qquad\text{(28)}
```

The ``\lim`` on the left and the remaining ``\omega`` on the right are both source notation; no replacement by a reviewer-derived asymptotic symbol is made. Section VI describes removing this low-frequency inductive coefficient before dividing by two. The additional Dwight-based term leads to (36) on p. 7; the original Dwight summations are not reproduced and were not inspected.

The mode factors (30) and (33) are repeated above in their printed form. The full cable-mode formulas (29) and (32) also contain separate pipe, insulation, and core-skin terms. Those terms are documented in the paper's separate low-frequency pipe records, not redefined as part of the proximity correction.

## Notation map

Notation is unchanged, including the collision between conductor label ``k``, the source's prose ``k=1`` Dwight-series truncation index, and the separately defined factor ``k_1``. No new LCM ID or normalized propagation constant is introduced.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``\Delta Z_{prox,km}`` | unchanged | Pair proximity impedance increment | Per length; source describes impedance in ``m`` due to eddy currents in ``k`` |
| ``\Delta Z_{prox,kk}`` | unchanged | Earlier formula's diagonal proximity increment | Per length; (27) is a secondary parent, not corrected (36) |
| ``\Delta Z_{prox,1},\Delta Z_{prox,0}`` | unchanged | Differential/common-mode proximity increments | Source's symmetrical three-core multipliers |
| ``r_{1k},r_{1m}`` | unchanged | Solid-core radii | Length; not outer insulation radius |
| ``d_{km}`` | unchanged | Center-to-center conductor spacing | Length; not spacing between conductor surfaces |
| ``z_k,z_m`` | unchanged | Conductor skin arguments from (26) | Dimensionless; branch not stated |
| ``\sigma_{ck},\sigma_{cm}`` | unchanged | Core conductivities | ``\mathrm{S/m}``; corresponding ``m`` version in (27) |
| ``\mu_0,\mu_{rk},\mu_{rm}`` | unchanged | Reference permeability and relative core permeabilities | ``\mathrm{H/m}``, dimensionless; final (36) fixes ``\mu_{rk}=1`` |
| ``g`` | unchanged | Squared core-radius/axis-spacing ratio | Dimensionless; no author-specified numerical truncation threshold |
| ``k_1`` | unchanged | Additional solid-conductor factor | Dimensionless, ``g I_2(z_k)/I_0(z_k)`` |
| ``I_n,I_{n-1},I_0,I_2`` | unchanged | Modified Bessel functions, first kind, indicated orders | Not defined in this manuscript; definition image-corroborated in cited Kane original p. 1647 |
| ``n`` | unchanged | Summation order | Positive integers, no upper finite cutoff in (36) |
| ``k,m,P`` | unchanged | Conductor labels and total number of cores | ``m\ne k`` in diagonal parent (27) |
| ``\omega,j`` | unchanged | Angular frequency and imaginary unit | ``\mathrm{rad/s}``, ``j^2=-1``; positive ``j\omega`` factors |

## Evidence and approximation sources

- The title, byline, Fig. 1, and Table 1 were image-checked. NTNU's publication listing corroborates the title, year, and DOI but does not establish equation-version equality.
- The approximation is substantive: the author proposes subtracting an inductive term, applying a factor of two, and adding a Dwight-based term. It is not a notation-only reproduction of Kane et al.
- The source attributes (25) to Kane, Ahmad and Auriol, “Multiwire shielded cable parameter computation,” IEEE Transactions on Magnetics 31(3), 1646–1649 (1995). It explicitly doubts that paper's further attribution to Schelkunoff because the cited 1934 work does not address the same core-to-core proximity case. No original-priority claim is made here from that secondary chain.
- Dwight, “Proximity Effect in Wires and Thin Tubes,” Transactions AIEE (1923), pp. 850–859, is reference [11]. Only Høidalen's description of retaining one additional term is used here. A same-title *Journal of the AIEE* version has different pagination and is not substituted for the cited Transactions edition.
- Author-stated physical shortcoming: predecessor mutual terms differ under interchange for unequal radii/permeabilities, and a low-frequency inductive coefficient survives. The correction is not claimed to solve the remaining unequal-cross-section asymmetry (p. 8).
- Numerical evidence: comparison with FEM over the source's 1 Hz–1 MHz test sweep for a symmetrical pipe cable, and mode-dependent discussion in pp. 7–8. The reported trends do not define a universal error bound. Pipe material tests do not introduce pipe permeability into (36).
- The paper's finite-pipe low-frequency, surface/connection, and mode terms have a [separate record](../finite-pipe-low-frequency-surface-and-loop-terms/Hoidalen2013.md); the [infinite-wall low-frequency analysis](../infinite-pipe-low-frequency-logarithmic-terms/Hoidalen2013.md) is kept distinct. The earlier Kane, Brown–Rocamora, and Ametani formulas retain their original attributions.

## Limitations and discrepancies

- **Conversion defects:** Markdown loses the equality and lower summation bound in (25), interleaves unrelated columns in (25)–(28), and damages (27)–(28) with repeated placeholder-like LaTeX. Its (36) keeps the main modulus-square/log-root structure but loses the equality/lower bound and mixes surrounding figure/prose fragments. All displayed mathematics here follows the PDF image instead.
- **Definition supplied by the parent:** Høidalen does not explicitly define ``I_n``; the cited Kane original p. 1647 identifies the modified first-kind family. This resolves the Bessel-function kind without asserting full algebraic equivalence of Høidalen's parent expressions to Kane's original.
- **Author-stated model limitation:** unequal-core mutual asymmetry remains after (36); high-frequency inductance is not fully reproduced. No symmetrization or high-frequency correction is inserted.
- **Published asymptotic notation:** (28) writes a limit with a frequency-dependent right-hand side. It is retained as the source's low-frequency expression, not repaired.
- **Dwight dependency:** the factor in (36) is checked against Høidalen's manuscript, but the approximation-order derivation remains attributed through Høidalen rather than verified from the cited Transactions edition.
- **No fabricated singularity diagnosis:** the isolated factors ``1/z_k`` and ``1/(1-k_1)`` are preserved. This record does not infer a physical singularity from either without evaluating the complete expression in its stated domain.

## Numerical interpretation and implementation

`PipeImpedance.core_proximity(Val(:Hoidalen2013), ...)` evaluates
(36) independently of the chosen wall approximation. Selecting
`proximity=:Hoidalen2013` in a pipe formula adds this increment
through physical core radii, temperature-corrected resistivities,
and centre separations. These radii are not the insulation radii.

The source's three-core mode factors determine the supported
symmetric increment:

```math
\Delta\mathbf Z_{\mathrm{core}}
=\Delta Z_{prox,km}
\begin{bmatrix}
2&1&1\\
1&2&1\\
1&1&2
\end{bmatrix}.
```

Its differential eigenvalue is ``\Delta Z_{prox,km}``, and its
common eigenvalue is ``4\Delta Z_{prox,km}``. Matrix assembly is
restricted to three equal nonmagnetic solid cores with equal pair
separations. The implementation does not infer a matrix for unequal
cross-sections from these mode factors.

The static logarithm is subtracted analytically, term by term, using
``I_n(z)/(zI_{n-1}(z))-1/(2n)``. This preserves the complete
formula while retaining its small low-frequency difference in finite
arithmetic. Finite sums use numerical convergence controls, not a
source-wide error bound. The author's high-frequency inductance
limitation remains applicable.
