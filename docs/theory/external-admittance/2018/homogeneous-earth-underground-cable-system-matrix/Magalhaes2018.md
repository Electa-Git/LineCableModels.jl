# Magalhães–Barros–Lima homogeneous-earth underground-cable admittance matrix

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Insulated circular cables represented by their axes; the outer insulation radius enters the self term. Each core/sheath pair receives identical earth-return entries. |
| Calculated quantities | Single-cable and multiconductor per-unit-length earth-return admittance matrices; self and mutual entries of the intermediate matrix ``\Lambda-T``; source-provided series combination with insulation admittance |
| Earth structure | Homogeneous medium 1 beneath homogeneous medium 2, separated by a plane interface. |
| Model and approximation | Integral and matrix representation under the source's quasi-TEM reduction of the cited full-wave model. The paper does not reproduce the parent longitudinal-propagation prescription. |
| Main source | A. P. C. Magalhães et al. (2018). |
| Citation key(s) | `:Magalhaes2018` |
| Evidence status | Original publication and accepted-manuscript page images checked. |

**Description.** Per-unit-length earth-return admittance of one or more insulated cables buried in a homogeneous lossy half-space below another homogeneous medium. The source forms self and mutual entries of ``\Lambda`` and ``T``, assembles those entries into matrices, and inverts the complete matrix ``\Lambda-T``.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent ``\Gamma`` remains in equations (4)–(13). The source says these equations are a quasi-TEM approximation of a full-wave model with two propagation constants, but does not restate the numerical prescription by which the longitudinal constant was eliminated. | Stated and unresolved — discussion preceding (2), p. 663; conclusion, p. 668. |
| Air propagation constant ``γ_air`` | Medium 2 uses ``\gamma_2=\sqrt{j\omega\mu_0(\sigma_2+j\omega\varepsilon_2)}`` and ``u_2=\sqrt{\lambda^2+\gamma_2^2}``; the displayed equation permits a lossy/dispersive upper medium. | Stated — definitions below (6), p. 663. |
| Earth propagation constant ``γ_earth`` | Medium 1 uses ``\gamma_1=\sqrt{j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)}`` and ``u_1=\sqrt{\lambda^2+\gamma_1^2}``. | Stated — definitions below (6), p. 663. |
| Earth permittivity and displacement current | Retained through ``\sigma_1+j\omega\varepsilon_1`` in both ``\gamma_1`` and the prefactor of ``Y_{\mathrm{ext}}``. The authors emphasize that neglecting ``\omega\varepsilon_1`` in comparison with ``\sigma_1`` does not justify neglecting the admittance correction. | Stated — (4), definitions below (6), and following paragraph, p. 663. |
| Range of validity | Quasi-TEM upper limit stated as 10 MHz for one cable and at least a few MHz in the paper's tested cases. Mutual terms additionally require ``\lambda(f_{\max})\ge x_{\max}``; equation (18) and Table I give a soil-dependent estimate that the authors do not present as a final bound. | Stated — pp. 663–664 and 668–669, especially (18) and Table I. |
| Earth permeability ``μ_earth`` | Fixed to free-space permeability: ``\mu_0=\mu_1=\mu_2``. | Stated — below (4), p. 663. |
| Arrangement | Underground; one cable or a system of parallel cables at depths ``h_i`` with horizontal separations ``x_{ij}``; self and mutual terms. | Stated — Figs. 1–2 and (5)–(13), pp. 663–664. |
| Earth structure | Homogeneous medium 1 beneath homogeneous medium 2, separated by a plane interface. | Stated — Figs. 1–2 and definitions below (6), p. 663. |
| Conductor and insulation geometry | Insulated circular cable. Only the outer insulation radius ``r`` enters the earth-return term; ``r_0`` is the conductor radius in the separately printed coaxial insulation term. Each physical cable's core/sheath pair receives a ``2\times2`` block of identical earth-return entries. | Stated — Fig. 1, (1), paragraph below (6), and matrix text below (9), pp. 663–664. |
| Constitutive and field assumptions | Linear homogeneous media; quasi-TEM approximation; arbitrary ``\sigma_i,\varepsilon_i`` in the displayed two-media expressions; equal nonmagnetic permeabilities. Internal conductor effects and insulation fields are separate matrices. | Stated — discussion preceding (2), (4)–(6), and (9)–(13), pp. 663–664. |
| Conventions | ``j\omega`` phasor convention as printed; ``y`` is downward in Figs. 1–2, depths ``h_i>0``, conductors run along ``z``, and ``x_{ij}`` is signed in ``e^{-jx_{ij}\lambda}``. Cable voltage is measured with respect to ground, with nonzero ground-surface potential retained. Output is per unit length. | Equation-implied and stated — Figs. 1–2, paragraph preceding (2), and (4)–(13), pp. 663–664; appendix discussion, p. 669. |

**Expression.** Single-cable and multiconductor earth-return admittance matrices, source equations (4)–(8) and (11)–(13), printed pp. 663–664.

```math
Y_{\mathrm{ext}}=2\pi(\sigma_1+j\omega\varepsilon_1)[\Lambda-T]^{-1}.
\qquad\text{(4,12)}
```

For a single cable,

```math
\begin{aligned}
\Lambda&=K_0(r\gamma_1)-K_0(d\gamma_1) \\
d&=\sqrt{4h^2+r^2},
\end{aligned}\qquad\text{(5)}
```

```math
T=\int_{-\infty}^{\infty}
\frac{u_2}{u_1}
\frac{e^{-hu_1}-e^{-2hu_1}}{n^2u_1+u_2}
e^{-jr\lambda}\,d\lambda.
\qquad\text{(6)}
```

For cables ``i,j``,

```math
\begin{aligned}
\Lambda_{ij}&=K_0(d_{ij}\gamma_1)-K_0(D_{ij}\gamma_1) \\
d_{ij}&=\sqrt{(h_i-h_j)^2+x_{ij}^2} \\
D_{ij}&=\sqrt{(h_i+h_j)^2+x_{ij}^2},
\end{aligned}\qquad\text{(7)}
```

```math
T_{ij}=\int_{-\infty}^{\infty}
\frac{u_2}{u_1}
\frac{e^{-(h_i+h_j)u_1/2}-e^{-(h_i+h_j)u_1}}
{n^2u_1+u_2}
e^{-jx_{ij}\lambda}\,d\lambda.
\qquad\text{(11)}
```

The complete definitions are

```math
\begin{aligned}
u_1&=\sqrt{\lambda^2+\gamma_1^2} \\
u_2&=\sqrt{\lambda^2+\gamma_2^2} \\
n&=\frac{\gamma_2}{\gamma_1},
\end{aligned}
```

```math
\begin{aligned}
\gamma_1&=\sqrt{j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)} \\
\gamma_2&=\sqrt{j\omega\mu_0(\sigma_2+j\omega\varepsilon_2)} \\
\mu_0&=\mu_1=\mu_2.
\end{aligned}
```

Each scalar entry is expanded to the source's ``2\times2`` all-ones block for the core and sheath of a cable. After assembling every self and mutual block, the source applies the matrix inverse in (12); it does not define isolated admittance entries as ``1/(\Lambda_{ij}-T_{ij})``. Equivalently, only as a naming aid for the already printed relation, ``P_{\mathrm{ext}}=[\Lambda-T]/[2\pi(\sigma_1+j\omega\varepsilon_1)]`` gives ``Y_{\mathrm{ext}}=P_{\mathrm{ext}}^{-1}``.

The source then combines insulation and earth-return admittances as complete matrices:

```math
Y=\left(Y_d^{-1}+Y_{\mathrm{ext}}^{-1}\right)^{-1},
\qquad\text{(3,13)}
```

where the one-cable insulation formula printed in (1) is

```math
Y_d=2\pi j\omega\varepsilon_d\left(\ln\frac{r}{r_0}\right)^{-1}.
\qquad\text{(1)}
```

**Approximation.** The integrals and matrix assembly are not an additional analytical approximation within the source's quasi-TEM model. Quasi-TEM itself is a physical approximation of the cited full-wave parent. The paper does not reproduce the parent derivation or the explicit longitudinal-propagation prescription, so that reduction is not reconstructed here.

**Limitations.** The model has one earth half-space and one upper half-space, fixes both permeabilities to ``\mu_0``, omits conductor skin/proximity from this family, and does not model finite cable length or nonparallel routing. The source reports possible abnormal mutual-impedance behavior near 10 MHz for ``\rho_1\le100\ \Omega\,\mathrm m`` and ``\varepsilon_{r1}\ge40`` (p. 664). The mutual validity condition depends on maximum horizontal separation, and equation (18) is not presented as a final bound. The square-root branches are not stated adjacent to the definitions. The exact full-wave-to-quasi-TEM longitudinal reduction remains unresolved. Matrix inversion is required; entrywise reciprocation would change the published formulation.

**Reference.** [Magalhaes2018](@cite).  A. P. C. Magalhães, M. T. Correia de Barros, A. C. S. de Lima, P. E. D. Rocha, and R. A. Meyberg, “Earth Return Admittance Effect on Underground Cable System Modeling,” *IEEE Transactions on Power Delivery*, 33(2), 662–670 (2018), DOI `10.1109/TPWRD.2017.2741600`, equations (1)–(13), pp. 663–664, and equation (18), p. 669.

**Transcription source.** Original final-publication page images. Every sign, exponential, integration limit, Bessel argument, geometry definition, material factor, block-matrix statement, and inverse in (1)–(13) was checked on printed pp. 663–664. The accepted manuscript was inspected as an independent layout witness only.

## Source transcription

In source order, equation (4) prints the companion earth-return impedance

```math
Z_{\mathrm{ext}}=\frac{j\omega\mu_0}{2\pi}[\Lambda+S],
```

with ``S`` from (6) or ``S_{ij}`` from (8). The paper calls this essentially identical to Pollaczek and claims novelty for the admittance expression instead; it is therefore retained here as context and not counted as a distinct 2018 external-impedance contribution.

```math
\begin{aligned}
S&=\int_{-\infty}^{\infty}\frac{e^{-2hu_1}}{u_1+u_2}e^{-jr\lambda}\,d\lambda \\
S_{ij}&=\int_{-\infty}^{\infty}\frac{e^{-(h_i+h_j)u_1}}{u_1+u_2}e^{-jx_{ij}\lambda}\,d\lambda.
\end{aligned}\qquad\text{(6,8)}
```

The appendix prints a heuristic maximum-frequency expression,

```math
f_{\max}=\mu_0\frac{\pi c^2}{\rho_1}
\left[\varepsilon_r\left(\varepsilon_r+\mu_0\varepsilon_r\varepsilon_0(2\pi c)^2\right)\right]^{-1/2},
\qquad\text{(18)}
```

under an infinite-soil approximation and the requirements that the field wavelength exceed both penetration depth and ``x_{\max}``. This token sequence is preserved as printed; its dimensional appearance is not repaired.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Y_{\mathrm{ext}}`` | unchanged | Assembled earth-return admittance matrix | ``\mathrm S/\mathrm m`` |
| ``\Lambda,T`` | unchanged | Assembled Bessel-difference and Sommerfeld-integral matrices before inversion | source-normalized; ``[\Lambda-T]`` is dimensionless |
| ``\Lambda_{ij},T_{ij}`` | unchanged | Self/mutual scalar entries replicated into cable blocks | source indices ``i,j`` |
| ``K_0`` | unchanged | Modified Bessel function of the second kind, order zero | dimensionless argument |
| ``r_0,r`` | unchanged | Conductor radius and outer insulation radius | metres |
| ``h,h_i,h_j`` | unchanged | Positive downward burial depths | metres |
| ``x_{ij}`` | unchanged | Signed horizontal separation | metres |
| ``d,d_{ij},D_{ij}`` | unchanged | Self image distance and direct/image pair distances | metres |
| ``\lambda`` | unchanged | Horizontal spectral variable | ``\mathrm m^{-1}`` |
| ``u_1,u_2`` | unchanged | Vertical spectral roots in media 1 and 2 | ``\mathrm m^{-1}`` |
| ``\gamma_1,\gamma_2`` | unchanged | Bulk propagation constants in earth and upper medium | ``\mathrm m^{-1}`` |
| ``n`` | unchanged | Ratio ``\gamma_2/\gamma_1`` | dimensionless |
| ``\sigma_i,\varepsilon_i`` | unchanged | Conductivity and permittivity of medium ``i`` | ``\mathrm S/m`` and ``\mathrm F/m`` |
| ``Y_d`` | unchanged | Insulation-layer admittance matrix/one-cable formula | ``\mathrm S/\mathrm m`` |
| ``P_{\mathrm{ext}}`` | added only in explanatory prose | Name for the complete source-normalized inverse-admittance operator | not a source symbol; no entrywise interpretation |

No source symbol was renamed in the transcribed equations.

## Evidence and approximation sources

1. Cable geometry, insulation formula, voltage reference, and quasi-TEM/full-wave relationship: p. 663, Fig. 1 and (1)–(3).
2. Single-cable earth expressions and full material definitions: (4)–(6), p. 663.
3. Mutual geometry, kernels, all-ones cable blocks, and complete matrix operations: (7)–(13), p. 664.
4. The paper explicitly says the impedance repeats Pollaczek's result and the admittance is novel: pp. 663 and 668.
5. Validity evidence and unresolved bound: p. 664 and appendix, pp. 668–669, especially (18) and Table I.
6. Ground-referenced voltage and the nonzero surface-potential consequence: discussion below (17), p. 669.

## Limitations and discrepancies

- Equation (18) is preserved literally. Its printed prefactor and bracket raise a dimensional question, classified here as a suspected published defect; no corrected form is supplied.
- The source says the mutual conductance/capacitance sign under its ground-referenced voltage convention can differ from a zero-surface-potential construction. That is a convention-dependent result, not evidence of active behavior or permission to alter matrix signs.

## Numerical interpretation

The potential evaluator retains the ground-surface voltage reference and nonmagnetic half-spaces. The exponential difference is evaluated without subtracting nearly equal numbers. It is algebraically identical to Xue2018b's ground-surface potential selection. The companion series kernel is retained under Xue2018b; no second impedance implementation is registered. The complete potential matrix is assembled before inversion.
