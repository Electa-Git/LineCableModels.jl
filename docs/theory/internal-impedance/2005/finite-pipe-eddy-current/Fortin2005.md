# Fortin–Yang finite-pipe eddy-current impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Arbitrary core locations within pipe. Circular pipe, core axes; coatings do not carry pipe current. |
| Calculated quantities | Pipe-type cable self and mutual impedance from return and eddy currents in a finite wall |
| Earth structure | None in displayed contribution. |
| Model and approximation | Infinite harmonic sums require truncation; vector potential is two-dimensional and pipe material is unsaturated. |
| Main source | Simon Fortin, Y. Yang, J. Ma, and F. P. Dawalibi (2005) |
| Citation key(s) | `:Fortin2005` |
| Evidence status | Original-page image verified |

**Description.** Vector-potential integration of finite-pipe return and eddy currents, yielding direct terminal self/mutual voltage drops.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in two-dimensional pipe cross-section. | Formulation. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal pipe formula. |
| Earth propagation constant ``γ_earth`` | Earth term is separate. | Scope. |
| Earth permittivity and displacement current | Not part of pipe term. | Scope. |
| Range of validity | Multiple cables inside a circular finite-thickness pipe. | Title/model. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal formula. |
| Arrangement | Arbitrary core locations within pipe. | Figure 1. |
| Earth structure | None in displayed contribution. | Scope. |
| Conductor and insulation geometry | Circular pipe, core axes; coatings do not carry pipe current. | Geometry. |
| Constitutive and field assumptions | Linear nonsaturated magnetic or nonmagnetic pipe. | Text before §3. |
| Conventions | Unit source current may be used for impedance extraction. | §2. |

**Expression.** With the vector potential from (3)–(11),

```math
Z_{self}=j\omega[A(a,0)-A(b,0)]-E_z(a,0),
\qquad\text{(12)}
```

```math
Z_{mutual}=j\omega[A(a,0)-A(r_j,\alpha)]-E_z(a,0).
\qquad\text{(13)}
```

The finite-wall current density uses ``p=\sqrt{j\omega\mu\sigma}`` and the Bessel/harmonic expansion in (1)–(5).

**Approximation.** Infinite harmonic sums require truncation; vector potential is two-dimensional and pipe material is unsaturated.

**Limitations.** External earth return and cable-core internal impedance are appended separately.

**Reference.** [Fortin2005](@cite), equations (1)–(13).

**Transcription source.** Original-author conference PDF images; signs and evaluation points verified.

## Source transcription

The positive-order vector potentials are printed as:

```math
A_1=\frac{\rho_s}{j\omega}\sum_{n=1}^{\infty}
[F_nI_n(pr)+C_nK_n(pr)]\cos(n\theta),
```

```math
A_2=\frac{\mu_0}{4\pi}\ln(\rho^2+b^2-2\rho b\cos\theta)
+\sum_{n=1}^{\infty}L_{2n}\frac{\rho^n}{n}\cos(n\theta),
```

```math
A_3=\frac{\mu_0}{4\pi}\ln(\rho^2+b^2-2\rho b\cos\theta)
+\sum_{n=1}^{\infty}L_{3n}\frac{1}{n\rho^n}\cos(n\theta).
\qquad\text{(6)}
```

Here ``a,c`` are the inner and outer pipe radii, ``b`` is the source-axis offset, ``\rho,\theta`` are observation coordinates, and ``\rho_s`` is pipe resistivity. The source uses ``r`` in the wall functions and ``\rho`` in the other regions. The boundary conditions are:

```math
\begin{aligned}
\frac{\partial A_1(a,\varphi)}{\partial\varphi}
&=\frac{\partial A_2(a,\varphi)}{\partial\varphi}, \\
\left.\frac{1}{\mu_s}\frac{\partial A_1(r,\varphi)}{\partial r}\right|_{r=a}
&=\left.\frac{\partial A_2(r,\varphi)}{\partial r}\right|_{r=a}, \\
\frac{\partial A_1(c,\varphi)}{\partial\varphi}
&=\frac{\partial A_3(c,\varphi)}{\partial\varphi}, \\
\left.\frac{1}{\mu_s}\frac{\partial A_1(r,\varphi)}{\partial r}\right|_{r=c}
&=\left.\frac{\partial A_3(r,\varphi)}{\partial r}\right|_{r=c}.
\end{aligned}
\qquad\text{(7)}
```

``\mu_s`` is the relative pipe permeability in these boundary conditions. Each positive angular order therefore requires four scalar coefficients, not a discretization of the core cross-sections. Solving these four equations gives the finite-wall harmonic response also used by the [Da Silva Method 3 record](../../2006/finite-pipe-without-core-proximity-auxiliary-model/DaSilva2006.md#identification-and-source). The numerical comparison uses continuity of each Fourier amplitude and of its permeability-weighted radial derivative, with the same positive pipe-return resistance.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``A`` | unchanged | longitudinal magnetic vector potential | ``Wb/m`` |
| ``E_z`` | unchanged | longitudinal pipe electric field | ``V/m`` |
| ``p`` | unchanged | pipe diffusion root | ``m^{-1}`` |

## Limitations and discrepancies

The inspected conference publication gives the year as 2005.
