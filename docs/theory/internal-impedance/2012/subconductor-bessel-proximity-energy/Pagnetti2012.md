# Pagnetti et al. subconductor/Bessel internal-impedance method

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Circular solid cores and annular circular screen; arbitrary eccentric core position within the shown screen geometry. |
| Calculated quantities | Per-unit-length internal resistance and inductance of interacting solid conductors and a core inside a hollow screen, including proximity effect |
| Earth structure | None; this is an internal conductor formulation. |
| Model and approximation | Source conductors are discretized into constant-current filaments and cylindrical harmonics are truncated. In the analytic target-screen phase the hollow screen's outer radius is taken as infinite; the exact zeroth-order annular skin term (30) is substituted to reduce that error. |
| Main source | A. Pagnetti, A. Xémard, F. Paladian, and C. A. Nucci (2012) |
| Citation key(s) | `:Pagnetti2012` |
| Evidence status | Original publication page images checked |

**Description.** Semi-analytical proximity method that subdivides a source conductor into axial filaments, solves cylindrical Bessel expansions for a target solid or hollow conductor, and evaluates internal resistance and inductance from the resulting volume current and magnetic field.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Zero; fields/currents are invariant along the parallel conductors. | Stated — formulation in §§II–IV. |
| Air propagation constant ``γ_air`` | Not applicable; the external dielectric is treated magnetostatically for the internal calculation. | Stated — (15)–(20). |
| Earth propagation constant ``γ_earth`` | Not applicable. | Scope stated in §I and §V. |
| Earth permittivity and displacement current | Not applicable to the conductor diffusion model. | Stated by ``\xi=\sqrt{j\omega\mu\sigma}``. |
| Range of validity | Hollow-screen source solution approximates the outer radius as infinite; source says this is good when skin depth is smaller than screen thickness. Numerical subdivision and harmonic truncation control accuracy. | Stated — text before (29), p. 2066. |
| Earth permeability ``μ_earth`` | Not applicable; conductor relative permeabilities are retained. | Stated — definitions after (22), (32). |
| Arrangement | Two nearby solid conductors, or one solid core eccentrically inside a hollow cylindrical screen. | Stated — Figs. 1, 3, 4. |
| Earth structure | None; this is an internal conductor formulation. | Scope stated. |
| Conductor and insulation geometry | Circular solid cores and annular circular screen; arbitrary eccentric core position within the shown screen geometry. | Stated — Figs. 3–4. |
| Constitutive and field assumptions | Homogeneous isotropic conductors, harmonic diffusion, longitudinal current, circular boundaries; external field represented by filament sources. | Stated — §§II–IV. |
| Conventions | ``\xi_c=\sqrt{j\omega\mu_c\sigma_c}``, ``\xi_s=\sqrt{j\omega\mu_s\sigma_s}``; internal ``R`` and ``L`` are energy/loss definitions. | Stated — (11), (29), (40)–(42). |

**Expression.** The hollow-screen current is represented as

```math
\vec J_s(r_s,\phi_s)=\sum_{n=0}^{\infty}
[o_n\cos(n\phi_s)+p_n\sin(n\phi_s)]K_n(\xi_s r_s)\,\hat z,
\qquad\text{(29)}
```

with its zeroth term replaced by the exact annular skin-effect term

```math
\vec J_{skin}(r,\phi)=\frac{\xi_sI}{2\pi c_1}
\frac{I_0(\xi_s r)K_1(\xi_s c_2)+I_1(\xi_s c_2)K_0(\xi_s r)}
{I_1(\xi_s c_2)K_1(\xi_s c_1)-I_1(\xi_s c_1)K_1(\xi_s c_2)}\,\hat z.
\qquad\text{(30)}
```

Once the coupled coefficients have been solved, the source defines

```math
\begin{aligned}
R_{int}&=\frac{1}{\sigma|I|^2}\int_{cond}|\vec J|^2\,dS \\
L_{int}&=\frac{1}{|I|^2}\int_{cond}\mu|\vec H|^2\,dS,
\end{aligned}\qquad\text{(40,41)}
```

```math
I=\int_{cond}\vec J\cdot\hat n\,dS.
\qquad\text{(42)}
```

**Approximation.** Source conductors are discretized into constant-current filaments and cylindrical harmonics are truncated. In the analytic target-screen phase the hollow screen's outer radius is taken as infinite; the exact zeroth-order annular skin term (30) is substituted to reduce that error.

**Limitations.** The method is two-dimensional and assumes circular homogeneous conductors. The screen approximation is frequency/thickness dependent. It provides internal ``R`` and ``L`` after solving source-specific coefficient systems (27)–(28), (31)–(39); it is not a standalone scalar closed form.

**Reference.** [Pagnetti2012](@cite).  A. Pagnetti et al., *IEEE Transactions on Power Delivery* 27(4), 2012, DOI `10.1109/TPWRD.2012.2212466`, equations (11)–(44), printed pp. 2064–2067.

**Transcription source.** Original IEEE page images. The ``K_n`` hollow expansion, annular Bessel numerator/denominator, energy integrals, normalization, integration domains and current definition were visually verified.

## Source transcription

For the core in the shield, the source gives

```math
\vec J_1(r_1,\phi_1)=
\frac{I_1\xi_c}{2\pi a_1}\frac{I_0(\xi_c r_{1i})}{I_1(\xi_c a_1)}
+\sum_{n=1}^{\infty}[g_{1n}\cos(n\phi_1)+h_{1n}\sin(n\phi_1)]I_n(\xi_c r_1)\,\hat z,
\qquad\text{(43)}
```

with coefficients obtained from the printed filament-continuity systems (37)–(38). The arguments ``r_{1i}``/``r_1`` are copied as printed rather than harmonized.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\vec J_c,\vec J_s`` | unchanged | core and screen current densities | ``\mathrm{A/m^2}`` |
| ``I_n,K_n`` | unchanged | modified Bessel functions | source convention |
| ``a_1,c_1,c_2`` | unchanged | core radius, screen inner/outer radii | m |
| ``o_n,p_n,g_{1n},h_{1n}`` | unchanged | coupled angular-harmonic coefficients | source normalized |
| ``R_{int},L_{int}`` | unchanged | internal p.u.l. resistance/inductance | ``\Omega/\mathrm m``, ``\mathrm H/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

The paper contrasts its filament/Bessel construction with symmetric-current approximations and explicitly proposes it for proximity-aware internal impedances. Equations (40)–(42) are the source's physical output definitions; the coefficient system is a necessary numerical dependency.

## Limitations and discrepancies

- Equation (43) visibly uses ``r_{1i}`` in the zeroth Bessel factor and ``r_1`` in the higher harmonics; this token difference is preserved.
