# Zhang buried-conductor ground-admittance integral and evaluator

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | One buried insulated conductor; cable self geometry uses burial depth and outer radius. |
| Calculated quantities | Ground-admittance integral and asymptotic-extraction numerical evaluator |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The evaluator splits the spectral integral at a source-defined threshold, integrates spline moments on the finite interval, and evaluates the extracted tail with exponential integrals. It retains earth displacement current. |
| Main source | Boyuan Zhang, Jun Zou, Xuelong Du, Jaebok Lee, and Mun-No Ju (2017) |
| Citation key(s) | Primary: `:Zhang2017`; parent integral: `:Papadopoulos2010b`; Vance reduction: `:Vance1978` |
| Evidence status | Original PDF equations (11)–(18) and Appendix A checked |

**Description.** The paper gives a scalar ground-admittance integral for one buried insulated conductor and an asymptotic-extraction evaluator for its slowly damped, oscillatory tail.

**Integral.** Equations (11)–(16) give

```math
Y_g'=\frac{j\omega}{2\pi(\sigma_g+j\omega\varepsilon_g)}
\int_0^\infty[F(\lambda)+G(\lambda)]\cos(\lambda r)\,d\lambda,
\qquad\text{(11)}

F(\lambda)=\frac{1-e^{-2u_1h}}{u_1}
+\frac{2e^{-2u_1h}}{u_1+u_0},
\qquad\text{(12)}

G(\lambda)=
\frac{2u_1(\gamma_1^2-\gamma_0^2)e^{-2u_1h}}
{(u_1+u_0)(u_1\gamma_0^2+u_0\gamma_1^2)},
\qquad\text{(13)}

u_0=\sqrt{\lambda^2+\gamma_0^2+k_{x,0}^2},
\qquad
u_1=\sqrt{\lambda^2+\gamma_1^2+k_{x,1}^2},
\qquad\text{(14)}

\gamma_0^2=-\omega^2\mu_0\varepsilon_0,
\qquad
\gamma_1^2=j\omega\mu_0(\sigma_g+j\omega\varepsilon_0\varepsilon_{rg}),
\qquad\text{(15)}

k_{x,0}^2=\omega^2\varepsilon_0\mu_0,
\qquad
k_{x,1}^2=\omega^2\varepsilon_0\varepsilon_{rg}\mu_0.
\qquad\text{(16)}
```

**Asymptotic extraction.** Appendix A splits

```math
Y_g'=Y_T+Y_\infty,
\qquad
T=10\max_{i=0,1}\left|\sqrt{\gamma_i^2+k_{x,i}^2}\right|,
\qquad\text{(A.1,A.4)}
```

with the original integrand integrated on ``[0,T]`` and the tail replaced by

```math
F_{asy}(\lambda)=\frac{1}{\lambda}
-\frac{\gamma_1^2+k_{x,1}^2}{2\lambda^3}
+\frac{(\gamma_1^2+k_{x,1}^2)e^{-2\lambda h}}{2\lambda^3},
\qquad\text{(A.6)}

G_{asy}(\lambda)=\frac{\gamma_1^2-\gamma_0^2}{\gamma_1^2+\gamma_0^2}
\left[
\frac{e^{-2\lambda h}}{\lambda}
+\frac{\gamma_1^2-k_{x,1}^2}{2}
\frac{e^{-2\lambda h}}{\lambda^3}
\right].
\qquad\text{(A.7)}
```

Equations (A.9)–(A.12) evaluate the tail through generalized exponential integrals ``\mathop{\mathrm{Ei}}(n,z)``. A direct numerical tail integral of (A.6)–(A.7) is mathematically equivalent when the special function is unavailable.

The paper also records the logarithmic impedance and Vance relation

```math
Z_g'=\frac{j\omega\mu_0}{2\pi}
\ln\!\left(\frac{1+\gamma_1r}{\gamma_1r}\right),
\qquad
Y_g'=\frac{\gamma_1^2}{Z_g'}.
\qquad\text{(17--18)}
```

**Implementation.** Select square-root branches with nonnegative real parts. Integrate (11) directly or split at ``T``. For the split evaluator, use the original kernel on ``[0,T]`` and the extracted kernel above ``T``; the source uses piecewise cubic-spline moments for the finite interval.

**Limitations.** The geometry contains one buried conductor and homogeneous earth. The Vance reduction in (18) requires burial depth much larger than earth skin depth; the integral in (11) does not impose that reduction.

**Reference.** [Zhang2017](@cite), equations (11)–(18), printed pp. 895–896, and Appendix A, printed p. 900. Parent integral [Papadopoulos2010b](@cite); Vance reduction [Vance1978](@cite).

## Evidence and approximation sources

- Equations (11)–(16) were checked against the original page image.
- Appendix equations (A.1)–(A.12) were checked against the original page image.
- Equations (17)–(18) distinguish the logarithmic/Vance approximation from the parent integral.

