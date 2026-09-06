# Wait generalized overhead-wire series impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite thin circular conductor of radius ``a`` at height ``h`` parallel to a planar interface. |
| Calculated quantities | Generalized full-wave series impedance and its Carson qTEM reduction |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | Equation (24) retains the solved modal propagation constant and the direct/image Bessel term plus the interface spectral integral. Equations (30)–(36) give the qTEM limit. |
| Main source | James R. Wait (1972) |
| Citation key(s) | Primary: `:Wait1972a`; earlier result: `:Kikuchi1956`; validity analysis: `:Pogorzelski1977`; modal witnesses: `:Kuester1978`, `:Olsen1978` |
| Evidence status | Original PDF equations (23)–(36) checked; the p.u.l. decomposition is explicit |

**Description.** Wait rewrites the modal equation as a product of a generalized p.u.l. series impedance and shunt admittance. This record contains the series term.

**Expression.** The mode equation and generalized series impedance are

```math
Z+\frac{\beta^2}{Y}=0,
\qquad
i\beta=(ZY)^{1/2},
\qquad\text{(23)}

Z=\frac{i\mu_0\omega}{2\pi}\,[A+2(Q-iP)],
\qquad\text{(24)}
```

where

```math
A=K_0\!\left(i\sqrt{k_1^2-\beta^2}\,a\right)
-K_0\!\left(i\sqrt{k_1^2-\beta^2}\sqrt{4h^2+a^2}\right),
\qquad\text{(26)}

Q-iP=\int_0^\infty
\frac{\exp\!\left[-u_1\sqrt{4h^2+a^2}\right]}
{u_1+u_2}\cos(\lambda a)\,d\lambda,
\qquad\text{(27)}

u_1=\sqrt{\lambda^2+\beta^2-k_1^2},
\qquad
u_2=\sqrt{\lambda^2+\beta^2-k_2^2}.
```

The qTEM reduction is

```math
Z_e=\frac{i\mu_0\omega}{2\pi}
\left[\ln\!\left(\frac{2h}{a}\right)-J_c\right],
\qquad\text{(31)}

J_c=\frac{2}{k_2^2}\int_0^\infty
(u-\lambda)e^{-2\lambda h}\,d\lambda,
\qquad
u=\sqrt{\lambda^2-k_2^2}.
\qquad\text{(33)}
```

**Implementation.** For the generalized mode, solve (23) while reevaluating (24), (26), and (27) at the trial ``\beta``. Select square-root branches with nonnegative real parts. For the qTEM limit, evaluate (31)–(33).

**Limitations.** The wire is thin and infinite. The generalized expression uses two homogeneous half-spaces and a planar interface. The qTEM form further assumes small electrical height and radius under the conditions stated before equation (29).

**Reference.** [Wait1972a](@cite), equations (23)–(36), printed p. 678. See [Kikuchi1956](@cite), [Pogorzelski1977](@cite), [Kuester1978](@cite), and [Olsen1978](@cite) for the index-listed related analyses.

## Evidence and approximation sources

- Equations (23)–(28) define the generalized impedance–admittance decomposition.
- Equations (29)–(34) define the Carson qTEM limit.
- Equations (35)–(36) restate the generalized result as an implicit mode equation.

