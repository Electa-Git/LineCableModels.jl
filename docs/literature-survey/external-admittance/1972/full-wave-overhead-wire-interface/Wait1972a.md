# Wait generalized overhead-wire shunt admittance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite thin circular conductor of radius ``a`` at height ``h`` parallel to a planar interface. |
| Calculated quantities | Generalized full-wave shunt admittance and its Carson qTEM reduction |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | Equation (25) retains the solved modal propagation constant and the interface spectral term; the qTEM reduction follows from equations (30)–(36). |
| Main source | James R. Wait (1972) |
| Citation key(s) | Primary: `:Wait1972a`; earlier result: `:Kikuchi1956`; validity analysis: `:Pogorzelski1977`; excitation witness: `:Kuester1978` |
| Evidence status | Original PDF equations (23)–(36) checked; the shunt-admittance definition is explicit |

**Description.** Wait's generalized line equation assigns a shunt-admittance term to the same full-wave mode used by the companion series-impedance record.

**Expression.** The generalized admittance is

```math
Y=i2\pi\epsilon_1\omega\,[A+2(N-iM)]^{-1},
\tag{25}
```

with

```math
A=K_0\!\left(i\sqrt{k_1^2-\beta^2}\,a\right)
-K_0\!\left(i\sqrt{k_1^2-\beta^2}\sqrt{4h^2+a^2}\right),
\tag{26}

N-iM=\int_0^\infty
\frac{\exp\!\left[-u_1\sqrt{4h^2+a^2}\right]}
{u_1+(k_1/k_2)^2u_2}\cos(\lambda a)\,d\lambda,
\tag{28}

u_1=\sqrt{\lambda^2+\beta^2-k_1^2},
\qquad
u_2=\sqrt{\lambda^2+\beta^2-k_2^2}.
```

It enters the implicit mode equation

```math
Z+\frac{\beta^2}{Y}=0,
\qquad i\beta=(ZY)^{1/2}.
\tag{23}
```

Under the qTEM conditions, equations (30)–(32) give

```math
i\beta_0\simeq(Z_eY_e)^{1/2},
\qquad
Y_e=i2\pi\epsilon_1\omega
\left[\ln\!\left(\frac{2h}{a}\right)\right]^{-1}.
\tag{30,32}
```

**Implementation.** For the generalized mode, evaluate (25), (26), and (28) at each trial ``\beta`` and solve (23) together with the companion impedance. For the qTEM limit, evaluate (32) directly.

**Limitations.** The conductor is infinite and thin, and the interface is planar. Equation (32) is the qTEM reduction and does not retain the full interface admittance correction present in (25) and (28).

**Reference.** [Wait1972a](@cite), equations (23)–(36), printed p. 678. The index also cites [Kikuchi1956](@cite), [Pogorzelski1977](@cite), and [Kuester1978](@cite).

## Evidence and approximation sources

- Equations (23), (25), (26), and (28) define the generalized shunt term.
- Equations (30)–(36) define the qTEM reduction and its relation to the modal root.

