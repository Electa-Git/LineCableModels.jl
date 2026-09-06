# Gustavsen armored three-core cable 2D FEM impedance model

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Three-core cable with helically pitched phase conductors, sheaths, and a surrounding layer of steel armor wires. |
| Calculated quantities | Two-dimensional FEM impedance matrix, including skin, proximity, armor losses, and positive- and zero-sequence impedances |
| Earth structure | Earth-return impedance is appended separately and is not part of this formula. |
| Model and approximation | Equal armor-wire currents are enforced. Pitching is represented by an energy-matched complex permeability in a fictitious nonconductive region, with the source's effective positive- and zero-sequence pitches. |
| Main source | Bjørn Gustavsen (2023) |
| Citation key(s) | `:Gustavsen2023` |
| Evidence status | Original PDF formulation and validation cases checked |

**Description.** The method solves a two-dimensional magnetic-vector-potential problem and modifies the armor region to represent the three-dimensional pitch. It returns the conductor impedance matrix, induced currents, and losses.

**Field equation.** The finite-element problem is

```math
\mathop{\mathrm{div}}\!\left(\frac{1}{\mu}\mathop{\mathrm{grad}}A\right)
-j\omega\sigma A=-J_s.
\qquad\text{(1)}
```

The one-step FEM system includes nodal vector potentials and conductor voltage gradients. Its terminal relation is

```math
-\frac{d\mathbf v}{dz}=\mathbf Z\mathbf i.
\qquad\text{(3)}
```

The conductor loss integral is

```math
\begin{aligned}
p&=\int_S\rho\,|J_s+J_e|^2\,dS \\
J_e&=-j\omega\sigma A.
\end{aligned}\qquad\text{(4)}
```

For ``h=\exp(j2\pi/3)``, the positive-sequence excitation and impedance are

```math
\begin{aligned}
\mathbf i_+&=[1\;h^2\;h]^T \\
Z_+&=\frac{1}{3}[1\;h\;h^2]\,[v_1\;v_2\;v_3]^T.
\end{aligned}\qquad\text{(6--7)}
```

The source groups phase, sheath, and armor conductors into ``\mathbf Z_{3\times3}``, bonds sheath and armor through

```math
\mathbf Z_{2\times2}
=\left(\mathbf P\mathbf Z_{3\times3}^{-1}\mathbf P^T\right)^{-1},
\qquad
\mathbf P=\begin{bmatrix}1&0&0\\0&1&1\end{bmatrix},
\qquad\text{(11--12)}
```

and calculates

```math
Z_0=3\left(Z_{11}^{2\times2}
-\frac{Z_{12}^{2\times2}Z_{21}^{2\times2}}{Z_{22}^{2\times2}}\right).
\qquad\text{(13)}
```

The pitch angles and effective armor-wire angles are

```math
\begin{aligned}
\alpha&=\arctan\!\left(\frac{2\pi R_a}{P_c}\right) \\
\beta&=\arctan\!\left(\frac{2\pi R_a}{P_a}\right) \\
\gamma_+&=\alpha+\beta \\
\gamma_0&=\beta.
\end{aligned}\qquad\text{(15,23)}
```

The sign of ``\beta`` in ``\gamma_+`` follows the relative lay directions. The fictitious nonconductive material receives a complex permeability ``\mu^*`` chosen so that a local slab stores the same complex magnetic energy as the pitched-wire field.

**Implementation.** Build the two-dimensional conductor geometry, solve (1) with the one-step current constraints, enforce identical net current in all armor wires, determine ``\mu^*`` by the source's energy matching, and apply (6)–(14) for the requested sequence quantities. Add earth-return impedance only after the FEM reduction.

**Limitations.** The validation covers single-layer steel-wire armor with small gaps. The model approximates helical fields in two dimensions and assumes equal armor-wire currents. The source calls for further assessment for two-layer armor and wide wire spacing.

**Reference.** [Gustavsen2023](@cite), equations (1)–(30), printed pp. 3011–3018.

## Evidence and approximation sources

- Equations (1)–(14) define the FEM and terminal-impedance assembly.
- Equations (15)–(25) define the helical-field and energy-matching reduction.
- The positive- and zero-sequence cases were checked against the source's three-dimensional FEM comparisons.

