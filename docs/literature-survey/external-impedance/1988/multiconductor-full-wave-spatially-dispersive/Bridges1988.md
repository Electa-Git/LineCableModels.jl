# Bridges multiconductor spatially dispersive external impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | ``N`` thin circular conductors above dissipative earth. |
| Calculated quantities | Spatially dispersive self and mutual external-impedance matrix, assembled with the conductor surface-impedance matrix |
| Earth structure | Homogeneous lossy earth below air. |
| Model and approximation | The longitudinal spectral variable is retained. Circular averaging gives Bessel factors and the earth contribution remains a Sommerfeld integral; modal roots are calculated only after the impedance matrix is assembled. |
| Main source | G. Bridges, O. Aboul-Atta, and L. Shafai (1988) |
| Citation key(s) | Primary: `:Bridges1988`; precursor: `:Wait1977` |
| Evidence status | Original publication equations checked; the matrix is an explicit p.u.l. formulation |

**Description.** The formulation retains the longitudinal wavenumber ``k_z`` and assembles a full external-impedance matrix for multiple overhead conductors. Modal propagation constants are zeros of the assembled matrix determinant.

**Expression.** For equal air and earth permeabilities, source equation (7) gives

```math
Z^e_{mn}(k_z)=A_n\left\{
\tau_e^2\left[K_0(\tau_e\rho_{mn})-K_0(\tau_e\rho^*_{mn})\right]
-k_e^2J(\tau_e,\rho^*_{mn})
+k_z^2G(\tau_e,\rho^*_{mn})
\right\},
\tag{7a}

A_n=\left(-\frac{j\omega\mu_e}{2\pi k_e^2}\right)
\frac{1}{(\tau_ea_n)K_1(\tau_ea_n)},
\tag{7b}

J=\int_{-\infty}^{\infty}
\frac{e^{j\lambda|x_m-x_n|-U_e(y_m+y_n)}}{U_e+U_g}\,d\lambda,
\qquad
G=\int_{-\infty}^{\infty}
\frac{e^{j\lambda|x_m-x_n|-U_e(y_m+y_n)}}{n^2U_e+U_g}\,d\lambda,
\tag{7c}

U_e=\sqrt{\lambda^2+\tau_e^2},
\quad
U_g=\sqrt{\lambda^2+\tau_g^2},
\quad
\tau_e^2=k_z^2-k_e^2,
\quad
\tau_g^2=k_z^2-k_g^2,
```

with ``\operatorname{Re}U_e,\operatorname{Re}U_g,\operatorname{Re}\tau_e,\operatorname{Re}\tau_g\ge0`` and

```math
\rho_{mn}=\sqrt{(x_m-x_n)^2+(y_m-y_n)^2},
\qquad
\rho^*_{mn}=\sqrt{(x_m-x_n)^2+(y_m+y_n)^2}.
```

The spectral system uses

```math
\mathbf Z(k_z)=\mathbf Z^w(k_z)-\mathbf Z^e(k_z),
\qquad
\det\!\left[\mathbf Z^w(k_z)-\mathbf Z^e(k_z)\right]=0.
\tag{5,10}
```

``\mathbf Z^w`` is the conductor surface-impedance matrix. Equation (8) gives its diagonal solid-conductor form.

**Implementation.** For each trial ``k_z``, evaluate every matrix element in (7), add the selected conductor surface term through (5), and solve (10). The Bessel factor in ``A_n`` represents circumferential averaging over conductor ``n``.

**Limitations.** Each conductor radius must be small compared with the free-space wavelength and all interconductor and image distances. The currents are axial and azimuthally uniform. The earth is one homogeneous half-space.

**Reference.** [Bridges1988](@cite), equations (1)–(10), printed pp. 429–430; precursor [Wait1977](@cite).

## Evidence and approximation sources

- Equation (7) and its geometric definitions were transcribed from the original page image.
- Equations (5), (8), and (10) define the matrix assembly, solid-conductor surface term, and modal root.

