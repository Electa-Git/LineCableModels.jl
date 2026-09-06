# Wedepohl–Wilcox buried-cable low-order earth-return impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried filamentary cable pair; self uses cable radius and mutual terms use conductor separation and burial depths. |
| Calculated quantities | Low-order closed self and mutual earth-return impedance approximations |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The Wedepohl–Wilcox reduction truncates the Pollaczek decomposition for small complex propagation-distance products. |
| Main source | L. M. Wedepohl and D. J. Wilcox (1973) |
| Citation key(s) | Primary: `:Wedepohl1973`; comparative equation witness: `:Guneri2018` |
| Evidence status | Original publication and comparative PDF equations checked |

**Description.** Wedepohl and Wilcox retain the leading terms of their Pollaczek series for buried cables.

**Expression.** For two cables ``i`` and ``j``, source equation (8) gives

```math
\begin{aligned}
z_{ji}&=\frac{j\omega\mu}{2\pi}
\left\{
-\ln\!\left(\frac{\gamma m s_{ji}}{2}\right)
+\frac{1}{2}-\frac{2}{3}m\ell
\right\} \\
\ell&=h_i+h_j.
\end{aligned}\qquad\text{(8)}
```

The self term follows from ``s_{ji}\mapsto r_4`` and ``\ell\mapsto2h``:

```math
z_7=\frac{j\omega\mu}{2\pi}
\left\{
-\ln\!\left(\frac{\gamma m r_4}{2}\right)
+\frac{1}{2}-\frac{4}{3}mh
\right\}.
\qquad\text{(7)}
```

Here

```math
\begin{aligned}
m&=\sqrt{\frac{j\omega\mu}{\rho}} \\
s_{ji}&=\sqrt{x_{ji}^2+(h_i-h_j)^2},
\end{aligned}
```

``\rho`` is earth resistivity. For numerical evaluation, the logarithmic constant is ``\gamma=e^{\gamma_E}=1.781072\ldots``, where ``\gamma_E=0.577215\ldots`` is the Euler–Mascheroni constant. This value follows from the small-argument expansion of ``K_0`` used in Appendix 8.3.

**Implementation.** Evaluate the principal complex logarithm and select the square root consistently with the harmonic convention. Assemble each scalar mutual term into all four entries of the off-diagonal two-conductor cable submatrix, as stated in section 2.3.5.

**Limitations.** The earth is homogeneous, displacement current is neglected, the cables are replaced by thin insulated conductors in the earth-return derivation, and the source requires ``|ms_{ji}|<0.25`` or ``|mr_4|<0.25``. Use the source's full equation (40) outside this small-product region.

**Reference.** [Wedepohl1973](@cite), equations (7)–(8), printed p. 255; comparative witness [Guneri2018](@cite).

## Evidence and approximation sources

- The original page image fixes the factors ``2m\ell/3`` and ``4mh/3`` and the self substitution.
- Appendix 8.3 supplies the parent integral and the series reduction.
