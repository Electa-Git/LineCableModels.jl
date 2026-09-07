# Wedepohl–Wilcox hollow-shell approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Homogeneous circular conducting annulus with inner radius a and outer radius b. |
| Calculated quantities | Inner, outer, and transfer surface impedances |
| Earth structure | Not applicable. |
| Model and approximation | Hyperbolic approximation to the cylindrical surface-impedance relations; not an exact low-frequency formula. |
| Main source | L. M. Wedepohl and D. J. Wilcox (1973) |
| Citation key(s) | Primary: `:Wedepohl1973`; equation witness: `:Ametani2021` |
| Evidence status | IET book page images checked, equations (A1.56)–(A1.58). |

**Description.** The shell approximation retains a curvature correction in the two surface terms and an arithmetic-mean radius in the transfer term.

**Expression.**

```math
Z_{iw}=\frac{\rho m}{2\pi a}\coth[m(b-a)]-\frac{\rho}{2\pi b(a+b)}.
\qquad\text{(A1.56)}
```

```math
Z_{mw}=\frac{\rho m}{\pi(a+b)\sinh[m(b-a)]}.
\qquad\text{(A1.57)}
```

```math
Z_{iw}=\frac{\rho m}{2\pi b}\coth[m(b-a)]+\frac{\rho}{2\pi b(a+b)}.
\qquad\text{(A1.58)}
```

Here m²=jωμ/ρ, with the decaying root. The book prints Z_iw again in (A1.58); this is the outer-surface term, as specified by its radius and the surrounding derivation.

**Limitations.** The shell must be homogeneous, circular, and hollow. The retained approximation is not equivalent to Schelkunoff's full Bessel ratios or to the distinct Zhao approximation. It supplies no core-proximity correction.

**Reference.** [Wedepohl1973](@cite), as reproduced in [Ametani2021](@cite) (Appendix A1.4.2, pp. 67–68).

**Numerical interpretation.** Scaled exponentials evaluate the hyperbolic factors without overflow. The zero-frequency values are the limits of the approximate expressions, not a replacement by exact annular dc resistance.
