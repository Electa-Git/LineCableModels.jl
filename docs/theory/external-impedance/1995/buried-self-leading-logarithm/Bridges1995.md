# Bridges buried-wire leading logarithm

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | One thin bare or concentrically insulated buried wire, sampled at its external radius r. |
| Calculated quantities | Approximate external self impedance |
| Earth structure | Homogeneous conducting nonmagnetic earth. |
| Model and approximation | Leading radial logarithm, with displacement current and the explicit interface-depth correction omitted. |
| Main source | G. E. J. Bridges (1995), as reproduced by Ametani et al. (2021) |
| Citation key(s) | Attributed source: `:Bridges1995`; equation witness: `:Ametani2021` |
| Evidence status | IET book page image checked, equation (2.51). |

**Description.** The retained TL reduction is a self coefficient, not the full plane-wave scattering problem of the cited paper.

**Expression.**

```math
Z_e=-\frac{j\omega\mu_0}{2\pi}
\ln\!\left(\frac{e_c\gamma_1r}{2}\right),
\qquad e_c=1.7811,\quad \gamma_1^2=j\omega\mu_0\sigma_1.
\qquad\text{(2.51)}
```

**Limitations.** The logarithm is a small-propagation-radius approximation. No mutual term or independently depth-dependent correction is supplied. The original paper's antenna-excitation setting does not prevent this explicitly reproduced TL reduction from supplying an external self coefficient.

**Reference.** [Bridges1995](@cite), as reproduced in [Ametani2021](@cite) (Section 2.5.3.3, p. 18).

**Numerical interpretation.** The printed rounded factor 1.7811 is retained. This is not the complete Wedepohl expression, which also retains constant and depth-dependent terms.
