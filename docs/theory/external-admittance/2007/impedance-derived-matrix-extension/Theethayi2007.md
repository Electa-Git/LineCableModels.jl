# Theethayi impedance-derived potential recipe

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Buried self or horizontally separated mutual conductor pair. |
| Calculated quantities | Potential coefficients from the selected logarithmic-exponential impedance |
| Earth structure | Homogeneous conducting nonmagnetic earth below air. |
| Model and approximation | Vance scalar conversion with the Theethayi impedance approximation; matrix extension follows the IET book. |
| Main source | N. Theethayi, R. Thottappillil, M. Paolone, C. A. Nucci, and F. Rachidi (2007) |
| Citation key(s) | `:Theethayi2007`; scalar attribution: `:Vance1978`; matrix extension: `:Ametani2021` |
| Evidence status | Scalar source equations (9)–(10) and IET book Section 2.6.3 checked. |

**Description.** This recipe combines an existing impedance approximation with the Vance conversion. It is not another independent material law.

**Expression.**

```math
P_{e,ij}=\frac{j\omega Z_{e,ij}}{\gamma_1^2},
\qquad \gamma_1^2=j\omega\mu_0(\sigma_1+j\omega\varepsilon_1).
```

At the matrix level, the selected extension gives Y_e=γ₁²Z_e⁻¹. The implementation first assembles P_e; it does not take elementwise reciprocals to produce Y_e. Internal dielectric terms are assembled separately before the total matrix inverse.

**Limitations.** The original scalar relation alone does not establish mutual coupling. The matrix extension is explicitly discussed in the later book and inherits the selected impedance approximation and its horizontal-spacing restriction. It is distinct from the book's further Xue expression, which uses a complex exponential instead of the Theethayi attenuation magnitude.

**Reference.** [Theethayi2007](@cite) (equations (9)–(10)), [Vance1978](@cite), and [Ametani2021](@cite) (Section 2.6.3, equation (2.72)).
