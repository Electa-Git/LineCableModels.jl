# Alvarado–Betancourt corrected complex image

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Thin overhead self or mutual conductor pair with height sum H and horizontal separation x. |
| Calculated quantities | Closed-form overhead earth-return impedance |
| Earth structure | Homogeneous conducting nonmagnetic earth below lossless air. |
| Model and approximation | Complex-image approximation with a cubic residual correction; conduction-only earth. |
| Main source | F. L. Alvarado and R. Betancourt (1983) |
| Citation key(s) | Primary: `:Alvarado1983`; equation witness: `:Papadopoulos2020` |
| Evidence status | Published comparison page image checked, equation (8) and the stated self substitution. |

**Description.** A cubic residual improves the homogeneous complex-image approximation.

**Expression.**

```math
J_m=\frac12\ln\!\left[
\frac{(1+2p_g/H)^2+(x/H)^2}{1+(x/H)^2}
\right]
-\frac1{24}\sum_{s\in\{-1,1\}}
\left[1+\frac{H}{2p_g}(1+sjx/H)\right]^{-3}.
\qquad\text{(8)}
```

```math
Z_e=\frac{j\omega\mu_0}{2\pi}\left[\ln(D/d)+J_m\right],
\qquad p_g=(j\omega\mu_0\sigma_g)^{-1/2}.
```

For mutual terms, H=h_i+h_j, D²=x²+H², and d²=x²+(h_i−h_j)². For self terms, x=0, H=2h, and d=r. The wire radius is not a lateral separation in the correction.

**Limitations.** The approximation retains neither displacement current nor magnetic-earth contrast. It is not identical to the uncorrected Dubanton image or the exact Carson integral.

**Reference.** [Alvarado1983](@cite), with the explicit equation and self prescription reproduced in [Papadopoulos2020](@cite) (Section 2.1.3, p. 4).
