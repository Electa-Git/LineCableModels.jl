# Ametani classical transmission-line potential reference

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Thin overhead conductors or buried coaxial units with separately assembled radial insulation. |
| Calculated quantities | Constant overhead space potential; zero exterior potential for buried classical TL units |
| Earth structure | Equipotential reference plane; no finite-conductivity potential correction. |
| Model and approximation | Classical TL reference, not a lossy-earth potential solution. |
| Main source | A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad (2021) |
| Citation key(s) | `:Ametani2021` |
| Evidence status | IET book page image checked, Table 2.1 and equations (2.76)–(2.78), p. 24. |

**Description.** The classical TL comparison combines an earth-return series impedance with constant overhead space potential or buried-cable insulation potential.

**Expression.**

```math
P_{e,ij}^{00}=\frac{1}{2\pi\varepsilon_0}\ln(D_{ij}/d_{ij}),
\qquad P_{e,ij}^{11}=P_{e,ij}^{01}=0.
```

```math
\mathbf P=\mathbf P_i+\mathbf P_e,
\qquad \mathbf Y=j\omega\mathbf P^{-1}.
```

The self overhead distances are D=2h and d=r. Buried insulation potential remains in P_i; zero P_e does not mean zero total cable capacitance.

**Limitations.** The zero mixed correction is a block-separated classical reference, not evidence that physical mixed-earth coupling vanishes. Select a mixed-potential formulation for that coupling. This boundary convention is not attributed to Pollaczek's 1926 impedance paper.

**Reference.** [Ametani2021](@cite) (Section 2.6.5, Table 2.1).
