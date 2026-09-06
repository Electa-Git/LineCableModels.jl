# Petrache infinite-earth logarithmic impedance approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried cable in an unbounded conducting medium; self uses cable radius and mutual terms use conductor distance. |
| Calculated quantities | Logarithmic approximation for self and mutual earth-return impedance |
| Earth structure | Infinite homogeneous earth. |
| Model and approximation | A closed logarithm in the complex earth propagation-distance product approximates the infinite-earth return term. |
| Main source | Marian Petrache, Florin Rachidi, Mario Paolone, Carlo Alberto Nucci, Vladimir A. Rakov, and M. A. Uman (2005) |
| Citation key(s) | Primary attribution: `:Petrache2005`; equation witness: `:Guneri2018` |
| Evidence status | Equation checked in the accessible comparative publication; publication identity and DOI verified |

**Description.** The approximation gives self and mutual earth-return impedances as elementary logarithms of the full lossy-earth propagation-distance product.

**Expression.** The mutual and self terms are

```math
Z_{g,m}=\frac{j\omega\mu_0}{2\pi}
\ln\!\left(\frac{1+\gamma_gd}{\gamma_gd}\right),
\tag{P1}

Z_{g,s}=\frac{j\omega\mu_0}{2\pi}
\ln\!\left(\frac{1+\gamma_gR}{\gamma_gR}\right),
\tag{P2}

\gamma_g=\sqrt{j\omega\mu(\sigma+j\omega\epsilon)}.
\tag{P3}
```

Use conductor distance ``d`` for a mutual element and outer cable radius ``R`` for a self element.

**Implementation.** Select the decay-consistent square root for ``\gamma_g`` and a continuous complex-logarithm branch over the frequency grid. Apply (P1) or (P2) directly.

**Limitations.** The formula does not contain burial depth or an explicit air–earth interface. It is an infinite-homogeneous-earth approximation and should remain separate from the source's half-space cable assembly.

**Reference.** [Petrache2005](@cite); equation witness [Guneri2018](@cite), section 3.4.

## Evidence and approximation sources

- The accessible comparative publication prints the mutual and self forms with their material definition.
- The primary publication identity and DOI establish the attribution used by the index.

