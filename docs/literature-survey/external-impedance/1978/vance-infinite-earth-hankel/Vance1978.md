# Vance infinite-earth Hankel-ratio impedance approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried cable in an unbounded conducting medium; self uses cable radius and mutual terms use conductor distance. |
| Calculated quantities | Hankel-ratio approximation for self and mutual earth-return impedance |
| Earth structure | Infinite homogeneous earth. |
| Model and approximation | The result uses the ratio of first-kind Hankel functions at the complex earth propagation-distance product and does not represent the air–earth boundary separately. |
| Main source | E. F. Vance (1978) |
| Citation key(s) | Attributed source: `:Vance1978`; equation witness: `:Guneri2018` |
| Evidence status | Equation checked in the accessible comparative publication; original book attribution retained |

**Description.** Vance's approximation treats the cable as embedded in an unbounded lossy medium. The self and mutual forms differ only by their radial distance.

**Expression.** The comparative source prints

```math
Z_{g,m}=\frac{\omega\mu_0}{2\pi\gamma_gd}
\frac{H_0^{(1)}(j\gamma_gd)}{H_1^{(1)}(j\gamma_gd)},
\tag{V1}

Z_{g,s}=\frac{\omega\mu_0}{2\pi\gamma_gR}
\frac{H_0^{(1)}(j\gamma_gR)}{H_1^{(1)}(j\gamma_gR)},
\tag{V2}

\gamma_g=\sqrt{j\omega\mu(\sigma+j\omega\epsilon)}.
\tag{V3}
```

``H_0^{(1)}`` and ``H_1^{(1)}`` are first-kind Hankel functions. Use ``d`` for mutual separation and ``R`` for the self radius.

**Implementation.** Calculate ``\gamma_g`` with a decay-consistent square-root branch, evaluate the Hankel ratio, and apply the mutual or self radial argument. A scaled special-function implementation avoids overflow for large complex arguments.

**Limitations.** The approximation omits burial depth and the air–earth interface. It applies when the surrounding-earth approximation is adequate; the index keeps it separate from half-space formulas.

**Reference.** Attribution to [Vance1978](@cite); equation witness [Guneri2018](@cite), section 3.3.

## Evidence and approximation sources

- The accessible comparative publication supplies both scalar equations, the full lossy-earth propagation constant, and the source attribution.
- The record does not infer an air–earth image term that the formula does not contain.

