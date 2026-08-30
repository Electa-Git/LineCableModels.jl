The natural model is: dielectric loss enters through each layer’s complex admittivity. You do not modify \(C\), and you do not construct an aggregate \(G\) before reducing the layers.

For layer \(k\),

\[
\kappa_k=\sigma_{\mathrm{eff},k}+s\varepsilon_k,\qquad s=j\omega
\]

\[
y_k'=\frac{2\pi\kappa_k}{\ln(b_k/a_k)}
\]

and the Engine potential coefficient is

\[
p_k=\frac{s}{y_k'}
=\frac{\ln(b_k/a_k)}{2\pi}\frac{s}{\kappa_k}.
\]

For radial layers,

\[
p_{\mathrm{gap}}=\sum_k p_k,\qquad
y_{\mathrm{gap}}'=\frac{s}{p_{\mathrm{gap}}}.
\]

Therefore,

\[
G_{\mathrm{gap}}=\Re\{y_{\mathrm{gap}}'\},\qquad
C_{\mathrm{gap}}=\frac{\Im\{y_{\mathrm{gap}}'\}}{\omega}.
\]

That is the composite aggregate. \(G\) and \(C\) are observations of the reduced complex admittance, not independently reduced material properties.

## Implemented Engine path

The complete calculation now follows one common path:

- [`admittivity`](/home/amartins/Documents/KUL/LineCableModels-engine/src/engine/admittance.jl:19) evaluates \(\kappa=\sigma+s\varepsilon\) from the frequency-evaluated material.
- [`potential_coefficient`](/home/amartins/Documents/KUL/LineCableModels-engine/src/engine/admittance.jl:51) maps each physical annulus to \(p=s/y\).
- [`admittance!`](/home/amartins/Documents/KUL/LineCableModels-engine/src/engine/admittance.jl:116) adds the annular coefficients in radial series, assembles \(\mathbf P_{\mathrm{ins}}\), then adds the earth potential coefficients.
- [`compute.jl`](/home/amartins/Documents/KUL/LineCableModels-engine/src/engine/compute.jl:114) finally calculates

\[
\mathbf Y=s\mathbf P^{-1}.
\]

Consequently, once the correct layer admittivity reaches that loop, [`G(parameters)`](/home/amartins/Documents/KUL/LineCableModels-engine/src/engine/lineparameters/lineparameters.jl:281) already returns the resulting \(\Re\{\mathbf Y\}\).

The workspace now retains each complete dielectric `Material`, including
`kind` and `tan_delta`, so constitutive relations can evaluate all material
properties at each frequency.

Likewise, DataModel’s [`dielectric_circuit`](/home/amartins/Documents/KUL/LineCableModels-engine/src/datamodel/flatten.jl:725) uses only resistivity and permittivity.

## Formula structure I recommend

Keep the existing `InsulationAdmittance` family and the parallel
`SemiconAdmittance` family. No new workflow or problem type is needed.

The Coaxial Engine owns geometry and radial aggregation. Each registered
formula defines a material constitutive relation:

```julia
material = constitutive(route, material, frequency, temperature)
κ = admittivity(material, s)
p = potential_coefficient(r_in, r_ex, material, s)
```

The built-ins should have unambiguous semantics:

- `:Gustavsen2013`

  \[
  \kappa=s\varepsilon
  \]

- `:Marti2001` — preserve the former parallel-RC behavior and numerical baseline

  \[
  \kappa=\sigma_{\mathrm{dc}}+s\varepsilon
  \]

- A future literature formula for a reported \(\tan\delta\) that already
  includes conduction:

  \[
  \kappa=\omega\varepsilon\tan\delta+s\varepsilon
  \]

  The material resistivity is not added to the shunt loss.

- A future literature formula where \(\tan\delta\) describes polarization
  loss only:

  \[
  \kappa=\sigma_{\mathrm{dc}}+\omega\varepsilon\tan\delta+s\varepsilon.
  \]

This avoids a poisonous `includes_conduction=true` flag and makes double-counting impossible to overlook.

Later dielectric-spectroscopy, Debye, Cole–Cole, or other literature formulations can be registered as `:AuthorYear`. They return the complete \(\kappa(f,T)\). Constant loss tangent should remain explicitly documented as a frequency-axis approximation, not a causal \(s\)-domain model.

For mixed XLPE, semiconductor, and bedding layers, `material.kind` selects the
`insulation_admittance` or `semicon_admittance` relation before the common
layer operator. This allows conductive semiconductors and measured-loss XLPE
within the same cable without runtime registry lookup.

This is also exactly the direction supported by the sources: Gustavsen’s improved model gives each screen its own parallel \(G/C\) branch and notes that both properties may be frequency-dependent ([Gustavsen 2001](/home/amartins/Documents/KUL/LineCableModels/bibliography/Gustavsen_2001_Panel%20session%20on%20data%20for%20modeling%20system%20transients%20insulated%20cables.pdf), pp. 719–721). The newer paper directly uses complex permittivity per physical layer and combines their admittances radially in series ([EEEIC24 paper](/home/amartins/Documents/KUL/LineCableModels/bibliography/dielectric_losses/Attenuation%20of%20traveling%20waves%20in%20HVDC%20cable%20lines%20EEEIC24%20199_cr%20%28002%29.pdf), §II-B).

## Downstream power

There are two different \(P\)’s here:

- Engine \(\mathbf P\): potential-coefficient matrix.
- \(P'_{\mathrm{diel}}\): dissipated power per unit length.

For total shunt loss using RMS phase-domain voltages,

\[
P'_{\mathrm{shunt}}
=\Re\{\mathbf V^H\mathbf Y\mathbf V\}.
\]

Using the Hermitian conductance matrix,

\[
\mathbf G_H=\frac{\mathbf Y+\mathbf Y^H}{2},
\qquad
P'_{\mathrm{shunt}}=\mathbf V^H\mathbf G_H\mathbf V.
\]

For the reciprocal matrices produced here, this agrees with `real(Y)` up to numerical noise.

Insulation-only loss cannot generally be recovered by calculating `real(s*inv(Pin))`: earth and insulation potential coefficients are coupled before inversion.

For an exact insulation breakdown, retain the radial branch dissipation. If \(I_g\) is the current crossing dielectric gap \(g\),

\[
P'_{\mathrm{diel}}
=\sum_g |I_g|^2
\sum_{k\in g}\Re\!\left\{\frac{1}{y_k'}\right\}.
\]

This can also be retained as a dielectric dissipation matrix \(\mathbf G_{\mathrm{diel}}\), so downstream code still gets

\[
P'_{\mathrm{diel}}
=\mathbf V^H\mathbf G_{\mathrm{diel}}\mathbf V
\]

without preserving every voltage vector from the solve. It should be optional result detail, computed before bundle/Kron reduction and transformed by the corresponding voltage congruence.

Finally, inferring `rho` from a 50/60 Hz \(\tan\delta\) is a legitimate reference-frequency calibration, but `:Marti2001` then freezes \(G\) while the actual constant-\(\tan\delta\) relation has \(G\propto\omega\). That inferred resistivity should therefore be advertised as a single-frequency equivalent, not a broadband material property.
