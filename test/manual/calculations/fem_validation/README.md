# FEM validation

Run current native FEM rules against the `:unified` reference for overhead,
buried and mixed pairs, or against homogeneous solid-cylinder closed forms.
The tool does not refine meshes or change production sources.

```sh
julia --project=test test/manual/calculations/fem_validation/run.jl /tmp/fem-validation all
```

Groups are `main`, `homogeneous` and `all`. The main grid has 18 frequencies
from 1 mHz to 300 MHz, three placements and four earth resistivities
(0.1, 1, 100 and 1000 Ω·m). Each conductor has radius 0.0425 m and the centres
are 2 m apart, at heights or depths of 1 m. Homogeneous checks include overhead
and buried single cylinders in air with prescribed longitudinal propagation,
and buried cylinders in dielectric earth. Exact cutoff is excluded from the
closed forms because it is singular.

An optional third argument is a JSON settings file:

```json
{
  "options": {"mesh_size_factor": 2.5},
  "frequencies": [50, 10000],
  "earth": [{"rho": 100, "eps_r": 12}],
  "timeout_seconds": 2700,
  "max_dofs": 4000000,
  "list_only": false
}
```

Without overrides all mesh controls use current production defaults. Runs are
sequential. GetDP preprocessing records the actual DOF count and skips solves
above `max_dofs`. For independent processes sharing an output directory,
`case_indices` selects one-based entries from the full case list and `record_file`
names each process's result ledger; keep the selections disjoint. `list_only` emits the requested cases without references, meshes
or solves. Saved results are reused only when the exported numerical inputs
and native sources, comparison formulas and budget calculations match. Use a
fresh output directory after source changes.

Each case keeps its exported bundle, mesh and native logs, reference CSV, derived
native mesh values, qualification flags and result JSON. Errors include signed
conductance and susceptance components, diagonal-scaled Z and Y, reciprocity, DOFs and native seconds.
Numerically zero reference components have no relative component error.
Native failures and timeouts are recorded without diagnostic reruns.

For homogeneous earth, only the scratch bundle's native air conductivity,
permittivity and permeability are changed to match the earth. Native source
files remain unchanged. Homogeneous-air checks use the public lossless-earth
input. The overhead closed form uses the interface voltage reference and the
buried closed form the reference at infinity. The overhead impedance reference uses the voltage-reference correction
`Gamma^2 * (1/Y_h - 1/Y_deep)` from the native impedance chain.

The analytical comparison uses a mean-field, single-line-source receiver and
is within its declared scope for `abs(kappa_m*r_p)` about or below 0.1.
Thick receivers near an interface, non-passive prescribed propagation, exact
cutoff and incoming-sheet references require separate interpretation.


Engineering error budgets use `e_Z = max(abs(Z-Z_ref))/max(abs(diag(Z_ref)))`
and the corresponding `e_Y`. For conductance entries at least 0.001 times the
largest reference admittance diagonal, the conductance gate is a 5% relative
error; smaller entries remain reported without a gate. Sign canaries use a
0.000001 diagonal-scale threshold. All-entry sign agreement is also retained.
A reciprocity excess above the reference of 0.3 percentage points is a canary.
Nonzero prescribed propagation with `imag(Gamma) < k0`, where
`k0 = 2*pi*f*sqrt(mu_air*epsilon_air)`, is a fast-wave robustness check,
reported without an accuracy gate. Zero prescribed propagation remains in scope.
Receiver transverse sizes above 0.1 and exact transverse cutoff are reported
outside the reference scope, rather than used for acceptance. Budget fields
are measurements; selecting a candidate remains a separate comparison against
the production baseline. The homogeneous checks include both impedance and admittance references, using
the exported conductor material values for the internal impedance.

Native observations retain attenuation exponents, directional absorption, interval counts, extent-cap and cutoff flags, and earth sizing/layer values. Exact transverse cutoff warns that the transverse problem is singular. Only the earth-resistive sizing ceiling marks results unqualified; below-target PML attenuation does not emit a qualification warning.
