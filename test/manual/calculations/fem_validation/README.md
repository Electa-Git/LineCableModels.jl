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
  "mesh_options": {"mesh_size_factor": 2.5},
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
and native sources match. Use a fresh output directory after source changes.

Each case keeps its exported bundle, mesh and native logs, reference CSV, derived
native mesh values, qualification flags and result JSON. Errors include signed
conductance and susceptance components, diagonal-scaled Z and Y, reciprocity, DOFs and native seconds.
Numerically zero reference components have no relative component error.
Native failures and timeouts are recorded without diagnostic reruns.

For homogeneous earth, only the scratch bundle's native air conductivity,
permittivity and permeability are changed to match the earth. Native source
files remain unchanged. Homogeneous-air checks use the public lossless-earth
input. The overhead closed form uses the interface voltage reference and the
buried closed form the reference at infinity. The tool reports only admittance for the overhead closed form.

The analytical comparison uses a mean-field, single-line-source receiver and
is within its declared scope for `abs(kappa_m*r_p)` about or below 0.1.
Thick receivers near an interface, non-passive prescribed propagation, exact
cutoff and incoming-sheet references require separate interpretation.
