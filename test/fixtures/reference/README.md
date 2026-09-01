# Numerical references

Files in this directory contain static numerical references calculated
independently of LineCableModels. Each record states its units, input values,
working precision, source, and generation method. Tests read these files but
do not update them.

`coaxial_capacitance.toml` evaluates the analytical capacitance per unit length of a
coaxial geometry at 256-bit precision using the CODATA 2018 vacuum permittivity. The
calculation was performed without loading LineCableModels. The record tests
cross-precision convergence of the generic dielectric kernel and records the
corresponding 50 Hz lossless impedance and potential-coefficient values.
The same record includes the independent 256-bit tubular-inductance value for the
documented 4π×10⁻⁷ H/m permeability convention.
The record also includes the conductance of the same geometry at 2×10¹¹ Ω·m for the
Ametani2004 complex-permittivity formulation.

LineCableModels regression output is not stored here because it would repeat
the implementation under test instead of supplying an independent reference.

The exception is `cable_model_v1_preservation.json`. It records the complete
cable-design, line-system, normalized native-input, and numerical-result
boundary at commit `36dfe57d330eef8dcdb00f35b905954be84401e5`. It is a
deliberate refactor-preservation lock, not an independent physics reference,
and must not be rewritten to accommodate a changed engine payload.
