# Internal impedance of conductors

Internal impedance describes longitudinal current diffusion and magnetic
energy within metal, including skin effect and, where the selected model
provides it, proximity effect. Surface and transfer coefficients populate
the metallic terms of the [matrix formulation](matrix_formulation.md);
pipe-wall coefficients also depend on the positions of the enclosed
conductors.

## Cylindrical and equivalent-section conductors

- [Default cylindrical-conductor surface impedances](internal-impedance/default.md)
  (`:schelkunoff1934` is the explicit literature identifier for the same route)

## Pipe walls and conductor proximity

- [Default analytical pipe-type treatment](internal-impedance/pipe-default.md)

## Field and constitutive formulations

These records describe discretized fields or additional constitutive inputs.
Their presence in the theory collection does not imply support by the
coaxial backend.


[Back to Contents](contents.md)
