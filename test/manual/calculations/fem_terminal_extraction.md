> Historical investigation record. Equations, observations and limitations are
> retained at their recorded source versions. Referenced prototypes and campaign
> launchers may have been retired; their commands and pending-work statements
> are not current execution instructions. See the [cleanup record](fem_development_cleanup.md)
> and [current PML controls](fem_fixed_pml_controls.md).

# Mixed-pair terminal extraction experiment

This records the historical averaged-path comparison. Its solve entrypoint is retired; collecting saved native outputs remains available. See the [native voltage contract](../../../docs/plans/fem-native-voltage-extraction-plan.md).

The production series impedance already uses the GetDP global quantity `U`,
printed `OnRegion Terminal~{response_terminal}`. The terminal scalar potential
likewise uses global `V`. Neither involves a contour average of the terminal
unknown. This is the global-region operation described in the
[GetDP manual](https://getdp.info/doc/texinfo/getdp.html#Types-for-PostOperation).

Averaging instead enters the aerial surface-reference samples and the vertical
vector-potential paths. Bare aerial disks use normalized contour quadrature;
buried receivers use the lowest contour node. The source quantities and the
field equations are unchanged by these postprocessing operations.

## Isolated comparisons

The retired solve mode of `run_fem_terminal_extraction.py` reused a retained PML case's mesh, material
values, native executable, source snapshots, excitations, and solver controls.
It copied the native sources to its own directory and instrumented only
postprocessing. Every solved source column yielded:

| CSV quantity | Definition |
|---|---|
| `Z` | Existing global-group series impedance |
| `P` | Existing scalar group minus aerial surface reference plus vector path |
| `Ppoint` | Same definition with one lowest-contour-node receiver and no averaging |
| `Pgroup` | Scalar global group alone, relative to the zero-Dirichlet boundary |
| `Psurface_scalar` | Scalar group minus the existing surface reference, without vector paths |
| `Y`, `Ypoint`, `Ygroup`, `Ysurface_scalar` | Matrix inverses of those four P definitions |

The separate scalar, surface and line terms are retained too. P has units
ohm m, Y has units S/m, and Z has units ohm/m. The scalar-only inverse in
quasi-fw is a diagnostic: removing the vector term changes the measured
voltage. It is not silently adopted as the physical shunt admittance.

Both formulations retain `v=0` on the entire exterior boundary, including
its air side. In quasi-fw, the transverse vector-potential boundary/gauge
constraints are also unchanged. Thus these comparisons separate extraction
effects from changing the upper boundary condition.

The native experiment uses two 0.0425 m copper disks at (0,1) and (1,-1) m,
rho=100 ohm m, and the 24-skin-depth/192-layer PML preset with mesh factors
3 and 8. The selected frequencies are 0.1 Hz, 1 kHz, 100 kHz and 1 MHz.
The wider retained-data audit covers all 99 frequencies and rho=1/100/1000.

All eight selected native solves completed. Their original Z and P entries
reproduce the retained outputs exactly, and the saved scalar/reference/vector
terms reconstruct P exactly. Constraints, function spaces, weak forms, source
definitions and resolution blocks are byte-identical to the retained sources.
The new point measurements change none of the 32 sampled conductance signs.

## Reciprocity findings

Define the relative mutual asymmetry as
`abs(A12-A21)/max(abs(A12),abs(A21))`, using complex entries and an ordinary
transpose, not a conjugate transpose. This avoids hiding a mutual difference
behind much larger diagonal entries.

Across the retained 99-frequency mixed sweeps at all three resistivities,
the largest series-Z asymmetry is 5.40e-11 for quasi-fw and 2.51e-15 for
quasi-tem. The present series matrix is numerically reciprocal in these
fixtures. The following shunt results are for rho=100 ohm m at 1 MHz:

| Formulation / extraction | Y12, microSiemens/m | Y21, microSiemens/m | Relative asymmetry |
|---|---:|---:|---:|
| Quasi-fw, existing | -11.5873 + j33.4102 | -12.4991 - j30.7671 | 1.81502 |
| Quasi-fw, one point | -11.5864 + j33.4099 | -12.4989 - j30.7672 | 1.81505 |
| Quasi-fw, group only | -0.442027 - j66.8529 | +13.7670 - j82.0975 | 0.250346 |
| Quasi-tem, existing | +1.26757 + j9.99087 | -12.2026 - j30.8292 | 1.29644 |
| Quasi-tem, one point | +1.27148 + j9.98934 | -12.2024 - j30.8293 | 1.29643 |
| Quasi-tem, group only | -12.3092 - j30.7585 | -12.3092 - j30.7585 | 2.05e-16 |
| Unified analytical reference | -11.3546 + j32.9251 | -12.3320 - j30.2916 | 1.81533 |

Across all four frequencies, removing averaging changes complex Y by at most
0.002704% for quasi-fw and 0.04163% for quasi-tem, entry by entry; both maxima
occur at 1 MHz. This does not explain
the large mixed-direction asymmetry. Componentwise differences near zero
need their own assessment; these percentages use the complex-entry magnitude.

For the scalar quasi-tem model, the group-only potential matrix is symmetric
before changing the receiver reference. Subtracting the earth-surface potential
only from the aerial receiver row destroys that matrix symmetry. This is
observed directly from the same FEM solution, independently of the previous
boundary-element calculation. The largest group-only Y asymmetry over the
four frequencies is 7.34e-16. Quasi-fw's group-only asymmetry ranges from
8.97e-4 to 0.3602 at the same samples.

The unified reference itself has a strongly asymmetric mixed Y. Its documented
voltage convention uses a surface reference for aerial receivers and deep earth
for buried receivers. Agreement with that reference and symmetry of the
published Y are therefore different requirements in these retained fixtures.
This observation does not establish a physically nonreciprocal medium.

An explicit hypothesis is that the voltage and source/current coordinates are
not dual port coordinates after the receiver-dependent change of reference.
In the scalar experiment, if C maps absolute electrode potentials to currents
and M maps absolute potentials to the chosen receiver voltages, the measured
matrix is `C*inv(M)`, which need not be symmetric even when C is symmetric.
The corresponding bilinear reciprocal-coordinate change would transform the
currents too, yielding `transpose(inv(M))*C*inv(M)`. This is an explanation to
investigate, not a proposed replacement matrix or an instruction to force
symmetry. The terminal voltage/current and return-path definitions must be
settled before making such a change.

Quasi-fw couples scalar and vector potentials. Group-only extraction does not
restore symmetry, and its retained 99-frequency curves have strong jumps as
the meshes/gauge trees change. A direct scalar coefficient is insufficient to
validate the physical voltage. The ECE formulation paper defines terminal
voltages through the electric field and treats their reduction to scalar
terminal values under specific boundary/gauge conditions; those conditions
must be checked for this 2D reduction. See
[Ciuprina and Sabariego, 2024](https://link.springer.com/article/10.1186/s13362-024-00165-6),
sections 2 and 4. No alternative gauge was solved in this experiment.

The pure group extraction also fails to preserve the successful aerial and
buried quasi-fw results. At rho=100 ohm m and 1 MHz:

| Placement | Existing G11, S/m | Group-only G11, S/m | Analytical G11, S/m |
|---|---:|---:|---:|
| Both aerial | -2.92198e-7 | +1.00924e-6 | -2.88927e-7 |
| Both buried | +1.35419e-2 | +1.03819e-2 | +1.34107e-2 |

These checks use the already-retained group scalars with the current exterior
Dirichlet condition. They do not require another PDE solve.

## Series self impedances

Changing these P/Y measurements leaves Z unchanged. For each mixed receiver,
comparing its self Z against the same-height receiver in the corresponding
both-aerial or both-buried system gives the following maximum complex relative
differences over the 99-frequency quasi-fw sweeps:

| Soil resistivity, ohm m | Maximum self-Z change |
|---:|---:|
| 1 | 0.27245% |
| 100 | 0.08727% |
| 1000 | 0.06921% |

These are measured differences, not a convergence bound. The unified formula
also has small diagonal changes through its complete current-closure matrix;
its maximum corresponding changes are 0.01127%, 0.005163%, and 0.001082%.
The remaining FEM placement sensitivity needs a local geometry/mesh convergence
comparison with unchanged exterior controls. It is not evidence that the
admittance averaging changed the series equation.

A separate, verified low-frequency error is the coarse circular geometry.
At rho=100 ohm m and 0.1 Hz, each 4.25 cm-radius core has meshed area
0.00541875 m², versus pi*r²=0.00567450173055 m²: a 4.5070% area deficit from
its 12-sided polygon. With copper rho=1.7241e-8 ohm m, rho/area is
3.1817301e-6 ohm/m on this polygon and 3.0383284e-6 ohm/m on the circle.
The predicted resistance increase, 1.4340166e-7 ohm/m, closely accounts for
the observed self-R difference of 1.4316294e-7 ohm/m against the analytical
reference. This coarse-preset error occurs in all placements; it does not
explain a mixed-only effect. Preserve the PML while refining the conductor
contour/interior for the next series-accuracy test.

## Reproduction and artifacts

The preparation-dependent solve entrypoint has been retired. To rebuild the
summary of an existing study from its saved native outputs:

```sh
python3 test/manual/calculations/run_fem_terminal_extraction.py --collect \
  /tmp/fem-exterior-grading/dense/quasi_fw-air_earth-rho100 \
  /path/to/existing/extraction-study
```

Replace `quasi_fw` with `quasi_tem` for a saved scalar study. Collection reads
the existing raw TSV and command records and rewrites only the summary CSV.
It requires NumPy, without Gmsh or a solver. Ordinary plotting reads the CSVs;
the historical averaged and single-point labels are preserved. Reproducing the
retired solve requires its original source snapshot, not the current backend.

Evidence is under `/tmp/fem-terminal-extraction`, with a lightweight retained
copy under `.linecablemodels/fem/terminal-extraction-evidence`. It includes
signed CSVs, spectrum plots, reference comparisons and the analysis script.
The `plots/quasi_fw-G.png` and `plots/quasi_tem-G.png` files compare the
four-frequency extraction choices; `plots/self-Z-placement.png` shows the
99-frequency mixed-versus-same-layer self-Z differences. The
`group-spectrum/plots` subdirectory compares all 99 retained quasi-fw
frequencies with the scalar-only diagnostic at all three resistivities.
The initial study output used literal backslash-t separators; only those
diagnostic TSVs were corrected and recollected. Native Z/P outputs and solves
were preserved, and subsequent runs use actual tab separators.

No production source, boundary condition, detached-export runtime, or manual
two-wire runner was modified for this experiment.
