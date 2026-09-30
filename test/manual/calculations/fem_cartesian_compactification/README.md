# Cartesian compactification and absorption experiment

Historical experiment: all executable prototypes, including the independent
wave controls, were archived during the [development cleanup](../fem_development_cleanup.md).
The equations and results below are preserved as research notes; their commands
are not current entry points. The cable comparison used the retired
path-preparation implementation. See the [native voltage contract](../../../../docs/plans/fem-native-voltage-extraction-plan.md).

This is an isolated native GetDP experiment. It does not change the production
backend, the two-wire manual runner, its running processes, or retained runs.
The aim is to test a smaller computational exterior while preserving raw
conductance signs. Results live in `/tmp/fem-cartesian-compactification/`.

## Mathematical distinction

GetDP supports a real rectangular infinite-shell Jacobian, `VolRectShell`.
Its outer mesh boundary already represents physical infinity; attaching a
second physical layer *after* that endpoint is not a valid coordinate map.
Two valid constructions are a real map stopped at a finite distance followed
by a PML, or a single composed real/complex map tending to complex infinity.
This first experiment tests the latter.

For each outward coordinate, let u=d/T be the normalized mesh distance through
the shell. In the physical interior the map is the identity. In the shell use

```
s(u) = 1 + (1-i) A u^3/(1-u),                  0 <= u < 1
F(u) = T [u + (1-i) A H(u)],
H(u) = -log(1-u) - u - u^2/2 - u^3/3.
```

Then dF/dd=s, s(0)=1, and F tends to complex infinity at the outer mesh edge.
For exp(+i omega t), the negative imaginary coordinate damps outgoing waves.
The positive real part extends the represented distance for evanescent fields.
The coefficient has a simple pole rather than the stronger pole produced by
putting a polynomial absorber after a rational map. Discrete endpoint and corner
behavior still require convergence checks; analytic matching is not a finite
element error guarantee.

`buffered.pro` tests a second construction, closer to a real shell followed
by a PML. With b=1/2, w=max((u-b)/(1-b),0), and Q=12:

```
C = [(Q-1)*L/T-b]/H(b)
F(u) = T [u + C H(u) - i A H(w)]
s(u) = 1 + C u^3/(1-u) - i A/(1-b) w^3/(1-w).
```

Before u=b the map is purely real. At u=b its represented physical
half-width is Q*L; absorption begins there, before the singular endpoint.
Thus the original two-skin-depth mesh represents an absorption onset at
24 skin depths when the layout floor does not dominate. The real and imaginary
maps still both extend to infinity at u=1. The entire composed Jacobian is
included in the tensors below. Here Q=12 is a test control, not an optimized
or validated default. Compressing additional real wave periods still requires
enough elements to resolve them.

With J=diag(sx,sy,1), D=sx*sy, use the complete pullback

```
T = D J^-1 J^-T = diag(sy/sx, sx/sy, sx*sy)
mu_new = mu*T; kappa_new = (sigma+i*omega*epsilon)*T
nu_new = nu*inverse(T).
```

The axial mass terms, cross terms and vector-field transformations use this
same map. Ordinary `Vol` Jacobians remain in place. Adding `VolRectShell` as
well would apply another real transformation. In particular, GetDP 3.5's
`VolRectShell` is a directional radial scaling with off-diagonal Jacobian
terms, not the separable diagonal map tested here. Its regions and corner
partition cannot be interchanged with this tensor without a fresh derivation.

The x map is identical above and below y=0. The ground interface stays at y=0
and retains material transmission. The outer boundary retains homogeneous
Dirichlet constraints. Gauss points lie inside u<1; the prototype does not
evaluate field maps at the singular exterior boundary. Aerial voltage paths
stay in the unchanged physical region. Buried reference paths are not tested.

## Controls and scope

The dimensionless strengths are A=alpha/(k*T), using the air phase wavenumber
for the shared side and upper stretch, and the earth phase wavenumber below.
Alpha=1 and 4 are explicit experimental choices. For a normally outgoing air
plane wave, the magnitude decays asymptotically as (1-u)^alpha. There is no
finite normal-reflection-target claim for this infinite profile.

`run_waves.py` compares the existing finite cubic PML with this profile at
64/128 layers on identical saved meshes. It uses independent complex K0 fields
in the physical region, including corner probes, for low/high-frequency air
and a conducting medium with a deliberately strong shared side stretch.
Maxwell TM, Maxwell TE and the scalar operator are checked separately.

The retired `run_cable.py` copied the immutable GetDP sources of the earlier completed
two-wire aerial run `run-wjUAdy`. It replaces only the copied `pml.pro`, reuses
the original small-domain meshes and voltage paths, and runs both terminal
excitations. The default checks are rho=0.1 ohm m, r=0.0425 m, f=0.1 Hz and
4641.588833612777 Hz, quasi-fw. The small-domain controls are two skin depths,
128 PML layers, and physical mesh factor 1. Raw P is read from native TSV and
Y=inv(P); signed G and native cost are written to CSV. These output directories
are not represented as completed public API runs.

The independent wave script runs one GetDP process at a time with one solver
thread. Its controls remain available. The cable launcher depended on generated
voltage-path files and has been removed; the recorded cable results retain
their original source and measurement convention.

```sh
export LINECABLEMODELS_GETDP=/path/to/getdp
python3 test/manual/calculations/fem_cartesian_compactification/run_waves.py
```

Python needs numpy/scipy. `COMPACT_MESH_ROOT` selects a directory containing
the existing cylindrical benchmark `meshes.csv`; `COMPACT_OUTPUT` selects wave
output. The historical cable artifacts record their immutable source run,
output directory, alpha values and frequency indices. Native commands were captured as shell
text, field errors and conductances as CSV, and solver output as plain logs.
For separate buffered wave checks, set `COMPACT_PROFILES=buffered`,
`COMPACT_ALPHAS=4`, `COMPACT_LAYERS=128`, and select a separate `COMPACT_OUTPUT`.

`--split-shell` on the cable script, or `COMPACT_SPLIT_SHELL=1` on the wave
script, redistributes the existing 128 normal intervals into 64 before and
64 after absorption starts. The old geometric progression left only 13
intervals after that onset. This control changes only exterior node coordinates;
connectivity, node count and all physical-region coordinates are preserved.
It writes the modified mesh and a node-count record into the experiment's own
output directory. `COMPACT_CASES=air-low` selects the low-frequency wave check.

The redistribution is a controlled diagnostic, not a general exterior mesh
design: in the high-frequency wave control the largest phase step in the
purely real buffer is 4.63 radians on the old mesh and 8.26 radians on the
split mesh. Resolving absorption onset alone can therefore worsen resolution
of real propagating waves. A subsequent mesh design must account for the
represented physical phase and decay, including where the field is still
appreciable, rather than mesh-coordinate distance alone.

## Results, 2026-09-26

All planned native checks completed: 18 direct-profile wave cases, three
buffered wave cases, one redistributed low-frequency wave case, and seven
coupled cable case/profile/frequency checks. Each wave case compares three
operators at eleven physical probes (726 comparisons). None of the 28 raw
cable conductance entries acquired the required negative reference sign.
These candidates are therefore **not replacements for the working preset**.

The direct profile's maximum physical-probe relative wave errors were:

| Profile | Layers | Low-frequency air | High-frequency air | Conducting medium |
| --- | ---: | ---: | ---: | ---: |
| Existing finite cubic PML | 64 | 0.3079% | 0.4257% | 0.4585% |
| Existing finite cubic PML | 128 | 0.0440% | 0.3475% | 0.4542% |
| Compact, alpha=1 | 64 | 0.6537% | 2.2676% | 0.4568% |
| Compact, alpha=1 | 128 | 0.1514% | 1.1163% | 0.4539% |
| Compact, alpha=4 | 64 | 0.3900% | 0.5047% | 0.4574% |
| Compact, alpha=4 | 128 | 0.0665% | 0.3319% | 0.4540% |

Thus ordinary wave-field agreement does not establish resolution of the tiny
real part of the coupled terminal admittance. In the rho=0.1 ohm m, 0.1 Hz
two-wire case, where the analytical G11/G12 are -1.578139e-24/-5.953726e-25 S/m:

| Exterior | G11 [S/m] | G12 [S/m] | Native time | Native peak memory |
| --- | ---: | ---: | ---: | ---: |
| Earlier small-domain finite PML | +1.025312e-22 | +1.035091e-22 | prior run | prior run |
| Compact, alpha=1 | +6.609006e-23 | +6.706660e-23 | 108 s | 2218 MB |
| Compact, alpha=4 | +6.949722e-23 | +7.047390e-23 | 139 s | 2238 MB |
| Buffered, old shell grid | +1.900419e-21 | +1.901466e-21 | 275 s | 2223 MB |
| Buffered, redistributed shell | +6.423325e-24 | +7.397761e-24 | 361 s | 2224 MB |
| Working enlarged-domain preset | -1.339862e-24 | -3.570835e-25 | 539 s | 10516 MB |

The compactified cable meshes have 463,593 DOFs at 0.1 Hz. Redistribution keeps
all 117,026 nodes and 10,018 physical-region node coordinates, moving 107,008
exterior nodes. It reduces the buffered low-frequency wave error from 2.14%
to 0.348%, and changes the cable G error substantially, but still fails its sign.
The working preset row is the earlier one-thread measurement, not a rerun or
a full performance comparison of the two implementations.

At 4641.588833612777 Hz the direct alpha=1/4 G11 values are
+1.551217e-13/+1.444178e-13 S/m; the buffered result is +1.819195e-13.
The reference is -2.013703e-14. All four entries fail the sign check in each
case. The redistributed mesh was tested only at 0.1 Hz, not this frequency.

The buffered high-frequency air wave error reaches 393% on the old grid.
That grid contains an undamped phase step of 4.63 radians per element in the
real buffer. The 64/64 redistribution improves absorption resolution but would
increase this maximum step to 8.26 radians. Its high-frequency wave case was
not run, and no broadband success is inferred from its low-frequency change.

The measured evidence supports further work on the *represented-coordinate*
mesh, resolving both real wave phase before absorption and decay after its
onset. It does not support simply adding `VolRectShell` to the current PML,
or treating an infinite analytic transform as automatically resolved by a
small fixed mesh. The next candidate needs a deliberate allocation of nodes
in represented physical/complex distance, rather than another blind strength
adjustment. No change to the user's running/manual production setup was made.

Raw evidence under `/tmp/fem-cartesian-compactification/`:

- `waves.csv`: direct and finite-profile wave comparisons.
- `buffered-waves/waves.csv` and `buffered-split-waves/waves.csv`: staged controls.
- `cable/conductance.csv`, `buffered-cable/conductance.csv`, and
  `buffered-split-cable/conductance.csv`: all four raw G/P entries and costs.
- Per-case directories retain commands, solver sources, native logs and fields.

Analytic map derivatives agreed with finite differences to 9.1e-9 relative
at the sampled interior points; the determinant/tensor pullback identity was
also checked. Python syntax checks passed. These are structural checks, not
validation of the proposed exterior discretization.

## Sources

- [GetDP Jacobian documentation](https://getdp.info/doc/texinfo/getdp.html#Types-for-Jacobian).
- [Henrotte's directional rectangular-shell example](https://onelab.info/pipermail/getdp/2015/001774.html).
- [Hugonin and Lalanne, nonlinear complex coordinate transforms, 2005](https://doi.org/10.1364/JOSAA.22.001844).
- [Bermudez et al., PML with unbounded absorbing integral, 2007](https://doi.org/10.1016/j.jcp.2006.09.018).
- The user's local `MultiLayerZY/fem/templates/models-getdp/AcousticScattering/Acoustic2D_penetrable.pro`
  already includes reciprocal-distance and reciprocal-distance-squared absorbing
  profiles. Its phasor convention must not be copied blindly.

The cubic-onset logarithmic profile here is a derived test candidate, not a
claim to reproduce either paper's precise discretization or accuracy theorem.
