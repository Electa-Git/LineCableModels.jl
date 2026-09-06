# Patel–Triverio cable-hole homogeneous-ground MoM–SO impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Round conductors and circular holes; conductors may be solid or hollow in the stated extension. |
| Calculated quantities | Complete p.u.l. cable series matrix with skin, proximity, a cable hole/tunnel, and homogeneous air–ground return |
| Earth structure | Homogeneous ground half-space below air, with a homogeneous cable hole. |
| Model and approximation | Fourier boundary expansions and MoM quadrature are truncated. The physical air–ground Green function is a spectral integral; no complex-depth fit is substituted. |
| Main source | U. R. Patel and P. Triverio (2015) |
| Citation key(s) | `:Patel2015` |
| Evidence status | Original publication page images checked |

**Description.** Boundary-only cable-hole equivalence coupled directly to the homogeneous air–ground magnetic Green function, yielding one series matrix that includes conductor skin/proximity and earth return.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinal invariance; telegrapher-model end effects omitted. | Stated — §II. |
| Air propagation constant ``γ_air`` | ``k_0=\omega\sqrt{\mu_0\epsilon_0}`` in the Green function. | Stated — below (28). |
| Earth propagation constant ``γ_earth`` | ``k_g=\sqrt{\omega\mu_0(\omega\epsilon_0-j\sigma_g)}``. | Stated — §III-B and (27). |
| Earth permittivity and displacement current | Ground permittivity is set to ``\epsilon_0`` and displacement retained through ``k_g``. | Stated — Fig. 2, (27). |
| Range of validity | Harmonic truncations control numerical convergence; source reports 1 Hz–1 MHz examples, not a universal analytical bound. | Stated — §§III, VII. |
| Earth permeability ``μ_earth`` | Ground uses ``\mu_0``. | Stated — Fig. 2. |
| Arrangement | Arbitrary parallel round solid/hollow conductors in one or multiple circular holes. | Stated — §§II, VI. |
| Earth structure | Homogeneous ground half-space below air, with a homogeneous cable hole. | Stated — Fig. 2. |
| Conductor and insulation geometry | Round conductors and circular holes; conductors may be solid or hollow in the stated extension. | Stated — §§II, VI. |
| Constitutive and field assumptions | Linear isotropic media, 2-D harmonic fields, surface equivalence, Fourier/MoM discretization. | Stated — §§II–V. |
| Conventions | ``e^{j\omega t}`` implied by ``k=\sqrt{\omega\mu(\omega\epsilon-j\sigma)}``; p.u.l. ``R+j\omega L``. | Stated — (5), (31)–(36). |

**Expression.** The air–ground magnetic Green function is

```math
G_g(x,y,x',y')=\frac1{4\pi}\int_{-\infty}^{\infty}
\frac{e^{-j\beta_x(x-x')}}{\sqrt{\beta_x^2-k_g^2}}
\left[e^{-|y-y'|\sqrt{\beta_x^2-k_g^2}}+R_{TM}e^{(y+y')\sqrt{\beta_x^2-k_g^2}}\right]d\beta_x,
\tag{27}
```

```math
R_{TM}=\frac{\sqrt{\beta_x^2-k_g^2}-\sqrt{\beta_x^2-k_0^2}}
{\sqrt{\beta_x^2-k_g^2}+\sqrt{\beta_x^2-k_0^2}}.
\tag{28}
```

After cable-hole discretization,

```math
\widehat{\mathbf A}=-\mu_0(\mathbf1+\mu_0\mathbf G_g\widehat{\mathbf Y}_s)^{-1}\mathbf G_g\mathbf T\mathbf J,
\tag{30}
```

and the result is

```math
\mathbf R=\Re\{[\mathbf U^T(\mathbf1-j\omega\mathbf Y_s\mathbf\Psi)^{-1}\mathbf Y_s\mathbf U]^{-1}\},\quad
\mathbf L=\omega^{-1}\Im\{[\mathbf U^T(\mathbf1-j\omega\mathbf Y_s\mathbf\Psi)^{-1}\mathbf Y_s\mathbf U]^{-1}\}.
\tag{35,36}
```

**Approximation.** Fourier boundary expansions and MoM quadrature are truncated. The physical air–ground Green function is a spectral integral; no complex-depth fit is substituted.

**Limitations.** Flat homogeneous air–ground interface and longitudinally invariant circular holes/conductors. Full evaluation requires the hole operator/transformation matrices in (18)–(25) and ``\mathbf\Psi`` in (34).

**Reference.** [Patel2015](@cite), equations (1)–(36), especially (27)–(36), printed pp. 2113–2114.

**Transcription source.** Original IEEE page images. Spectral prefactor, exponents, reflection coefficient and final nested matrix inversions were visually verified.

## Source transcription

Equation (34) defines ``\mathbf\Psi`` by combining the cable-hole transformation, homogeneous-ground Green matrix and particular/general hole solutions. That dependency is retained by exact locator and is not simplified into a scalar earth-return correction.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``G_g,\mathbf G_g`` | unchanged | continuous/discretized air–ground magnetic Green function | source normalized |
| ``\widehat{\mathbf Y}_s,\mathbf T`` | unchanged | empty-hole surface operator and conductor-to-hole map | boundary operators |
| ``\mathbf\Psi`` | unchanged | complete exterior potential operator | source normalized |
| ``\mathbf R,\mathbf L`` | unchanged | p.u.l. resistance/inductance matrices | ``\Omega/\mathrm m``, ``\mathrm H/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

The direct Green-function coupling is the distinction from the 2014 additive ground correction. The paper labels the cable-hole surface representation as novel.

## Limitations and discrepancies

- Equation (27) uses a Fourier sign/exponential convention tied to the paper's coordinate placement; it must not be mixed with a different spectral convention without derivation.
- Insulation electric losses and shunt admittance are outside this series-only result.

