# Transmission line parameters

This catalogue describes the formulations implemented by LineCableModels for
the electromagnetic parameters of multiconductor lines and cables. It is a
modelling reference: each entry identifies the physical model, states the
implemented mathematical relation, and cites its literature source.

The catalogue is assembled from the scientific docstrings owned by the
formulation source files. Formula discovery, backend routing, and extension
methods are documented separately in the
[Computational engine](engine.md) developer guide.

Unless an entry states otherwise, ``\omega=2\pi f``, ``j^2=-1``, and the
propagation constant of medium ``m`` is

```math
\gamma_m^2=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m).
```

For conductor pair ``i,j``, ``y_{ij}`` is horizontal separation,
``d_{ij}`` is direct distance, ``D_{ij}`` is image distance, and
``H=h_i+h_j``. For self terms, ``d_{ii}`` is the conductor radius. ``I_n``,
``K_n``, ``Y_n``, and ``\mathbf H_n`` denote modified Bessel, Bessel of the
second kind, and Struve functions as appropriate. All impedances are
per-unit-length quantities. Earth-admittance formulations first produce a
potential-coefficient matrix ``\mathbf P``; the shunt admittance is then
``\mathbf Y=j\omega\mathbf P^{-1}``.

## Internal impedance

These formulations describe longitudinal current diffusion in a conductor.
For a hollow circular conductor, ``a`` and ``b`` are its inner and outer
radii, ``\rho`` is resistivity, ``\mu`` is permeability, and
``m=\sqrt{j\omega\mu/\rho}``.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Engine.InternalImpedance,
    "engine", "internalimpedance", "formulas",
)
```

## Insulation impedance

The magnetic field within a concentric insulation layer contributes a
longitudinal series impedance even when the dielectric is electrically
lossless. For a layer with inner radius ``a``, outer radius ``b``, and relative
permeability ``\mu_r``, the registered formulations use the annular magnetic
field solution.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Engine.InsulationImpedance,
    "engine", "insulationimpedance", "formulas",
)
```

## Insulation admittance

For a homogeneous concentric dielectric, define its complex constitutive
coefficient

```math
\kappa=\sigma+j\omega\varepsilon_0\varepsilon_r.
```

Several dielectric layers are combined radially in series through their
potential coefficients.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Engine.InsulationAdmittance,
    "engine", "insulationadmittance", "formulas",
)
```

## Semiconducting-screen admittance

Semiconducting screens contribute conduction and displacement currents to the
same radial dielectric network as the adjacent insulation layers.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Engine.SemiconAdmittance,
    "engine", "semiconadmittance", "formulas",
)
```

## Earth return impedance

These formulations give the earth-return contribution ``Z_{e,ij}`` to the
series-impedance matrix. A homogeneous model has air ``m=0`` and earth
``m=1``. Stratified models number soil layers downward from one.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Engine.EarthImpedance,
    "engine", "earthimpedance", "formulas",
)
```

## Earth return admittance

These formulations give earth or space potential coefficients. The complete
matrix ``\mathbf P`` is assembled before inversion; equations written for
``P_{e,ij}`` are therefore matrix entries, not elementwise admittances.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Engine.EarthAdmittance,
    "engine", "earthadmittance", "formulas",
)
```

## Modal decomposition

Modal formulations transform fully coupled phase-domain matrices. For current
eigenvectors ``\mathbf T_I`` and voltage eigenvectors ``\mathbf T_V``, the
transformed line parameters are

```math
\mathbf Z_m=\mathbf T_V\mathbf Z\mathbf T_I^{-1},\qquad
\mathbf Y_m=\mathbf T_I\mathbf Y\mathbf T_V^{-1}.
```

The frequency-tracked formulations solve the eigenproblem of
``\mathbf Y\mathbf Z`` and preserve mode identity between adjacent frequency
samples.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.Transforms,
    "transforms", "formulas",
)
```

## Frequency-dependent soil

These constitutive relations map a reference soil resistivity ``\rho_0`` to
frequency-dependent conductivity ``\sigma(f)`` and relative permittivity
``\varepsilon_r(f)``. Magnetic permeability is unchanged. Unless stated
otherwise, conductivity is in S/m and frequency is in hertz.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.EarthProps.FD,
    "earthprops", "fd", "formulas",
)
```

## Equivalent homogeneous earth

These formulations reduce a horizontally layered, nonmagnetic earth to one
frequency-dependent homogeneous material for overhead-line earth-return
calculations. The recursion starts at the bottommost soil layer and advances
upward. It changes constitutive properties only; it does not alter the
physical earth geometry stored in the problem.

```@eval
Main.formulation_catalogue(
    Main.LineCableModels.EarthProps.EHEM,
    "earthprops", "ehem", "formulas",
)
```
