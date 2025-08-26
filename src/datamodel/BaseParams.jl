"""
	LineCableModels.DataModel.BaseParams

The [`BaseParams`](@ref) submodule provides fundamental functions for determining the base electrical parameters (R, L, C, G) of cable components within the [`LineCableModels.DataModel`](@ref) module. This includes implementations of standard engineering formulas for resistance, inductance, and geometric parameters of various conductor configurations.

# Overview

- Implements basic electrical engineering formulas for calculating DC resistance and inductance of different conductor geometries (tubular, strip, wire arrays).
- Implements basic formulas for capacitance and dielectric losses in insulators and semiconductors.
- Provides functions for temperature correction of material properties.
- Calculates geometric mean radii for different conductor configurations.
- Includes functions for determining the effective length for helical wire arrangements.
- Calculates equivalent electrical parameters and correction factors for different geometries and configurations.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module BaseParams

# Export public API
export calc_equivalent_alpha
export calc_parallel_equivalent
export calc_helical_params
export calc_strip_resistance
export calc_temperature_correction
export calc_tubular_resistance
export calc_tubular_inductance
export calc_wirearray_coords
export calc_inductance_trifoil
export calc_wirearray_gmr
export calc_tubular_gmr
export calc_equivalent_mu
export calc_shunt_capacitance
export calc_shunt_conductance
export calc_equivalent_gmr
export calc_gmd
export calc_solenoid_correction
export calc_equivalent_rho
export calc_equivalent_eps
export calc_equivalent_lossfact
export calc_sigma_lossfact

# Load common dependencies
include("../utils/commondeps.jl")
using ...LineCableModels

# Module-specific dependencies
using Measurements
using ...Utils
using ..DataModel: ConductorGroup, WireArray
import ..DataModel: AbstractCablePart

"""
$(TYPEDSIGNATURES)

Calculates the equivalent temperature coefficient of resistance (`alpha`) when two conductors are connected in parallel, by cross-weighted-resistance averaging:

```math
\\alpha_{eq} = \\frac{\\alpha_1 R_2 + \\alpha_2 R1}{R_1 + R_2}
```
where ``\\alpha_1``, ``\\alpha_2`` are the temperature coefficients of the conductors, and ``R_1``, ``R_2`` are the respective resistances.

# Arguments

- `alpha1`: Temperature coefficient of resistance of the first conductor \\[1/°C\\].
- `R1`: Resistance of the first conductor \\[Ω\\].
- `alpha2`: Temperature coefficient of resistance of the second conductor \\[1/°C\\].
- `R2`: Resistance of the second conductor \\[Ω\\].

# Returns

- The equivalent temperature coefficient \\[1/°C\\] for the parallel combination.

# Examples

```julia
alpha_conductor = 0.00393  # Copper
alpha_new_part = 0.00403   # Aluminum
R_conductor = 0.5
R_new_part = 1.0
alpha_eq = $(FUNCTIONNAME)(alpha_conductor, R_conductor, alpha_new_part, R_new_part)
println(alpha_eq)  # Output: 0.00396 (approximately)
```
"""
function calc_equivalent_alpha(alpha1::T, R1::T, alpha2::T, R2::T) where {T<:REALSCALAR}
    return (alpha1 * R2 + alpha2 * R1) / (R1 + R2)
end

function calc_equivalent_alpha(alpha1, R1, alpha2, R2)
    T = resolve_T(alpha1, R1, alpha2, R2)
    return calc_equivalent_alpha(
        coerce_to_T(alpha1, T),
        coerce_to_T(R1, T),
        coerce_to_T(alpha2, T),
        coerce_to_T(R2, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the parallel equivalent of two impedances (or series equivalent of two admittances):

```math
Z_{eq} = \\frac{Z_1 Z_2}{Z_1 + Z_2}
```

This expression, when applied recursively to [`LineCableModels.DataModel.WireArray`](@ref) objects, implements the formula for the hexagonal wiring pattern described in CIGRE TB-345 [app14198982](@cite) [cigre345](@cite):

```math
\\frac{1}{R_{\\text{dc}}} = \\frac{\\pi d^2}{4 \\rho} \\left( 1 + \\sum_{1}^{n} \\frac{6n}{k_n} \\right)
```

```math
k_n = \\left[ 1 + \\left( \\pi \\frac{D_n}{\\lambda_n} \\right)^2 \\right]^{1/2}
```

where ``R_{\\text{dc}}`` is the DC resistance, ``d`` is the diameter of each wire, ``\rho`` is the resistivity, ``n`` is the number of layers following the hexagonal pattern, ``D_n`` is the diameter of the ``n``-th layer, and ``\\lambda_n `` is the pitch length of the ``n``-th layer, obtained using [`calc_helical_params`](@ref).

# Arguments

- `Z1`: The total impedance of the existing system \\[Ω\\].
- `Z2`: The impedance of the new layer being added \\[Ω\\].

# Returns

- The parallel equivalent impedance \\[Ω\\].

# Examples

```julia
Z1 = 5.0
Z2 = 10.0
Req = $(FUNCTIONNAME)(Z1, Z2)
println(Req) # Outputs: 3.3333333333333335
```

# See also

- [`calc_helical_params`](@ref)
"""
function calc_parallel_equivalent(Z1::T, Z2::T) where {T<:Union{REALSCALAR,COMPLEXSCALAR}}

    # Case 1: Inf / Inf -> NaN
    # The parallel combination of an open circuit (Inf) and any finite impedance is the finite impedance.
    if isinf(Z1)
        return Z2
    elseif isinf(Z2)
        return Z1
    end

    # Case 2: 0 / 0 -> NaN
    # The parallel combination of two short circuits (0) is a short circuit.
    # The standard formula works fine if only one is zero, but not if both are.
    if iszero(Z1) && iszero(Z2)
        return zero(T)
    end
    return (Z1 * Z2) / (Z1 + Z2)
end

function calc_parallel_equivalent(Z1, Z2)
    T = resolve_T(Z1, Z2)
    return calc_parallel_equivalent(
        coerce_to_T(Z1, T),
        coerce_to_T(Z2, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the mean diameter, pitch length, and overlength based on cable geometry parameters. The lay ratio is defined as the ratio of the pitch length ``L_p`` to the external diameter ``D_e``:

```math
\\lambda = \\frac{L_p}{D_e}
```
where ``D_e`` and ``L_p`` are the dimensions represented in the figure.

![](./assets/lay_ratio.svg)

# Arguments

- `radius_in`: Inner radius of the cable layer \\[m\\].
- `radius_ext`: Outer radius of the cable layer \\[m\\].
- `lay_ratio`: Ratio of the pitch (lay) length to the external diameter of the corresponding layer of wires \\[dimensionless\\].

# Returns

- `mean_diameter`: Mean diameter of the cable layer \\[m\\].
- `pitch_length`: The length over which the strands complete one full twist \\[m\\].
- `overlength`: Effective length increase resulting from the helical path \\[1/m\\].

# Notes

Reference values for `lay_ratio` are given under standard EN 50182 [CENELEC50182](@cite):

| Conductor type | Steel wires | Aluminum wires | Lay ratio - Steel | Lay ratio - Aluminum |
|---------------|----------------------|---------------------|----------------------|-------------------|
| AAAC 4 layers | - | 61 (1/6/12/18/24) | - | 15/13.5/12.5/11 |
| ACSR 3 layers | 7 (1/6) | 54 (12/18/24) | 19 | 15/13/11.5 |
| ACSR 2 layers | 7 (1/6) | 26 (10/16) | 19 | 14/11.5 |
| ACSR 1 layer | 7 (1/6) | 10 | 19 | 14 |
| ACCC/TW | - | 36 (8/12/16) | - | 15/13.5/11.5 |	

# Examples

```julia
radius_in = 0.01
radius_ext = 0.015
lay_ratio = 12

mean_diam, pitch, overlength = $(FUNCTIONNAME)(radius_in, radius_ext, lay_ratio)
# mean_diam ≈ 0.025 [m]
# pitch ≈ 0.3 [m]
# overlength > 1.0 [1/m]
```
"""
function calc_helical_params(radius_in::T, radius_ext::T, lay_ratio::T) where {T<:REALSCALAR}
    mean_diameter = 2 * (radius_in + (radius_ext - radius_in) / 2)
    pitch_length = lay_ratio * mean_diameter
    overlength = !isapprox(pitch_length, 0.0) ? sqrt(1 + (π * mean_diameter / pitch_length)^2) : 1

    return mean_diameter, pitch_length, overlength
end

function calc_helical_params(radius_in, radius_ext, lay_ratio)
    T = resolve_T(radius_in, radius_ext, lay_ratio)
    return calc_helical_params(
        coerce_to_T(radius_in, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(lay_ratio, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the DC resistance of a strip conductor based on its geometric and material properties, using the basic resistance formula in terms of the resistivity and cross-sectional area:

```math
R = \\rho \\frac{\\ell}{W T}
```
where ``\\ell`` is the length of the strip, ``W`` is the width, and ``T`` is the thickness. The length is assumed to be infinite in the direction of current flow, so the resistance is calculated per unit length.

# Arguments

- `thickness`: Thickness of the strip \\[m\\].
- `width`: Width of the strip \\[m\\].
- `rho`: Electrical resistivity of the conductor material \\[Ω·m\\].
- `alpha`: Temperature coefficient of resistivity \\[1/°C\\].
- `T0`: Reference temperature for the material properties \\[°C\\].
- `Top`: Operating temperature of the conductor \\[°C\\].

# Returns

- DC resistance of the strip conductor \\[Ω\\].

# Examples

```julia
thickness = 0.002
width = 0.05
rho = 1.7241e-8
alpha = 0.00393
T0 = 20
T = 25
resistance = $(FUNCTIONNAME)(thickness, width, rho, alpha, T0, T)
# Output: ~0.0001758 Ω
```

# See also  

- [`calc_temperature_correction`](@ref)
"""
function calc_strip_resistance(
    thickness::T,
    width::T,
    rho::T,
    alpha::T,
    T0::T,
    Top::T,
) where {T<:REALSCALAR}

    cross_section = thickness * width
    return calc_temperature_correction(alpha, Top, T0) * rho / cross_section
end

function calc_strip_resistance(thickness, width, rho, alpha, T0, Top)
    T = resolve_T(thickness, width, rho, alpha, T0, Top)
    return calc_strip_resistance(
        coerce_to_T(thickness, T),
        coerce_to_T(width, T),
        coerce_to_T(rho, T),
        coerce_to_T(alpha, T),
        coerce_to_T(T0, T),
        coerce_to_T(Top, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the temperature correction factor for material properties based on the standard linear temperature model [cigre345](@cite):

```math
k(T) = 1 + \\alpha (T - T_0)
```
where ``\\alpha`` is the temperature coefficient of the material resistivity, ``T`` is the operating temperature, and ``T_0`` is the reference temperature. 

# Arguments

- `alpha`: Temperature coefficient of the material property \\[1/°C\\].
- `T`: Current temperature \\[°C\\].
- `T0`: Reference temperature at which the base material property was measured \\[°C\\]. Defaults to T₀.

# Returns

- Temperature correction factor to be applied to the material property \\[dimensionless\\].

# Examples

```julia
    # Copper resistivity correction (alpha = 0.00393 [1/°C])
    k = $(FUNCTIONNAME)(0.00393, 75.0, 20.0)  # Expected output: 1.2161
```
"""
function calc_temperature_correction(alpha::T, Top::T, T0::T=T₀) where {T<:REALSCALAR}
    @assert abs(Top - T0) < ΔTmax """
    Temperature is outside the valid range for linear resistivity model:
    Top = $Top
    T0 = $T0
    ΔTmax = $ΔTmax
    |Top - T0| = $(abs(Top - T0))"""
    return 1 + alpha * (Top - T0)
end

function calc_temperature_correction(alpha, Top, T0=T₀)
    T = resolve_T(alpha, Top, T0)
    return calc_temperature_correction(
        coerce_to_T(alpha, T),
        coerce_to_T(Top, T),
        coerce_to_T(T0, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the DC resistance of a tubular conductor based on its geometric and material properties, using the resistivity and cross-sectional area of a hollow cylinder with radii ``r_{in}`` and ``r_{ext}``:

```math
R = \\rho \\frac{\\ell}{\\pi (r_{ext}^2 - r_{in}^2)}
```
where ``\\ell`` is the length of the conductor, ``r_{in}`` and ``r_{ext}`` are the inner and outer radii, respectively. The length is assumed to be infinite in the direction of current flow, so the resistance is calculated per unit length.

# Arguments

- `radius_in`: Internal radius of the tubular conductor \\[m\\].
- `radius_ext`: External radius of the tubular conductor \\[m\\].
- `rho`: Electrical resistivity of the conductor material \\[Ω·m\\].
- `alpha`: Temperature coefficient of resistivity \\[1/°C\\].
- `T0`: Reference temperature for the material properties \\[°C\\].
- `Top`: Operating temperature of the conductor \\[°C\\].

# Returns

- DC resistance of the tubular conductor \\[Ω\\].

# Examples

```julia
radius_in = 0.01
radius_ext = 0.02
rho = 1.7241e-8
alpha = 0.00393
T0 = 20
T = 25
resistance = $(FUNCTIONNAME)(radius_in, radius_ext, rho, alpha, T0, T)
# Output: ~9.10e-8 Ω
```

# See also

- [`calc_temperature_correction`](@ref)
"""
function calc_tubular_resistance(
    radius_in::T,
    radius_ext::T,
    rho::T,
    alpha::T,
    T0::T,
    Top::T,
) where {T<:REALSCALAR}
    cross_section = π * (radius_ext^2 - radius_in^2)
    return calc_temperature_correction(alpha, Top, T0) * rho / cross_section
end

function calc_tubular_resistance(radius_in, radius_ext, rho, alpha, T0, Top)
    T = resolve_T(radius_in, radius_ext, rho, alpha, T0, Top)
    return calc_tubular_resistance(
        coerce_to_T(radius_in, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(rho, T),
        coerce_to_T(alpha, T),
        coerce_to_T(T0, T),
        coerce_to_T(Top, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the inductance of a tubular conductor per unit length, disregarding skin-effects (DC approximation) [916943](@cite) [cigre345](@cite) [1458878](@cite):

```math
L = \\frac{\\mu_r \\mu_0}{2 \\pi} \\log \\left( \\frac{r_{ext}}{r_{in}} \\right)
```
where ``\\mu_r`` is the relative permeability of the conductor material, ``\\mu_0`` is the vacuum permeability, and ``r_{in}`` and ``r_{ext}`` are the inner and outer radii of the conductor, respectively.

# Arguments

- `radius_in`: Internal radius of the tubular conductor \\[m\\].
- `radius_ext`: External radius of the tubular conductor \\[m\\].
- `mu_r`: Relative permeability of the conductor material \\[dimensionless\\].

# Returns

- Internal inductance of the tubular conductor per unit length \\[H/m\\].

# Examples

```julia
radius_in = 0.01
radius_ext = 0.02
mu_r = 1.0
L = $(FUNCTIONNAME)(radius_in, radius_ext, mu_r)
# Output: ~2.31e-7 H/m
```

# See also

- [`calc_tubular_resistance`](@ref)
"""
function calc_tubular_inductance(radius_in::T, radius_ext::T, mu_r::T) where {T<:REALSCALAR}
    return mu_r * μ₀ / (2 * π) * log(radius_ext / radius_in)
end

function calc_tubular_inductance(radius_in, radius_ext, mu_r)
    T = resolve_T(radius_in, radius_ext, mu_r)
    return calc_tubular_inductance(
        coerce_to_T(radius_in, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(mu_r, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the center coordinates of wires arranged in a circular pattern.

# Arguments

- `num_wires`: Number of wires in the circular arrangement \\[dimensionless\\].
- `radius_wire`: Radius of each individual wire \\[m\\].
- `radius_in`: Inner radius of the wire array (to wire centers) \\[m\\].
- `C`: Optional tuple representing the center coordinates of the circular arrangement \\[m\\]. Default is (0.0, 0.0).

# Returns

- Vector of tuples, where each tuple contains the `(x, y)` coordinates \\[m\\] of the center of a wire.

# Examples

```julia
# Create a 7-wire array with 2mm wire radius and 1cm inner radius
wire_coords = $(FUNCTIONNAME)(7, 0.002, 0.01)
println(wire_coords[1]) # Output: First wire coordinates

# Create a wire array with custom center position
wire_coords = $(FUNCTIONNAME)(7, 0.002, 0.01, C=(0.5, 0.3))
```

# See also

- [`LineCableModels.DataModel.WireArray`](@ref)
"""
function calc_wirearray_coords(
    num_wires::U,
    radius_wire::T,
    radius_in::T,
    C::Tuple{T,T},
) where {T<:REALSCALAR,U<:Int}
    wire_coords = Tuple{T,T}[]  # Global coordinates of all wires
    lay_radius = num_wires == 1 ? 0 : radius_in + radius_wire

    # Calculate the angle between each wire
    angle_step = 2 * π / num_wires
    for i in 0:num_wires-1
        angle = i * angle_step
        x = C[1] + lay_radius * cos(angle)
        y = C[2] + lay_radius * sin(angle)
        push!(wire_coords, (x, y))  # Add wire center
    end
    return wire_coords
end

function calc_wirearray_coords(num_wires::Int, radius_wire, radius_in; C=nothing)
    T = C === nothing ? resolve_T(radius_wire, radius_in) : resolve_T(radius_wire, radius_in, C...)
    C_val = C === nothing ? coerce_to_T((0.0, 0.0), T) : coerce_to_T(C, T)
    return calc_wirearray_coords(
        num_wires,
        coerce_to_T(radius_wire, T),
        coerce_to_T(radius_in, T),
        C_val,
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the positive-sequence inductance of a trifoil-configured cable system composed of core/screen assuming solid bonding, using the formula given under section 4.2.4.3 of CIGRE TB-531:

```math
Z_d = \\left[Z_a - Z_x\\right] - \\frac{\\left( Z_m - Z_x \\right)^2}{Z_s - Z_x}
```
```math
L = \\mathfrak{Im}\\left(\\frac{Z_d}{\\omega}\\right)
```
where ``Z_a``, ``Z_s`` are the self impedances of the core conductor and the screen, and ``Z_m``, and ``Z_x`` are the mutual impedances between core/screen and between cables, respectively, as per sections 4.2.3.4, 4.2.3.5, 4.2.3.6 and 4.2.3.8 of the same document [cigre531](@cite).

# Arguments

- `r_in_co`: Internal radius of the phase conductor \\[m\\].
- `r_ext_co`: External radius of the phase conductor \\[m\\].
- `rho_co`: Electrical resistivity of the phase conductor material \\[Ω·m\\].
- `mu_r_co`: Relative permeability of the phase conductor material \\[dimensionless\\].
- `r_in_scr`: Internal radius of the metallic screen \\[m\\].
- `r_ext_scr`: External radius of the metallic screen \\[m\\].
- `rho_scr`: Electrical resistivity of the metallic screen material \\[Ω·m\\].
- `mu_r_scr`: Relative permeability of the screen conductor material \\[dimensionless\\].
- `S`: Spacing between conductors in trifoil configuration \\[m\\].
- `rho_e`: Soil resistivity \\[Ω·m\\]. Default: 100 Ω·m.
- `f`: Frequency \\[Hz\\]. Default: [`f₀`](@ref).

# Returns

- Positive-sequence inductance per unit length of the cable system \\[H/m\\].

# Examples

```julia
L = $(FUNCTIONNAME)(0.01, 0.015, 1.72e-8, 1.0, 0.02, 0.025, 2.83e-8, 1.0, S=0.1, rho_e=50, f=50)
println(L) # Output: Inductance value in H/m
```

# See also

- [`calc_tubular_gmr`](@ref)
"""
function calc_inductance_trifoil(
    r_in_co::T,
    r_ext_co::T,
    rho_co::T,
    mu_r_co::T,
    r_in_scr::T,
    r_ext_scr::T,
    rho_scr::T,
    mu_r_scr::T,
    S::T;
    rho_e::T=100.0,
    f::T=f₀,
) where {T<:REALSCALAR}

    ω = 2 * π * f
    C = μ₀ / (2π)

    # Compute simplified earth return depth
    DE = 659.0 * sqrt(rho_e / f)

    # Compute R'_E
    RpE = (ω * μ₀) / 8.0

    # Compute Xa
    GMRa = calc_tubular_gmr(r_ext_co, r_in_co, mu_r_co)
    Xa = (ω * C) * log(DE / GMRa)

    # Self impedance of a phase conductor with earth return
    Ra = rho_co / (π * (r_ext_co^2 - r_in_co^2))
    Za = RpE + Ra + im * Xa

    # Compute rs
    GMRscr = calc_tubular_gmr(r_ext_scr, r_in_scr, mu_r_scr)
    # Compute Xs
    Xs = (ω * C) * log(DE / GMRscr)

    # Self impedance of metal screen with earth return
    Rs = rho_scr / (π * (r_ext_scr^2 - r_in_scr^2))
    Zs = RpE + Rs + im * Xs

    # Mutual impedance between phase conductor and screen
    Zm = RpE + im * Xs

    # Compute GMD
    GMD = S # trifoil, for flat use: 2^(1/3) * S

    # Compute Xap
    Xap = (ω * C) * log(DE / GMD)

    # Equivalent mutual impedances between cables
    Zx = RpE + im * Xap

    # Formula from CIGRE TB-531, 4.2.4.3, solid bonding
    Z1_sb = (Za - Zx) - ((Zm - Zx)^2 / (Zs - Zx))

    # Likewise, but for single point bonding
    # Z1_sp = (Za - Zx)
    return imag(Z1_sb) / ω
end

function calc_inductance_trifoil(
    r_in_co,
    r_ext_co,
    rho_co,
    mu_r_co,
    r_in_scr,
    r_ext_scr,
    rho_scr,
    mu_r_scr,
    S;
    rho_e=100.0,
    f=f₀,
)
    T = resolve_T(r_in_co, r_ext_co, rho_co, mu_r_co, r_in_scr, r_ext_scr, rho_scr, mu_r_scr, S, rho_e, f)
    return calc_inductance_trifoil(
        coerce_to_T(r_in_co, T),
        coerce_to_T(r_ext_co, T),
        coerce_to_T(rho_co, T),
        coerce_to_T(mu_r_co, T),
        coerce_to_T(r_in_scr, T),
        coerce_to_T(r_ext_scr, T),
        coerce_to_T(rho_scr, T),
        coerce_to_T(mu_r_scr, T),
        coerce_to_T(S, T);
        rho_e=coerce_to_T(rho_e, T),
        f=coerce_to_T(f, T)
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the geometric mean radius (GMR) of a circular wire array, using formula (62), page 335, of the book by Edward Rosa [rosa1908](@cite):

```math
GMR = \\sqrt[a] {r n a^{n-1}}
```

where ``a`` is the layout radius, ``n`` is the number of wires, and ``r`` is the radius of each wire.

# Arguments

- `lay_rad`: Layout radius of the wire array \\[m\\].
- `N`: Number of wires in the array \\[dimensionless\\].
- `rad_wire`: Radius of an individual wire \\[m\\].
- `mu_r`: Relative permeability of the wire material \\[dimensionless\\].

# Returns

- Geometric mean radius (GMR) of the wire array \\[m\\].

# Examples

```julia
lay_rad = 0.05
N = 7
rad_wire = 0.002
mu_r = 1.0
gmr = $(FUNCTIONNAME)(lay_rad, N, rad_wire, mu_r)
println(gmr) # Expected output: 0.01187... [m]
```
"""
function calc_wirearray_gmr(lay_rad::T, N::Int, rad_wire::T, mu_r::T) where {T<:REALSCALAR}
    gmr_wire = rad_wire * exp(-mu_r / 4)
    log_gmr_array = log(gmr_wire * N * lay_rad^(N - 1)) / N
    return exp(log_gmr_array)
end

function calc_wirearray_gmr(lay_rad, N::Int, rad_wire, mu_r)
    T = resolve_T(lay_rad, rad_wire, mu_r)
    return calc_wirearray_gmr(
        coerce_to_T(lay_rad, T),
        N,
        coerce_to_T(rad_wire, T),
        coerce_to_T(mu_r, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the geometric mean radius (GMR) of a tubular conductor, using [6521501](@cite):

```math
\\log GMR = \\log r_2 - \\mu_r \\left[ \\frac{r_1^4}{\\left(r_2^2 - r_1^2\\right)^2} \\log\\left(\\frac{r_2}{r_1}\\right) - \\frac{3r_1^2 - r_2^2}{4\\left(r_2^2 - r_1^2\\right)} \\right]
```

where ``\\mu_r`` is the material magnetic permeability (relative to free space), ``r_1`` and ``r_2`` are the inner and outer radii of the tubular conductor, respectively. If ``r_2`` is approximately equal to ``r_1`` , the tube collapses into a thin shell, and the GMR is equal to ``r_2``. If the tube becomes infinitely thick (e.g., ``r_2 \\gg r_1``), the GMR diverges to infinity.

# Arguments

- `radius_ext`: External radius of the tubular conductor \\[m\\].
- `radius_in`: Internal radius of the tubular conductor \\[m\\].
- `mu_r`: Relative permeability of the conductor material \\[dimensionless\\].

# Returns

- Geometric mean radius (GMR) of the tubular conductor \\[m\\].

# Errors

- Throws `ArgumentError` if `radius_ext` is less than `radius_in`.

# Examples

```julia
radius_ext = 0.02
radius_in = 0.01
mu_r = 1.0
gmr = $(FUNCTIONNAME)(radius_ext, radius_in, mu_r)
println(gmr) # Expected output: ~0.0135 [m]
```
"""
function calc_tubular_gmr(radius_ext::T, radius_in::T, mu_r::T) where {T<:REALSCALAR}
    if (radius_ext < radius_in) || (radius_ext <= 0.0)
        throw(ArgumentError("Invalid parameters: radius_ext must be >= radius_in and positive."))
    end

    # Constants
    if isapprox(radius_in, radius_ext)
        # Tube collapses into a thin shell with infinitesimal thickness and the GMR is simply the radius
        gmr = radius_ext
    elseif abs(radius_in / radius_ext) < eps() && abs(radius_in) > TOL
        # Tube becomes infinitely thick up to floating point precision
        gmr = Inf
    else
        is_solid = isapprox(radius_in, 0.0)
        term1 = is_solid ? 0.0 :
                (radius_in^4 / (radius_ext^2 - radius_in^2)^2) * log(radius_ext / radius_in)
        term2 = (3 * radius_in^2 - radius_ext^2) / (4 * (radius_ext^2 - radius_in^2))
        Lin = (μ₀ * mu_r / (2 * π)) * (term1 - term2)

        # Compute the GMR
        gmr = exp(log(radius_ext) - (2 * π / μ₀) * Lin)
    end

    return gmr
end

function calc_tubular_gmr(radius_ext, radius_in, mu_r)
    T = resolve_T(radius_ext, radius_in, mu_r)
    return calc_tubular_gmr(
        coerce_to_T(radius_ext, T),
        coerce_to_T(radius_in, T),
        coerce_to_T(mu_r, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the relative permeability (`mu_r`) based on the geometric mean radius (GMR) and conductor dimensions, by executing the inverse of [`calc_tubular_gmr`](@ref), and solving for `mu_r`:

```math
\\log GMR = \\log r_2 - \\mu_r \\left[ \\frac{r_1^4}{\\left(r_2^2 - r_1^2\\right)^2} \\log\\left(\\frac{r_2}{r_1}\\right) - \\frac{3r_1^2 - r_2^2}{4\\left(r_2^2 - r_1^2\\right)} \\right]
```

```math
\\mu_r = -\\frac{\\left(\\log GMR - \\log r_2\\right)}{\\frac{r_1^4}{\\left(r_2^2 - r_1^2\\right)^2} \\log\\left(\\frac{r_2}{r_1}\\right) - \\frac{3r_1^2 - r_2^2}{4\\left(r_2^2 - r_1^2\\right)}}
```

where ``r_1`` is the inner radius and ``r_2`` is the outer radius.

# Arguments

- `gmr`: Geometric mean radius of the conductor \\[m\\].
- `radius_ext`: External radius of the conductor \\[m\\].
- `radius_in`: Internal radius of the conductor \\[m\\].

# Returns

- Relative permeability (`mu_r`) of the conductor material \\[dimensionless\\].

# Errors

- Throws `ArgumentError` if `radius_ext` is less than `radius_in`.

# Notes

Assumes a tubular geometry for the conductor, reducing to the solid case if `radius_in` is zero.

# Examples

```julia
gmr = 0.015
radius_ext = 0.02
radius_in = 0.01
mu_r = $(FUNCTIONNAME)(gmr, radius_ext, radius_in)
println(mu_r) # Expected output: ~1.7 [dimensionless]
```

# See also
- [`calc_tubular_gmr`](@ref)
"""
function calc_equivalent_mu(gmr::T, radius_ext::T, radius_in::T) where {T<:REALSCALAR}
    if (radius_ext < radius_in) || (radius_ext <= 0.0)
        throw(ArgumentError("Invalid parameters: radius_ext must be >= radius_in and positive."))
    end
    is_solid = isapprox(radius_in, 0.0) || isapprox(radius_in, radius_ext)
    term1 = is_solid ? 0.0 :
            (radius_in^4 / (radius_ext^2 - radius_in^2)^2) * log(radius_ext / radius_in)
    term2 = (3 * radius_in^2 - radius_ext^2) / (4 * (radius_ext^2 - radius_in^2))
    # Compute the log difference
    log_diff = log(gmr) - log(radius_ext)

    # Compute mu_r
    mu_r = -log_diff / (term1 - term2)

    return mu_r
end

function calc_equivalent_mu(gmr, radius_ext, radius_in)
    T = resolve_T(gmr, radius_ext, radius_in)
    return calc_equivalent_mu(
        coerce_to_T(gmr, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(radius_in, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the shunt capacitance per unit length of a coaxial structure, using the standard formula for the capacitance of a coaxial structure [cigre531](@cite) [916943](@cite) [1458878](@cite):

```math
C = \\frac{2 \\pi \\varepsilon_0 \\varepsilon_r}{\\log \\left(\\frac{r_{ext}}{r_{in}}\\right)}
```
where ``\\varepsilon_0`` is the vacuum permittivity, ``\\varepsilon_r`` is the relative permittivity of the dielectric material, and ``r_{in}`` and ``r_{ext}`` are the inner and outer radii of the coaxial structure, respectively.

# Arguments

- `radius_in`: Internal radius of the coaxial structure \\[m\\].
- `radius_ext`: External radius of the coaxial structure \\[m\\].
- `epsr`: Relative permittivity of the dielectric material \\[dimensionless\\].

# Returns

- Shunt capacitance per unit length \\[F/m\\].

# Examples

```julia
radius_in = 0.01
radius_ext = 0.02
epsr = 2.3
capacitance = $(FUNCTIONNAME)(radius_in, radius_ext, epsr)
println(capacitance) # Expected output: ~1.24e-10 [F/m]
```
"""
function calc_shunt_capacitance(radius_in::T, radius_ext::T, epsr::T) where {T<:REALSCALAR}
    return 2 * π * ε₀ * epsr / log(radius_ext / radius_in)
end

function calc_shunt_capacitance(radius_in, radius_ext, epsr)
    T = resolve_T(radius_in, radius_ext, epsr)
    return calc_shunt_capacitance(
        coerce_to_T(radius_in, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(epsr, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the shunt conductance per unit length of a coaxial structure, using the improved model reported in [916943](@cite) [Karmokar2025](@cite) [4389974](@cite):

```math
G = \\frac{2\\pi\\sigma}{\\log(\\frac{r_{ext}}{r_{in}})}
```
where ``\\sigma = \\frac{1}{\\rho}`` is the conductivity of the dielectric/semiconducting material, ``r_{in}`` is the internal radius, and ``r_{ext}`` is the external radius of the coaxial structure.

# Arguments

- `radius_in`: Internal radius of the coaxial structure \\[m\\].
- `radius_ext`: External radius of the coaxial structure \\[m\\].
- `rho`: Resistivity of the dielectric/semiconducting material \\[Ω·m\\].

# Returns

- Shunt conductance per unit length \\[S·m\\].

# Examples

```julia
radius_in = 0.01
radius_ext = 0.02
rho = 1e9
g = $(FUNCTIONNAME)(radius_in, radius_ext, rho)
println(g) # Expected output: 2.7169e-9 [S·m]
```
"""
function calc_shunt_conductance(radius_in::T, radius_ext::T, rho::T) where {T<:REALSCALAR}
    return 2 * π * (1 / rho) / log(radius_ext / radius_in)
end

function calc_shunt_conductance(radius_in, radius_ext, rho)
    T = resolve_T(radius_in, radius_ext, rho)
    return calc_shunt_conductance(
        coerce_to_T(radius_in, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(rho, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the equivalent geometric mean radius (GMR) of a conductor after adding a new layer, by recursive application of the multizone stranded conductor defined as [yang2008gmr](@cite):

```math
GMR_{eq} = {GMR_{i-1}}^{\\beta^2} \\cdot {GMR_{i}}^{(1-\\beta)^2} \\cdot {GMD}^{2\\beta(1-\\beta)}
```
```math
\\beta = \\frac{S_{i-1}}{S_{i-1} + S_{i}}
```
where:
- ``S_{i-1}`` is the cumulative cross-sectional area of the existing cable part, ``S_{i}`` is the total cross-sectional area after inclusion of the conducting layer ``{i}``.
- ``GMR_{i-1}`` is the cumulative GMR of the existing cable part, ``GMR_{i}`` is the GMR of the conducting layer ``{i}``.
- ``GMD`` is the geometric mean distance between the existing cable part and the new layer, calculated using [`calc_gmd`](@ref).

# Arguments

- `existing`: The existing cable part ([`AbstractCablePart`](@ref)).
- `new_layer`: The new layer being added ([`AbstractCablePart`](@ref)).

# Returns

- Updated equivalent GMR of the combined conductor \\[m\\].

# Examples

```julia
material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
conductor = Conductor(Strip(0.01, 0.002, 0.05, 10, material_props))
new_layer = WireArray(0.02, 0.002, 7, 15, material_props)
equivalent_gmr = $(FUNCTIONNAME)(conductor, new_layer)  # Expected output: Updated GMR value [m]
```

# See also

- [`calc_gmd`](@ref)
"""
function calc_equivalent_gmr(existing::T, new_layer::U) where {T<:AbstractCablePart,U<:AbstractCablePart}
    beta = existing.cross_section / (existing.cross_section + new_layer.cross_section)
    current_conductor = existing isa ConductorGroup ? existing.layers[end] : existing
    gmd = calc_gmd(current_conductor, new_layer)
    return existing.gmr^(beta^2) * new_layer.gmr^((1 - beta)^2) *
           gmd^(2 * beta * (1 - beta))
end

"""
$(TYPEDSIGNATURES)

Calculates the geometric mean distance (GMD) between two cable parts, by using the  definition described in Grover [grover1981inductance](@cite):

```math
\\log GMD = \\left(\\frac{\\sum_{i=1}^{n_1}\\sum_{j=1}^{n_2} (s_1 \\cdot s_2) \\cdot \\log(d_{ij})}{\\sum_{i=1}^{n_1}\\sum_{j=1}^{n_2} (s_1 \\cdot s_2)}\\right)
```

where:
- ``d_{ij}`` is the Euclidean distance between elements ``i`` and ``j``.
- ``s_1`` and ``s_2`` are the cross-sectional areas of the respective elements.
- ``n_1`` and ``n_2`` are the number of sub-elements in each cable part.

# Arguments

- `co1`: First cable part ([`AbstractCablePart`](@ref)).
- `co2`: Second cable part ([`AbstractCablePart`](@ref)).

# Returns

- Geometric mean distance between the cable parts \\[m\\].

# Notes

For concentric structures, the GMD converges to the external radii of the outermost element.

!!! info "Numerical stability"
	This implementation uses a weighted sum of logarithms rather than the traditional product formula ``\\Pi(d_{ij})^{(1/n)}`` found in textbooks. The logarithmic approach prevents numerical underflow/overflow when dealing with many conductors or extreme distance ratios, making it significantly more stable for practical calculations.

# Examples

```julia
material_props = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
wire_array1 = WireArray(0.01, 0.002, 7, 10, material_props)
wire_array2 = WireArray(0.02, 0.002, 7, 15, material_props)
gmd = $(FUNCTIONNAME)(wire_array1, wire_array2)  # Expected output: GMD value [m]

strip = Strip(0.01, 0.002, 0.05, 10, material_props)
tubular = Tubular(0.01, 0.02, material_props)
gmd = $(FUNCTIONNAME)(strip, tubular)  # Expected output: GMD value [m]
```

# See also

- [`calc_wirearray_coords`](@ref)
- [`calc_equivalent_gmr`](@ref)
"""
function calc_gmd(co1::T, co2::U) where {T<:AbstractCablePart,U<:AbstractCablePart}

    if co1 isa WireArray
        coords1 = calc_wirearray_coords(co1.num_wires, co1.radius_wire, co1.radius_in)
        n1 = co1.num_wires
        r1 = co1.radius_wire
        s1 = pi * r1^2
    else
        coords1 = [(0.0, 0.0)]
        n1 = 1
        r1 = co1.radius_ext
        s1 = co1.cross_section
    end

    if co2 isa WireArray
        coords2 = calc_wirearray_coords(co2.num_wires, co2.radius_wire, co2.radius_in)
        n2 = co2.num_wires
        r2 = co2.radius_wire
        s2 = pi * r2^2
    else
        coords2 = [(0.0, 0.0)]
        n2 = 1
        r2 = co2.radius_ext
        s2 = co2.cross_section
    end

    log_sum = 0.0
    area_weights = 0.0

    for i in 1:n1
        for j in 1:n2
            # Pair-wise distances
            x1, y1 = coords1[i]
            x2, y2 = coords2[j]
            d_ij = sqrt((x1 - x2)^2 + (y1 - y2)^2)
            if d_ij > eps()
                # The GMD is computed as the Euclidean distance from center-to-center
                log_dij = log(d_ij)
            else
                # This means two concentric structures (solid/strip or tubular, tubular/strip or tubular, strip/strip or tubular)
                # In all cases the GMD is the outermost radius
                max(r1, r2)
                log_dij = log(max(r1, r2))
            end
            log_sum += (s1 * s2) * log_dij
            area_weights += (s1 * s2)
        end
    end
    return exp(log_sum / area_weights)
end
"""
$(TYPEDSIGNATURES)

Calculates the solenoid correction factor for magnetic permeability in insulated cables with helical conductors ([`WireArray`](@ref)), using the formula from Gudmundsdottir et al. [5743045](@cite):

```math
\\mu_{r, sol} = 1 + \\frac{2 \\pi^2 N^2 (r_{ins, ext}^2 - r_{con, ext}^2)}{\\log(r_{ins, ext}/r_{con, ext})}
```

where:
- ``N`` is the number of turns per unit length.
- ``r_{con, ext}`` is the conductor external radius.
- ``r_{ins, ext}`` is the insulator external radius.

# Arguments

- `num_turns`: Number of turns per unit length \\[1/m\\].
- `radius_ext_con`: External radius of the conductor \\[m\\].
- `radius_ext_ins`: External radius of the insulator \\[m\\].

# Returns

- Correction factor for the insulator magnetic permeability \\[dimensionless\\].

# Examples

```julia
# Cable with 10 turns per meter, conductor radius 5 mm, insulator radius 10 mm
correction = $(FUNCTIONNAME)(10, 0.005, 0.01)  # Expected output: > 1.0 [dimensionless]

# Non-helical cable (straight conductor)
correction = $(FUNCTIONNAME)(NaN, 0.005, 0.01)  # Expected output: 1.0 [dimensionless]
```
"""
function calc_solenoid_correction(
    num_turns::T,
    radius_ext_con::T,
    radius_ext_ins::T,
) where {T<:REALSCALAR}
    if isnan(num_turns)
        return 1.0
    else
        return 1.0 +
               2 * num_turns^2 * pi^2 * (radius_ext_ins^2 - radius_ext_con^2) /
               log(radius_ext_ins / radius_ext_con)
    end
end

function calc_solenoid_correction(
    num_turns,
    radius_ext_con,
    radius_ext_ins,
)
    T = resolve_T(num_turns, radius_ext_con, radius_ext_ins)
    return calc_solenoid_correction(
        coerce_to_T(num_turns, T),
        coerce_to_T(radius_ext_con, T),
        coerce_to_T(radius_ext_ins, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the equivalent resistivity of a solid tubular conductor, using the formula [916943](@cite):

```math
\\rho_{eq} = R_{eq} S_{eff} = R_{eq} \\pi (r_{ext}^2 - r_{in}^2)
```

where ``S_{eff}`` is the effective cross-sectional area of the tubular conductor.

# Arguments

- `R`: Resistance of the conductor \\[Ω\\].
- `radius_ext_con`: External radius of the tubular conductor \\[m\\].
- `radius_in_con`: Internal radius of the tubular conductor \\[m\\].

# Returns

- Equivalent resistivity of the tubular conductor \\[Ω·m\\].

# Examples

```julia
rho_eq = $(FUNCTIONNAME)(0.01, 0.02, 0.01)  # Expected output: ~9.42e-4 [Ω·m]
```
"""
function calc_equivalent_rho(R::T, radius_ext_con::T, radius_in_con::T) where {T<:REALSCALAR}
    eff_conductor_area = π * (radius_ext_con^2 - radius_in_con^2)
    return R * eff_conductor_area
end

function calc_equivalent_rho(R, radius_ext_con, radius_in_con)
    T = resolve_T(R, radius_ext_con, radius_in_con)
    return calc_equivalent_rho(
        coerce_to_T(R, T),
        coerce_to_T(radius_ext_con, T),
        coerce_to_T(radius_in_con, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the equivalent permittivity for a coaxial cable insulation, using the formula [916943](@cite):

```math
\\varepsilon_{eq} = \\frac{C_{eq} \\log(\\frac{r_{ext}}{r_{in}})}{2\\pi \\varepsilon_0}
```

where ``\\varepsilon_0`` is the permittivity of free space.

# Arguments

- `C_eq`: Equivalent capacitance of the insulation \\[F/m\\].
- `radius_ext`: External radius of the insulation \\[m\\].
- `radius_in`: Internal radius of the insulation \\[m\\].

# Returns

- Equivalent relative permittivity of the insulation \\[dimensionless\\].

# Examples

```julia
eps_eq = $(FUNCTIONNAME)(1e-10, 0.01, 0.005)  # Expected output: ~2.26 [dimensionless]
```

# See also
- [`ε₀`](@ref)
"""
function calc_equivalent_eps(C_eq::T, radius_ext::T, radius_in::T) where {T<:REALSCALAR}
    return (C_eq * log(radius_ext / radius_in)) / (2 * pi) / ε₀
end

function calc_equivalent_eps(C_eq, radius_ext, radius_in)
    T = resolve_T(C_eq, radius_ext, radius_in)
    return calc_equivalent_eps(
        coerce_to_T(C_eq, T),
        coerce_to_T(radius_ext, T),
        coerce_to_T(radius_in, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the equivalent loss factor (tangent) of a dielectric material:

```math
\\tan \\delta = \\frac{G_{eq}}{\\omega \\cdot C_{eq}}
```

where ``\\tan \\delta`` is the loss factor (tangent).

# Arguments

- `G_eq`: Equivalent conductance of the material \\[S·m\\].
- `C_eq`: Equivalent capacitance of the material \\[F/m\\].
- `ω`: Angular frequency \\[rad/s\\].

# Returns

- Equivalent loss factor of the dielectric material \\[dimensionless\\].

# Examples

```julia
loss_factor = $(FUNCTIONNAME)(1e-8, 1e-10, 2π*50)  # Expected output: ~0.0318 [dimensionless]
```
"""
function calc_equivalent_lossfact(G_eq::T, C_eq::T, ω::T) where {T<:REALSCALAR}
    return G_eq / (ω * C_eq)
end

function calc_equivalent_lossfact(G_eq, C_eq, ω)
    T = resolve_T(G_eq, C_eq, ω)
    return calc_equivalent_lossfact(
        coerce_to_T(G_eq, T),
        coerce_to_T(C_eq, T),
        coerce_to_T(ω, T),
    )
end

"""
$(TYPEDSIGNATURES)

Calculates the effective conductivity of a dielectric material from the known conductance (related to the loss factor ``\\tan \\delta``) via [916943](@cite) [Karmokar2025](@cite) [4389974](@cite):

```math
\\sigma_{eq} = \\frac{G_{eq}}{2\\pi} \\log(\\frac{r_{ext}}{r_{in}})
```
where ``\\sigma_{eq} = \\frac{1}{\\rho_{eq}}`` is the conductivity of the dielectric/semiconducting material, ``G_{eq}`` is the shunt conductance per unit length, ``r_{in}`` is the internal radius, and ``r_{ext}`` is the external radius of the coaxial structure.

# Arguments

- `G_eq`: Equivalent conductance of the material \\[S·m\\].
- `radius_in`: Internal radius of the coaxial structure \\[m\\].
- `radius_ext`: External radius of the coaxial structure \\[m\\].

# Returns

- Effective material conductivity per unit length \\[S·m\\].

# Examples

```julia
Geq = 2.7169e-9
sigma_eq = $(FUNCTIONNAME)(G_eq, radius_in, radius_ext)
```
"""
function calc_sigma_lossfact(G_eq::T, radius_in::T, radius_ext::T) where {T<:REALSCALAR}
    return G_eq * log(radius_ext / radius_in) / (2 * pi)
end

function calc_sigma_lossfact(G_eq, radius_in, radius_ext)
    T = resolve_T(G_eq, radius_in, radius_ext)
    return calc_sigma_lossfact(
        coerce_to_T(G_eq, T),
        coerce_to_T(radius_in, T),
        coerce_to_T(radius_ext, T),
    )
end

end # module BaseParams