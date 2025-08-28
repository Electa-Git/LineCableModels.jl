
"""
$(TYPEDEF)

Represents a [`CableComponent`](@ref), i.e. a group of [`AbstractCablePart`](@ref) objects, with the equivalent geometric and material properties:

$(TYPEDFIELDS)

!!! info "Definition & application"
	Cable components operate as containers for multiple cable parts, allowing the calculation of effective electromagnetic (EM) properties (``\\sigma, \\varepsilon, \\mu``). This is performed by transforming the physical objects within the [`CableComponent`](@ref) into one equivalent coaxial homogeneous structure comprised of one conductor and one insulator, each one represented by effective [`Material`](@ref) types stored in `conductor_props` and `insulator_props` fields.

	The effective properties approach is widely adopted in EMT-type simulations, and involves locking the internal and external radii of the conductor and insulator parts, respectively, and calculating the equivalent EM properties in order to match the previously determined values of R, L, C and G [916943](@cite) [1458878](@cite).

	In applications, the [`CableComponent`](@ref) type is mapped to the main cable structures described in manufacturer datasheets, e.g., core, sheath, armor and jacket.
"""
mutable struct CableComponent{T<:REALSCALAR}
    "Cable component identification (e.g. core/sheath/armor)."
    id::String
    "The conductor group containing all conductive parts."
    conductor_group::ConductorGroup{T}
    "Effective properties of the equivalent coaxial conductor."
    conductor_props::Material{T}
    "The insulator group containing all insulating parts."
    insulator_group::InsulatorGroup{T}
    "Effective properties of the equivalent coaxial insulator."
    insulator_props::Material{T}

    @doc """
    $(TYPEDSIGNATURES)

    Initializes a [`CableComponent`](@ref) object based on its constituent conductor and insulator groups. The constructor performs the following sequence of steps:

    1.  Validate that the conductor and insulator groups have matching radii at their interface.
    2.  Obtain the lumped-parameter values (R, L, C, G) from the conductor and insulator groups, which are computed within their respective constructors.
    3.  Calculate the correction factors and equivalent electromagnetic properties of the conductor and insulator groups:


    | Quantity | Symbol | Function |
    |----------|--------|----------|
    | Resistivity (conductor) | ``\\rho_{con}`` | [`calc_equivalent_rho`](@ref) |
    | Permeability (conductor) | ``\\mu_{con}`` | [`calc_equivalent_mu`](@ref) |
    | Resistivity (insulator) | ``\\rho_{ins}`` | [`calc_sigma_lossfact`](@ref) |
    | Permittivity (insulation) | ``\\varepsilon_{ins}`` | [`calc_equivalent_eps`](@ref) |
    | Permeability (insulation) | ``\\mu_{ins}`` | [`calc_solenoid_correction`](@ref) |

    # Arguments

    - `id`: Cable component identification (e.g. core/sheath/armor).
    - `conductor_group`: The conductor group containing all conductive parts.
    - `insulator_group`: The insulator group containing all insulating parts.

    # Returns

    A [`CableComponent`](@ref) instance with calculated equivalent properties:

    - `id::String`: Cable component identification.
    - `conductor_group::ConductorGroup{T}`: The conductor group containing all conductive parts.
    - `conductor_props::Material{T}`: Effective properties of the equivalent coaxial conductor.
        * `rho`: Resistivity \\[Ω·m\\].
        * `eps_r`: Relative permittivity \\[dimensionless\\].
        * `mu_r`: Relative permeability \\[dimensionless\\].
        * `T0`: Reference temperature \\[°C\\].
        * `alpha`: Temperature coefficient of resistivity \\[1/°C\\].
    - `insulator_group::InsulatorGroup{T}`: The insulator group containing all insulating parts.
    - `insulator_props::Material{T}`: Effective properties of the equivalent coaxial insulator.
        * `rho`: Resistivity \\[Ω·m\\].
        * `eps_r`: Relative permittivity \\[dimensionless\\].
        * `mu_r`: Relative permeability \\[dimensionless\\].
        * `T0`: Reference temperature \\[°C\\].
        * `alpha`: Temperature coefficient of resistivity \\[1/°C\\].

    # Examples

    ```julia
    conductor_group = ConductorGroup(...)
    insulator_group = InsulatorGroup(...)
    cable = $(FUNCTIONNAME)("component_id", conductor_group, insulator_group)  # Create cable component with base parameters @ 50 Hz
    ```

    # See also

    - [`calc_equivalent_rho`](@ref)
    - [`calc_equivalent_mu`](@ref)
    - [`calc_equivalent_eps`](@ref)
    - [`calc_sigma_lossfact`](@ref)
    - [`calc_solenoid_correction`](@ref)
    """
    function CableComponent{T}(
        id::String,
        conductor_group::ConductorGroup{T},
        insulator_group::InsulatorGroup{T},
    ) where {T<:REALSCALAR}

        # Geometry interface check (exact or approximately equal)
        if !(conductor_group.radius_ext == insulator_group.radius_in ||
             isapprox(conductor_group.radius_ext, insulator_group.radius_in))
            throw(ArgumentError("Conductor outer radius must match insulator inner radius."))
        end

        # Radii
        r1 = conductor_group.radius_in
        r2 = conductor_group.radius_ext
        r3 = insulator_group.radius_ext

        # 2) Conductor equivalents
        ρ_con = calc_equivalent_rho(conductor_group.resistance, r2, r1)
        μ_con = calc_equivalent_mu(conductor_group.gmr, r2, r1)
        α_con = conductor_group.alpha
        θ_con = conductor_group.layers[1].temperature
        conductor_props = Material{T}(ρ_con, T(0), μ_con, θ_con, α_con)

        # 3) Insulator equivalents (use already-aggregated C and G)
        C_eq = insulator_group.shunt_capacitance
        G_eq = insulator_group.shunt_conductance
        ε_ins = calc_equivalent_eps(C_eq, r3, r2)
        σ_ins = calc_sigma_lossfact(G_eq, r2, r3)
        ρ_ins = inv(σ_ins)               # safe if σ_ins ≠ 0 (your tests cover zero?)
        μ_ins = calc_solenoid_correction(conductor_group.num_turns, r2, r3)
        θ_ins = insulator_group.layers[1].temperature
        insulator_props = Material{T}(ρ_ins, ε_ins, μ_ins, θ_ins, T(0))

        return new{T}(
            id,
            conductor_group,
            conductor_props,
            insulator_group,
            insulator_props,
        )
    end
end

"""
$(TYPEDSIGNATURES)

Weakly-typed constructor that infers the scalar type `T` from the two groups, coerces them if necessary, and calls the strict kernel.

# Arguments
- `id`: Cable component identification.
- `conductor_group`: The conductor group (any `ConductorGroup{S}`).
- `insulator_group`: The insulator group (any `InsulatorGroup{R}`).

# Returns
- A `CableComponent{T}` where `T` is the resolved scalar type.
"""
function CableComponent(
    id::String,
    conductor_group::ConductorGroup,
    insulator_group::InsulatorGroup,
)
    # Resolve target T from the two groups (honors Measurements, etc.)
    T = resolve_T(conductor_group, insulator_group)

    # Coerce groups to T (identity if already T)
    cgT = coerce_to_T(conductor_group, T)
    igT = coerce_to_T(insulator_group, T)

    return CableComponent{T}(id, cgT, igT)
end

include("cablecomponent/base.jl")

