"""
$(TYPEDEF)

Select the formula of each part of the coaxial computation of line parameters, and the
reductions that it applies to the matrices. Each slot takes a formula of its family. Two
equal earth slots are one formula, which prepares its state once per frequency.

$(TYPEDFIELDS)
"""
struct LineParametersFormulation{
    I <: InternalImpedance.Formula,
    D <: InsulationImpedance.Formula,
    K <: InsulationAdmittance.Formula,
    S <: SemiconAdmittance.Formula,
    Z <: EarthReturn.Formula,
    Y <: EarthReturn.Formula,
    U <: Shunt.Formula,
    R <: Resistivity.Formula,
    P <: EarthProperties.Formula,
    O <: FormulationOptions
} <: AbstractFormulation
    "Surface impedances of each equivalent conductor."
    internal::I
    "Magnetic term of the series impedance across each dielectric interval."
    insulation::D
    "Admittivity of the insulating layers."
    insulation_admittance::K
    "Admittivity of the semiconducting layers."
    semicon_admittance::S
    "Earth-return impedance of the cables."
    earth_impedance::Z
    "Earth-return potential coefficients of the cables."
    earth_admittance::Y
    "Potential coefficients inside each cable, from the admittivities of its layers."
    shunt::U
    "Resistivity of the conductors at the operating temperature."
    resistivity::R
    "Resistivity, permittivity and permeability of the air and the earth at each frequency."
    properties::P
    "Reductions of the matrices: `bundle`, `kron` and `transposition`."
    options::O
    function LineParametersFormulation(
            internal::I, insulation::D, insulation_admittance::K,
            semicon_admittance::S, earth_impedance::Z, earth_admittance::Y, shunt::U,
            resistivity::R, properties::P,
            options::O) where {
            I <: InternalImpedance.Formula,
            D <: InsulationImpedance.Formula,
            K <: InsulationAdmittance.Formula,
            S <: SemiconAdmittance.Formula,
            Z <: EarthReturn.Formula,
            Y <: EarthReturn.Formula,
            U <: Shunt.Formula,
            R <: Resistivity.Formula,
            P <: EarthProperties.Formula,
            O <: FormulationOptions{LineParametersFormulation}
    }
        validate(earth_impedance, impedance!)
        validate(earth_admittance, potential!)
        return new{I, D, K, S, Z, Y, U, R, P, O}(
            internal, insulation, insulation_admittance,
            semicon_admittance, earth_impedance, earth_admittance, shunt, resistivity,
            properties, options)
    end
end

function defaults(::Type{LineParametersFormulation})
    return (bundle = true, kron = true, transposition = false)
end

# A slot takes a formula of its family, or the identifier of the formula, alone or paired
# with its options.
select(::Type{F}, id::Symbol) where {F <: InternalImpedance.Formula} = F(Val(id))
function select(::Type{F},
        selection::Pair{Symbol, <:NamedTuple}) where {F <: InternalImpedance.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: InternalImpedance.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: InsulationImpedance.Formula} = F(Val(id))
function select(::Type{F},
        selection::Pair{Symbol, <:NamedTuple}) where {F <: InsulationImpedance.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: InsulationImpedance.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: InsulationAdmittance.Formula} = F(Val(id))
function select(::Type{F},
        selection::Pair{Symbol, <:NamedTuple}) where {F <: InsulationAdmittance.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: InsulationAdmittance.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: SemiconAdmittance.Formula} = F(Val(id))
function select(::Type{F},
        selection::Pair{Symbol, <:NamedTuple}) where {F <: SemiconAdmittance.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: SemiconAdmittance.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: EarthReturn.Formula} = F(Val(id))
function select(
        ::Type{F}, selection::Pair{Symbol, <:NamedTuple}) where {F <: EarthReturn.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: EarthReturn.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: Shunt.Formula} = F(Val(id))
function select(::Type{F}, selection::Pair{Symbol, <:NamedTuple}) where {F <: Shunt.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: Shunt.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: Resistivity.Formula} = F(Val(id))
function select(
        ::Type{F}, selection::Pair{Symbol, <:NamedTuple}) where {F <: Resistivity.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: Resistivity.Formula} = formula

select(::Type{F}, id::Symbol) where {F <: EarthProperties.Formula} = F(Val(id))
function select(::Type{F},
        selection::Pair{Symbol, <:NamedTuple}) where {F <: EarthProperties.Formula}
    return F(Val(first(selection)); last(selection)...)
end
select(::Type{F}, formula::F) where {F <: EarthProperties.Formula} = formula

"""
$(TYPEDEF)

Select the laws of the finite-element computation of line parameters, which it shares with
[`LineParametersFormulation`](@ref), and its options: the longitudinal wavenumber `Γ` and
the reductions. The field solver of the Gmsh extension takes the place of the other
formulas.

$(TYPEDFIELDS)
"""
struct FEM{
    K <: InsulationAdmittance.Formula,
    S <: SemiconAdmittance.Formula,
    R <: Resistivity.Formula,
    P <: EarthProperties.Formula,
    O <: FormulationOptions
} <: AbstractFormulation
    "Admittivity of the insulating layers."
    insulation_admittance::K
    "Admittivity of the semiconducting layers."
    semicon_admittance::S
    "Resistivity of the conductors at the operating temperature."
    resistivity::R
    "Resistivity, permittivity and permeability of the air and the earth at each frequency."
    properties::P
    "Options: the longitudinal wavenumber `Γ` \\[1/m\\] and the reductions."
    options::O
    function FEM(insulation_admittance::K, semicon_admittance::S, resistivity::R,
            properties::P,
            options::O) where {
            K <: InsulationAdmittance.Formula,
            S <: SemiconAdmittance.Formula,
            R <: Resistivity.Formula,
            P <: EarthProperties.Formula,
            O <: FormulationOptions{FEM}
    }
        return new{K, S, R, P, O}(
            insulation_admittance, semicon_admittance, resistivity, properties, options)
    end
end

function defaults(::Type{FEM})
    return (Γ = complex(0.0), bundle = true, kron = true, transposition = false)
end

"""
$(TYPEDSIGNATURES)

Build the finite-element formulation. Each law keyword takes the identifier of a formula of
its family, an identifier with the options of its formula, or a formula. The other keywords
are the options.

# Keywords

- `insulation_admittance`: admittivity of the insulating layers. Default: `:lossless`.
- `semicon_admittance`: admittivity of the semiconducting layers. Default: `:lossless`.
- `resistivity`: resistivity of the conductors. Default: `:linear`.
- `properties`: properties of the air and the earth. Default: `:constant`.
- `Γ`: longitudinal wavenumber \\[1/m\\], one value or one per frequency. Default: `0`.
- `bundle`, `kron`, `transposition`: the reductions. Default: `true`, `true`, `false`.

# Errors

- Throws `ArgumentError` for an unknown formula identifier or option.
- Throws `DomainError` when `Γ` is empty or not finite.
"""
function FEM(; insulation_admittance = :lossless, semicon_admittance = :lossless,
        resistivity = :linear, properties = :constant, kw...)
    return FEM(select(InsulationAdmittance.Formula, insulation_admittance),
        select(SemiconAdmittance.Formula, semicon_admittance),
        select(Resistivity.Formula, resistivity),
        select(EarthProperties.Formula, properties), FormulationOptions(FEM; kw...))
end

# Γ is one finite value, or one finite value per frequency.
function validate(options::FormulationOptions{FEM})
    Γ = options.Γ
    !isempty(Γ) && all(isfinite, Γ) || throw(DomainError(Γ,
        "the FEM option Γ must be finite: one value, or one value per frequency"))
    return options
end
