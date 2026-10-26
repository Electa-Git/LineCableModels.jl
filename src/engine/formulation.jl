"""
$(TYPEDEF)

Select the formula of each part of the coaxial computation of line parameters. The options
are the longitudinal wavenumber `Γ` and the reductions of the matrices. Each slot takes a
formula of its family. Equal earth slots are one formula, which prepares its state once per
frequency.

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
    "Options: the longitudinal wavenumber `Γ` \\[1/m\\] and the reductions."
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
    return (Γ = complex(0.0), bundle = true, kron = true, transposition = false)
end

"""
$(TYPEDSIGNATURES)

Build the coaxial formulation. Each slot keyword takes the identifier of a formula of its
family, an identifier with the options of its formula, or a formula. Equal earth slots are
one formula, which prepares its state once per frequency. Formulas of one identifier in the
earth slots must be equal. The other keywords are the options.

# Keywords

- `internal`: surface impedances of the conductors. Default: `:schelkunoff1934`.
- `insulation`: magnetic term across the dielectric intervals. Default: `:ametani1980`.
- `insulation_admittance`: admittivity of the insulating layers. Default: `:lossless`.
- `semicon_admittance`: admittivity of the semicon layers. Default: `:lossless`.
- `earth_impedance`: earth-return impedance. Default: `:unified`.
- `earth_admittance`: earth-return potential coefficients. Default: `:unified`.
- `shunt`: potential coefficients inside the cables. Default: `:series`.
- `resistivity`: resistivity of the conductors. Default: `:linear`.
- `properties`: properties of the air and the earth. Default: `:constant`.
- `Γ`: longitudinal wavenumber \\[1/m\\], one value or one per frequency. Default: `0`.
- `bundle`, `kron`, `transposition`: the reductions. Default: `true`, `true`, `false`.

# Errors

- Throws `ArgumentError` for an unknown formula identifier or option.
- Throws `ArgumentError` for different formulas of one identifier in the earth slots.
- Throws `DomainError` for an earth formula that does not give the quantity of its slot.
- Throws `DomainError` when `Γ` is empty or not finite.
"""
function LineParametersFormulation(;
        internal = :schelkunoff1934, insulation = :ametani1980,
        insulation_admittance = :lossless, semicon_admittance = :lossless,
        earth_impedance = :unified, earth_admittance = :unified, shunt = :series,
        resistivity = :linear, properties = :constant, kw...)
    impedance_formula = select(EarthReturn.Formula, earth_impedance)
    admittance_formula = select(EarthReturn.Formula, earth_admittance)
    validate((impedance_formula, admittance_formula), LineParametersFormulation)
    return LineParametersFormulation(select(InternalImpedance.Formula, internal),
        select(InsulationImpedance.Formula, insulation),
        select(InsulationAdmittance.Formula, insulation_admittance),
        select(SemiconAdmittance.Formula, semicon_admittance), impedance_formula,
        impedance_formula == admittance_formula ? impedance_formula : admittance_formula,
        select(Shunt.Formula, shunt), select(Resistivity.Formula, resistivity),
        select(EarthProperties.Formula, properties),
        FormulationOptions(LineParametersFormulation; kw...))
end

"""
$(TYPEDSIGNATURES)

Check the earth slots of a coaxial formulation whose formulas share an identifier. A formula
of that identifier gives the impedance and the potential coefficients from one set of
currents on the circumferences of the cables. Both slots must take the same formula.

# Errors

- Throws `ArgumentError` when the two formulas differ.
"""
function validate(earth::Tuple{EarthReturn.Formula{ID}, EarthReturn.Formula{ID}},
        ::Type{LineParametersFormulation}) where {ID}
    first(earth) == last(earth) || throw(ArgumentError(
        "the earth slots hold two different :$ID formulas; a formula that gives Z and P " *
        "together is one formula in both slots"))
    return earth
end

# Formulas of two identifiers give each the quantity of its own slot.
function validate(earth::Tuple{EarthReturn.Formula, EarthReturn.Formula},
        ::Type{LineParametersFormulation})
    return earth
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
function validate(options::FormulationOptions{<:Union{LineParametersFormulation, FEM}})
    Γ = options.Γ
    !isempty(Γ) && all(isfinite, Γ) || throw(DomainError(Γ,
        "the option Γ must be finite: one value, or one value per frequency"))
    return options
end

"""
$(TYPEDSIGNATURES)

Check that `formulation` can compute `problem`: a vector `Γ` has one value per frequency.

# Errors

- Throws `DimensionMismatch` for a `Γ` vector of another length than the frequencies.
"""
function validate(problem::LineParametersProblem,
        formulation::Union{LineParametersFormulation, FEM})
    Γ = formulation.options.Γ
    ndims(Γ) == 0 || length(Γ) == length(problem.frequencies) ||
        throw(DimensionMismatch(
            "the option Γ has $(length(Γ)) values for $(length(problem.frequencies)) " *
            "frequencies"))
    return problem
end
