"""
$(TYPEDEF)

Select the formulas for a coaxial line-parameters computation, with the longitudinal
wavenumber `Γ` \\[1/m\\] and matrix reductions as options.
Equal earth slots share one formula and its state at each frequency.

$(TYPEDFIELDS)
"""
struct CoaxialFormulation{
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
    "Formula for the surface impedances of each equivalent conductor."
    internal::I
    "Formula for the magnetic series-impedance term across each dielectric interval."
    insulation::D
    "Admittivity law for insulator layers."
    insulation_admittance::K
    "Admittivity law for semicon layers."
    semicon_admittance::S
    "Formula for the earth-return impedances of the cables."
    earth_impedance::Z
    "Formula for the earth-return potential coefficients of the cables."
    earth_admittance::Y
    "Formula for potential coefficients inside each cable."
    shunt::U
    "Temperature law for conductor resistivity."
    resistivity::R
    "Frequency law for air and earth material properties."
    properties::P
    "Options: the longitudinal wavenumber `Γ` \\[1/m\\] and the reductions."
    options::O
    function CoaxialFormulation(
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
            O <: FormulationOptions{CoaxialFormulation}
    }
        validate(earth_impedance, impedance!)
        validate(earth_admittance, potential!)
        validate((earth_impedance, earth_admittance), CoaxialFormulation)
        return new{I, D, K, S, Z, Y, U, R, P, O}(
            internal, insulation, insulation_admittance,
            semicon_admittance, earth_impedance, earth_admittance, shunt, resistivity,
            properties, options)
    end
end

function defaults(::Type{CoaxialFormulation})
    return (Γ = complex(0.0), bundle = true, kron = true, transposition = false)
end

"""
$(TYPEDSIGNATURES)

Build a coaxial formulation. Each formula keyword accepts `:id`,
`:id => (key = value, ...)`, or a formula of its family. Other keywords set the formulation
options. Earth slots with the same identifier must hold equal formulas.

# Keywords

- `internal`: conductor surface-impedance formula. Default: `:schelkunoff1934`.
- `insulation`: magnetic dielectric-interval formula. Default: `:ametani1980`.
- `insulation_admittance`: admittivity law for insulator layers. Default: `:lossless`.
- `semicon_admittance`: admittivity law for semicon layers. Default: `:lossless`.
- `earth_impedance`: earth-return impedance formula. Default: `:unified`.
- `earth_admittance`: earth-return potential-coefficient formula. Default: `:unified`.
- `shunt`: formula for potential coefficients inside cables. Default: `:series`.
- `resistivity`: conductor resistivity law. Default: `:linear`.
- `properties`: air and earth properties law. Default: `:constant`.
- `Γ`: longitudinal wavenumber \\[1/m\\], one value or one per frequency. Default: `0`.
- `bundle`: merge conductors with the same positive connection index. Default: `true`.
- `kron`: remove conductors with connection index zero by [`eliminate!`](@ref). Default:
  `true`.
- `transposition`: apply ideal cyclic transposition. Default: `false`.

# Returns

- A [`CoaxialFormulation`](@ref) with typed formula slots and options.

# Errors

- Throws `ArgumentError` for an unknown formula identifier or option, or invalid integration
  controls.
- Throws `ArgumentError` for different formulas of one identifier in the earth slots.
- Throws `DomainError` for an earth formula that does not give the quantity of its slot.
- Throws `DomainError` when `Γ` is empty or nonfinite.
"""
function CoaxialFormulation(;
        internal = :schelkunoff1934, insulation = :ametani1980,
        insulation_admittance = :lossless, semicon_admittance = :lossless,
        earth_impedance = :unified, earth_admittance = :unified, shunt = :series,
        resistivity = :linear, properties = :constant, kw...)
    return CoaxialFormulation(select(InternalImpedance.Formula, internal),
        select(InsulationImpedance.Formula, insulation),
        select(InsulationAdmittance.Formula, insulation_admittance),
        select(SemiconAdmittance.Formula, semicon_admittance),
        select(EarthReturn.Formula, earth_impedance),
        select(EarthReturn.Formula, earth_admittance),
        select(Shunt.Formula, shunt), select(Resistivity.Formula, resistivity),
        select(EarthProperties.Formula, properties),
        FormulationOptions(CoaxialFormulation; kw...))
end

"""
$(TYPEDSIGNATURES)

Check that earth slots with one identifier share one formula. Its impedance and
potential coefficients depend on the same circumference currents.

# Arguments

- `earth`: pair of earth-return formulas with the same identifier.
- `CoaxialFormulation`: formulation type that consumes the slots.

# Returns

- `earth`, unchanged.

# Errors

- Throws `ArgumentError` when the two formulas differ.
"""
function validate(earth::Tuple{EarthReturn.Formula{ID}, EarthReturn.Formula{ID}},
        ::Type{CoaxialFormulation}) where {ID}
    first(earth) == last(earth) || throw(ArgumentError(
        "the earth slots hold two different :$ID formulas; a formula that gives Z and P " *
        "together is one formula in both slots"))
    return earth
end

# Formulas of two identifiers give each the quantity of its own slot.
function validate(earth::Tuple{EarthReturn.Formula, EarthReturn.Formula},
        ::Type{CoaxialFormulation})
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

Select the material laws for a finite-element line-parameters computation, with the
longitudinal wavenumber `Γ` \\[1/m\\] and matrix reductions as options. The Gmsh extension
supplies the field solver.

$(TYPEDFIELDS)
"""
struct FEM{
    K <: InsulationAdmittance.Formula,
    S <: SemiconAdmittance.Formula,
    R <: Resistivity.Formula,
    P <: EarthProperties.Formula,
    O <: FormulationOptions
} <: AbstractFormulation
    "Admittivity law for insulator layers."
    insulation_admittance::K
    "Admittivity law for semicon layers."
    semicon_admittance::S
    "Temperature law for conductor resistivity."
    resistivity::R
    "Frequency law for air and earth material properties."
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

Build a finite-element formulation. Each formula keyword accepts `:id`,
`:id => (key = value, ...)`, or a formula of its family. Other keywords set the formulation
options.

# Keywords

- `insulation_admittance`: admittivity law for insulator layers. Default: `:lossless`.
- `semicon_admittance`: admittivity law for semicon layers. Default: `:lossless`.
- `resistivity`: conductor resistivity law. Default: `:linear`.
- `properties`: air and earth properties law. Default: `:constant`.
- `Γ`: longitudinal wavenumber \\[1/m\\], one value or one per frequency. Default: `0`.
- `bundle`: merge conductors with the same positive connection index. Default: `true`.
- `kron`: remove conductors with connection index zero by [`eliminate!`](@ref). Default:
  `true`.
- `transposition`: apply ideal cyclic transposition. Default: `false`.

# Returns

- An [`FEM`](@ref) formulation with typed material laws and options.

# Errors

- Throws `ArgumentError` for an unknown formula identifier or option.
- Throws `DomainError` when `Γ` is empty or nonfinite.
"""
function FEM(; insulation_admittance = :lossless, semicon_admittance = :lossless,
        resistivity = :linear, properties = :constant, kw...)
    return FEM(select(InsulationAdmittance.Formula, insulation_admittance),
        select(SemiconAdmittance.Formula, semicon_admittance),
        select(Resistivity.Formula, resistivity),
        select(EarthProperties.Formula, properties), FormulationOptions(FEM; kw...))
end

# Γ is one finite value, or one finite value per frequency.
function validate(options::FormulationOptions{<:Union{CoaxialFormulation, FEM}})
    Γ = options.Γ
    !isempty(Γ) && all(isfinite, Γ) || throw(DomainError(Γ,
        "the option Γ must be finite: one value, or one value per frequency"))
    return options
end

"""
$(TYPEDSIGNATURES)

Check that a vector longitudinal wavenumber has one entry per frequency.

# Arguments

- `problem`: line-parameters problem.
- `formulation`: formulation with a scalar or vector `Γ` option \\[1/m\\].

# Returns

- `problem`, unchanged.

# Errors

- Throws `DimensionMismatch` when a `Γ` vector and the frequency vector have different
  lengths.
"""
function validate(problem::LineParametersProblem,
        formulation::Union{CoaxialFormulation, FEM})
    Γ = formulation.options.Γ
    ndims(Γ) == 0 || length(Γ) == length(problem.frequencies) ||
        throw(DimensionMismatch(
            "the option Γ has $(length(Γ)) values for $(length(problem.frequencies)) " *
            "frequencies"))
    return problem
end
