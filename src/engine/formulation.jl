"""
$(TYPEDEF)

Select the formulas for a coaxial line-parameters computation, with the longitudinal
wavenumber `Γ` \\[1/m\\] and matrix reductions as options.
Equal earth formulas share one block at each frequency.

$(TYPEDFIELDS)
"""
struct CoaxialFormulation{
    I <: InternalImpedance.Formula,
    D <: InsulationImpedance.Formula,
    K <: InsulationAdmittance.Formula,
    S <: SemiconAdmittance.Formula,
    Z <: NamedTuple{(:air, :earth, :mixed), <:Tuple{Vararg{EarthReturn.Formula}}},
    Y <: NamedTuple{(:air, :earth, :mixed), <:Tuple{Vararg{EarthReturn.Formula}}},
    U <: Shunt.Formula,
    R <: Resistivity.Formula,
    P <: EarthProperties.Formula,
    E <: EarthProperties.Formula,
    L <: Tuple{Vararg{AbstractFormulation}},
    O <: FormulationOptions
} <: AbstractFormulation
    "Formula for the surface impedances of each equivalent conductor."
    internal_impedance::I
    "Formula for the magnetic series-impedance term across each dielectric interval."
    insulation_impedance::D
    "Admittivity law for insulator layers."
    insulation_admittance::K
    "Admittivity law for semicon layers."
    semicon_admittance::S
    "Earth-return impedance formula of each domain of a pair, `(air, earth, mixed)`."
    earth_impedance::Z
    "Earth-return potential-coefficient formula of each domain of a pair, `(air, earth, mixed)`."
    earth_admittance::Y
    "Formula for potential coefficients inside each cable."
    shunt::U
    "Temperature law for conductor resistivity."
    resistivity::R
    "Frequency law for air and earth material properties."
    earth_properties::P
    "Model of the equivalent earth of the layers below the air."
    equivalent_earth::E
    "Distinct formulas of the slots, computed once at construction."
    formulas::L
    "Options: the longitudinal wavenumber `Γ` \\[1/m\\] and the reductions."
    options::O
    function CoaxialFormulation(
            internal_impedance::I, insulation_impedance::D, insulation_admittance::K,
            semicon_admittance::S, earth_impedance::Z, earth_admittance::Y, shunt::U,
            resistivity::R, earth_properties::P, equivalent_earth::E,
            options::O) where {
            I <: InternalImpedance.Formula,
            D <: InsulationImpedance.Formula,
            K <: InsulationAdmittance.Formula,
            S <: SemiconAdmittance.Formula,
            Z <: NamedTuple{(:air, :earth, :mixed), <:Tuple{Vararg{EarthReturn.Formula}}},
            Y <: NamedTuple{(:air, :earth, :mixed), <:Tuple{Vararg{EarthReturn.Formula}}},
            U <: Shunt.Formula,
            R <: Resistivity.Formula,
            P <: EarthProperties.Formula,
            E <: EarthProperties.Formula,
            O <: FormulationOptions{CoaxialFormulation}
    }
        map(formula -> validate(formula, impedance!), earth_impedance)
        map(formula -> validate(formula, potential!), earth_admittance)
        validate((Tuple(earth_impedance)..., Tuple(earth_admittance)...), CoaxialFormulation)
        validate(earth_properties, EarthLayer)
        validate(equivalent_earth, EarthModel)
        formulas = foldl(
            (found, formula) -> any(f -> f == formula, found) ? found : (found..., formula),
            (internal_impedance, insulation_impedance,
                insulation_admittance, semicon_admittance,
                Tuple(earth_impedance)..., Tuple(earth_admittance)..., shunt, resistivity,
                earth_properties, equivalent_earth);
            init = ())
        return new{I, D, K, S, Z, Y, U, R, P, E, typeof(formulas), O}(
            internal_impedance, insulation_impedance, insulation_admittance, semicon_admittance,
            earth_impedance, earth_admittance, shunt, resistivity, earth_properties,
            equivalent_earth, formulas, options)
    end
end

function defaults(::Type{CoaxialFormulation})
    return (Γ = complex(0.0), bundle = true, kron = true, transposition = false)
end

"""
$(TYPEDSIGNATURES)

Build a coaxial formulation. Each formula keyword accepts `:id`,
`:id => (key = value, ...)`, a formula of its family, or `:default`, the default of its slot.
Other keywords set the formulation options. Earth slots with the same identifier must hold
equal formulas.

# Keywords

- `internal_impedance`: conductor surface-impedance formula. Default: `:schelkunoff1934`.
- `insulation_impedance`: magnetic dielectric-interval formula. Default: `:ametani1980`.
- `insulation_admittance`: admittivity law for insulator layers. Default: `:lossless`.
- `semicon_admittance`: admittivity law for semicon layers. Default: `:lossless`.
- `earth_impedance`: earth-return impedance formula, or a record `(air, earth, mixed)` of them.
  A domain left out takes the default. Default: `:martins2026`.
- `earth_admittance`: earth-return potential-coefficient formula, or a record
  `(air, earth, mixed)` of them. A domain left out takes the default. Default: `:martins2026`.
- `shunt`: formula for potential coefficients inside cables. Default: `:series`.
- `resistivity`: conductor resistivity law. Default: `:linear`.
- `earth_properties`: air and earth properties law. Default: `:constant`.
- `equivalent_earth`: model of the equivalent earth of the layers below the air. Default:
  `:bottommost`.
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
- Throws `DomainError` for a model of the equivalent earth in `earth_properties`, or a law in
  `equivalent_earth`.
- Throws `DomainError` when `Γ` is empty or nonfinite.
"""
function CoaxialFormulation(;
        internal_impedance = :default, insulation_impedance = :default,
        insulation_admittance = :default, semicon_admittance = :default,
        earth_impedance = :default, earth_admittance = :default, shunt = :default,
        resistivity = :default, earth_properties = :default, equivalent_earth = :default, kw...)
    return CoaxialFormulation(
        select(CoaxialFormulation, Val(:internal_impedance), internal_impedance),
        select(CoaxialFormulation, Val(:insulation_impedance), insulation_impedance),
        select(CoaxialFormulation, Val(:insulation_admittance), insulation_admittance),
        select(CoaxialFormulation, Val(:semicon_admittance), semicon_admittance),
        select(CoaxialFormulation, Val(:earth_impedance), earth_impedance),
        select(CoaxialFormulation, Val(:earth_admittance), earth_admittance),
        select(CoaxialFormulation, Val(:shunt), shunt),
        select(CoaxialFormulation, Val(:resistivity), resistivity),
        select(CoaxialFormulation, Val(:earth_properties), earth_properties),
        select(CoaxialFormulation, Val(:equivalent_earth), equivalent_earth),
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

# Places of the earth slots with one identifier must store equal formulas.
function validate(places::NTuple{6, EarthReturn.Formula}, ::Type{CoaxialFormulation})
    for a in places, b in places

        validate((a, b), CoaxialFormulation)
    end
    return places
end

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
    earth_properties::P
    "Options: the longitudinal wavenumber `Γ` \\[1/m\\] and the reductions."
    options::O
    function FEM(insulation_admittance::K, semicon_admittance::S, resistivity::R,
            earth_properties::P,
            options::O) where {
            K <: InsulationAdmittance.Formula,
            S <: SemiconAdmittance.Formula,
            R <: Resistivity.Formula,
            P <: EarthProperties.Formula,
            O <: FormulationOptions{FEM}
    }
        validate(earth_properties, EarthLayer)
        return new{K, S, R, P, O}(
            insulation_admittance, semicon_admittance, resistivity,
            earth_properties, options)
    end
end

function defaults(::Type{FEM})
    return (Γ = complex(0.0), bundle = true, kron = true, transposition = false)
end

"""
$(TYPEDSIGNATURES)

Build a finite-element formulation. Each formula keyword accepts `:id`,
`:id => (key = value, ...)`, a formula of its family, or `:default`, the default of its slot.
Other keywords set the formulation options.

# Keywords

- `insulation_admittance`: admittivity law for insulator layers. Default: `:lossless`.
- `semicon_admittance`: admittivity law for semicon layers. Default: `:lossless`.
- `resistivity`: conductor resistivity law. Default: `:linear`.
- `earth_properties`: air and earth properties law. Default: `:constant`.
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
function FEM(; insulation_admittance = :default, semicon_admittance = :default,
        resistivity = :default, earth_properties = :default, kw...)
    return FEM(select(FEM, Val(:insulation_admittance), insulation_admittance),
        select(FEM, Val(:semicon_admittance), semicon_admittance),
        select(FEM, Val(:resistivity), resistivity),
        select(FEM, Val(:earth_properties), earth_properties), FormulationOptions(FEM; kw...))
end

# A Symbol lifts to its tag and a pair to its tag with options. A formula stands.
function select(
        ::Type{F}, slot::Val, id::Symbol) where {F <: Union{CoaxialFormulation, FEM}}
    return select(F, slot, Val(id))
end
function select(::Type{F}, slot::Val,
        selection::Pair{Symbol, <:NamedTuple}) where {F <: Union{CoaxialFormulation, FEM}}
    return select(F, slot, Val(first(selection)); last(selection)...)
end
function select(::Type{F}, ::Val,
        formula::AbstractFormulation) where {F <: Union{CoaxialFormulation, FEM}}
    return formula
end
# An identifier lifts through the family of the slot, which the struct declares.
function select(::Type{F}, ::Val{S}, ::Val{ID};
        kw...) where {F <: Union{CoaxialFormulation, FEM}, S, ID}
    return fieldtype(F, S)(Val(ID); kw...)
end

# The default formula of each slot, the one place it is read and changed.
function select(::Type{CoaxialFormulation}, ::Val{:internal_impedance}, ::Val{:default}; kw...)
    return InternalImpedance.Formula(Val(:schelkunoff1934); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:insulation_impedance}, ::Val{:default}; kw...)
    return InsulationImpedance.Formula(Val(:ametani1980); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:insulation_admittance}, ::Val{:default}; kw...)
    return InsulationAdmittance.Formula(Val(:lossless); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:semicon_admittance}, ::Val{:default}; kw...)
    return SemiconAdmittance.Formula(Val(:lossless); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_impedance}, ::Val{:default}; kw...)
    return select(CoaxialFormulation, Val(:earth_impedance), Val(:martins2026); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_admittance}, ::Val{:default}; kw...)
    return select(CoaxialFormulation, Val(:earth_admittance), Val(:martins2026); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:shunt}, ::Val{:default}; kw...)
    return Shunt.Formula(Val(:series); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:resistivity}, ::Val{:default}; kw...)
    return Resistivity.Formula(Val(:linear); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_properties}, ::Val{:default}; kw...)
    return EarthProperties.Formula(Val(:constant); kw...)
end
function select(::Type{CoaxialFormulation}, ::Val{:equivalent_earth}, ::Val{:default}; kw...)
    return EarthProperties.Formula(Val(:bottommost); kw...)
end

function select(::Type{FEM}, ::Val{:insulation_admittance}, ::Val{:default}; kw...)
    return InsulationAdmittance.Formula(Val(:lossless); kw...)
end
function select(::Type{FEM}, ::Val{:semicon_admittance}, ::Val{:default}; kw...)
    return SemiconAdmittance.Formula(Val(:lossless); kw...)
end
function select(::Type{FEM}, ::Val{:resistivity}, ::Val{:default}; kw...)
    return Resistivity.Formula(Val(:linear); kw...)
end
function select(::Type{FEM}, ::Val{:earth_properties}, ::Val{:default}; kw...)
    return EarthProperties.Formula(Val(:constant); kw...)
end

# An earth slot stores one formula per domain. One selection fills the three domains. A
# selection per domain leaves the others to the default.
function select(::Type{CoaxialFormulation}, ::Val{:earth_impedance}, ::Val{ID}; kw...) where {ID}
    formula = EarthReturn.Formula(Val(ID); kw...)
    return (air = formula, earth = formula, mixed = formula)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_impedance}, formula::EarthReturn.Formula)
    return (air = formula, earth = formula, mixed = formula)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_impedance}, selection::NamedTuple)
    return merge(select(CoaxialFormulation, Val(:earth_impedance), Val(:default)),
        map(domain -> select(CoaxialFormulation, Val(:earth_impedance), domain).air, selection))
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_admittance}, ::Val{ID}; kw...) where {ID}
    formula = EarthReturn.Formula(Val(ID); kw...)
    return (air = formula, earth = formula, mixed = formula)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_admittance}, formula::EarthReturn.Formula)
    return (air = formula, earth = formula, mixed = formula)
end
function select(::Type{CoaxialFormulation}, ::Val{:earth_admittance}, selection::NamedTuple)
    return merge(select(CoaxialFormulation, Val(:earth_admittance), Val(:default)),
        map(domain -> select(CoaxialFormulation, Val(:earth_admittance), domain).air, selection))
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
