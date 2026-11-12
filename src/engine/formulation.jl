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
    "Frequency law for the properties of each earth layer."
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
    return (Γ = complex(0.0), bundle = true, kron = true, transposition = false, order = 0)
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
- `earth_properties`: frequency law of each earth layer. Default: `:constant`.
- `equivalent_earth`: model of the equivalent earth of the layers below the air. Default:
  `:bottommost`.
- `Γ`: longitudinal wavenumber \\[1/m\\], one value or one per frequency. Default: `0`.
- `bundle`: merge conductors with the same positive connection index. Default: `true`.
- `kron`: remove conductors with connection index zero by [`eliminate!`](@ref). Default:
  `true`.
- `transposition`: apply ideal cyclic transposition. Default: `false`.
- `order`: highest harmonic order of the field around each member of an interior that the
  recipes keep. Default: `0`.

# Returns

- A [`CoaxialFormulation`](@ref) with typed formula slots and options.

# Notes

At `order = 0` the recipes keep the field that is uniform around each conductor. That field
leaves out the proximity of the members to each other and to their bore. A system of
concentric cables does not read `order`.

# Errors

- Throws `ArgumentError` for an unknown formula identifier or option, or invalid integration
  controls.
- Throws `ArgumentError` for different formulas of one identifier in the earth slots.
- Throws `DomainError` for an earth formula that does not give the quantity of its slot.
- Throws `DomainError` for a model of the equivalent earth in `earth_properties`, or a law in
  `equivalent_earth`.
- Throws `DomainError` when `Γ` is empty or nonfinite.
- Throws `DomainError` when `order` is negative.
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
    E <: EarthProperties.Formula,
    L <: Tuple{Vararg{AbstractFormulation}},
    O <: FormulationOptions
} <: AbstractFormulation
    "Admittivity law for insulator layers."
    insulation_admittance::K
    "Admittivity law for semicon layers."
    semicon_admittance::S
    "Temperature law for conductor resistivity."
    resistivity::R
    "Frequency law for the properties of each earth layer."
    earth_properties::P
    "Model of the equivalent earth of the layers below the air."
    equivalent_earth::E
    "Formulas of the slots, which a computation validates against its problem."
    formulas::L
    "Options: the longitudinal wavenumber `Γ` \\[1/m\\] and the reductions."
    options::O
    function FEM(insulation_admittance::K, semicon_admittance::S, resistivity::R,
            earth_properties::P, equivalent_earth::E,
            options::O) where {
            K <: InsulationAdmittance.Formula,
            S <: SemiconAdmittance.Formula,
            R <: Resistivity.Formula,
            P <: EarthProperties.Formula,
            E <: EarthProperties.Formula,
            O <: FormulationOptions{FEM}
    }
        validate(earth_properties, EarthLayer)
        validate(equivalent_earth, EarthModel)
        formulas = (insulation_admittance, semicon_admittance, resistivity,
            earth_properties, equivalent_earth)
        return new{K, S, R, P, E, typeof(formulas), O}(
            insulation_admittance, semicon_admittance, resistivity,
            earth_properties, equivalent_earth, formulas, options)
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
- `earth_properties`: frequency law of each earth layer. Default: `:constant`.
- `equivalent_earth`: model of the equivalent earth of the layers below the air. Default:
  `:bottommost`.
- `Γ`: longitudinal wavenumber \\[1/m\\], one value or one per frequency. Default: `0`.
- `bundle`: merge conductors with the same positive connection index. Default: `true`.
- `kron`: remove conductors with connection index zero by [`eliminate!`](@ref). Default:
  `true`.
- `transposition`: apply ideal cyclic transposition. Default: `false`.

# Returns

- An [`FEM`](@ref) formulation with typed material laws and options.

# Errors

- Throws `ArgumentError` for an unknown formula identifier or option.
- Throws `DomainError` for a model of the equivalent earth in `earth_properties`, or a law in
  `equivalent_earth`.
- Throws `DomainError` when `Γ` is empty or nonfinite.
"""
function FEM(; insulation_admittance = :default, semicon_admittance = :default,
        resistivity = :default, earth_properties = :default, equivalent_earth = :default, kw...)
    return FEM(select(FEM, Val(:insulation_admittance), insulation_admittance),
        select(FEM, Val(:semicon_admittance), semicon_admittance),
        select(FEM, Val(:resistivity), resistivity),
        select(FEM, Val(:earth_properties), earth_properties),
        select(FEM, Val(:equivalent_earth), equivalent_earth), FormulationOptions(FEM; kw...))
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
function select(::Type{FEM}, ::Val{:equivalent_earth}, ::Val{:default}; kw...)
    return EarthProperties.Formula(Val(:bottommost); kw...)
end

# An earth slot stores one formula per domain. One selection fills the three domains. A
# selection per domain leaves the others to the default.
function select(::Type{CoaxialFormulation},
        ::Union{Val{:earth_impedance}, Val{:earth_admittance}}, ::Val{ID}; kw...) where {ID}
    formula = EarthReturn.Formula(Val(ID); kw...)
    return (air = formula, earth = formula, mixed = formula)
end
function select(::Type{CoaxialFormulation},
        ::Union{Val{:earth_impedance}, Val{:earth_admittance}}, formula::EarthReturn.Formula)
    return (air = formula, earth = formula, mixed = formula)
end
function select(::Type{CoaxialFormulation},
        slot::Union{Val{:earth_impedance}, Val{:earth_admittance}}, selection::NamedTuple)
    return merge(select(CoaxialFormulation, slot, Val(:default)),
        map(place -> select(CoaxialFormulation, slot, place).air, selection))
end

# Γ is one finite value, or one finite value per frequency.
function validate(options::FormulationOptions{<:Union{CoaxialFormulation, FEM}})
    Γ = options.Γ
    !isempty(Γ) && all(isfinite, Γ) || throw(DomainError(Γ,
        "the option Γ must be finite: one value, or one value per frequency"))
    return options
end

# The highest order is a nonnegative integer, and the method of every formulation checks Γ.
function validate(options::FormulationOptions{CoaxialFormulation})
    options.order >= 0 || throw(DomainError(options.order,
        "the option order must be a nonnegative integer: the highest harmonic order of the field around a member that the recipes keep"))
    invoke(validate, Tuple{FormulationOptions{<:Union{CoaxialFormulation, FEM}}}, options)
    return options
end

"""
$(TYPEDSIGNATURES)

Check that a vector longitudinal wavenumber has one entry per frequency.

# Arguments

- `problem`: line-parameters problem.
- `options`: options of a formulation, with a scalar or vector `Γ` \\[1/m\\].

# Returns

- `problem`, unchanged.

# Errors

- Throws `DimensionMismatch` when a `Γ` vector and the frequency vector have different
  lengths.
"""
function validate(problem::LineParametersProblem,
        options::FormulationOptions{<:Union{CoaxialFormulation, FEM}})
    Γ = options.Γ
    ndims(Γ) == 0 || length(Γ) == length(problem.frequencies) ||
        throw(DimensionMismatch(
            "the option Γ has $(length(Γ)) values for $(length(problem.frequencies)) " *
            "frequencies"))
    return problem
end

"""
$(TYPEDSIGNATURES)

Check the options of the formulation against the frequencies, and the problem against each
formula of the formulation.

# Arguments

- `problem`: line-parameters problem.
- `formulation`: coaxial or finite-element formulation.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DimensionMismatch` when a `Γ` vector and the frequency vector have different lengths.
- Throws `DomainError` when a cable touches the interface for the earth-return formula
  `:martins2026`.
"""
function validate(problem::LineParametersProblem,
        formulation::Union{CoaxialFormulation, FEM})
    validate(problem, formulation.options)
    map(formula -> validate(problem, formula), formulation.formulas)
    return problem
end

"""
$(TYPEDSIGNATURES)

Check that the formula of each earth place treats the pairs of layers of the input, before the
loop. A pair of two conductors in the air takes the place `air`, a pair of two conductors below
the air takes `earth`, and any other pair takes `mixed`.

# Arguments

- `input`: numerical cable geometry with the layer of each cable.
- `formulation`: coaxial formulation with its earth-return formulas.

# Returns

- `input`, unchanged.

# Errors

- Throws `DomainError` for a place whose formula does not treat a pair of the input.
"""
function validate(input::CoaxialInput, formulation::CoaxialFormulation)
    present = sort!(unique(input.earth_layer))
    for places in (formulation.earth_impedance, formulation.earth_admittance),
        t in present, s in present
        if t == 1 && s == 1
            validate(places.air, t, s)
        elseif t > 1 && s > 1
            validate(places.earth, t, s)
        else
            validate(places.mixed, t, s)
        end
    end
    return input
end

"""
$(TYPEDSIGNATURES)

Check the interiors of the input for the recipes of the formulation at the frequencies of the
problem, before the loop. The check compares positions and radii on their `nominal` values.
Circles overlap when the distance of their centers is below the sum of their radii and not
`isapprox` to it.

The body of a member is its conductors and the layers of its interval out to the last one
whose medium differs from the filler's. A medium differs by its relative permeability, or by
its admittivity at a frequency of the problem. The body of the bore begins at the innermost
interior layer outside the filler whose medium differs, else at the bore's inner radius.
Above order 0, the expansions about the members' axes converge outside these bodies. It warns
once for each member with members of its own above order 0, whose interior it computes as
shielded by that member's conductor.

# Arguments

- `input`: numerical input of the cables, with the members of each assembly.
- `formulation`: coaxial formulation, with the option `order` and its formulas.
- `problem`: line-parameters problem, with its frequencies and temperature.

# Returns

- `input`, unchanged.

# Errors

- Throws `ArgumentError` when an assembly has members and a conductor, and no region around
  the members.
- Throws `ArgumentError` when `order ≥ 1` and the bodies of two members overlap, or a
  member's body overlaps the bore's.
- Throws `ArgumentError` when `order ≥ 1`, an interior has members, and the formula of the
  slot `internal_impedance` gives the surface impedances at order 0 alone.
"""
function validate(input::CoaxialInput, formulation::CoaxialFormulation,
        problem::LineParametersProblem)
    N = formulation.options.order
    θ = problem.temperature
    for (a, members) in enumerate(input.members)
        (isempty(members) || isempty(input.assemblies[a])) && continue
        bore = input.label[first(input.assemblies[a])]
        isempty(input.interior_ranges[a]) && throw(ArgumentError(
            "no region lies around the members of $bore: the coaxial formulation takes the medium of the gap between the members and that conductor from it"))
        N == 0 && continue
        labels = [input.label[last(input.assemblies[m])] for m in members]
        centers = [nominal.(input.positions[m] .- input.positions[a]) for m in members]
        validate(formulation.internal_impedance, impedance!)
        filler = input.dielectric_material[first(input.interior_ranges[a])]
        κ_f = [nominal(admittivity(formulation, filler, f, θ)) for f in problem.frequencies]
        distinct = [nominal(medium.mu_r) != nominal(filler.mu_r) ||
                    any(nominal(admittivity(formulation, medium, f, θ)) != κ
                    for (f, κ) in zip(problem.frequencies, κ_f))
                    for medium in input.dielectric_material]
        radii = [nominal(maximum(
                     (input.dielectric_r_ex[l]
                     for l in input.dielectric_ranges[last(input.assemblies[m])]
                     if distinct[l]);
                     init = input.r_ex[last(input.assemblies[m])])) for m in members]
        c = nominal(minimum(
            (input.dielectric_r_in[l]
            for l in input.interior_ranges[a][2:end] if distinct[l]);
            init = input.r_in[first(input.assemblies[a])]))
        for j in eachindex(members)
            for i in 1:(j - 1)
                d = hypot((centers[i] .- centers[j])...)
                d < radii[i] + radii[j] && !isapprox(d, radii[i] + radii[j]) &&
                    throw(ArgumentError(
                        "the body of the member $(labels[i]), of radius $(radii[i]) m, reaches the body of the member $(labels[j]), of radius $(radii[j]) m, at $d m between their axes in the interior of $bore: the expansions about their axes diverge above order 0; order = 0 computes them"))
            end
            reach = hypot(centers[j]...) + radii[j]
            reach > c && !isapprox(reach, c) &&
                throw(ArgumentError(
                    "the body of the member $(labels[j]), of radius $(radii[j]) m at $(hypot(centers[j]...)) m from the axis, reaches the body of its bore $bore, of inner radius $c m: the expansions about the axes diverge above order 0; order = 0 computes them"))
        end
        for m in members
            isempty(input.members[m]) ||
                @warn("the coaxial formulation computes the impedances above order 0 of the interior of $(input.label[first(input.assemblies[m])]) as shielded by that conductor from the interior around it: the field that crosses that conductor is left out")
        end
    end
    return input
end

# The quantities of a problem and of a formulation's Γ: the hull of each one's nominal values,
# empty where the problem has none.
function domain(problem::LineParametersProblem{T},
        options::FormulationOptions{<:Union{CoaxialFormulation, FEM}}) where {T <: Real}
    N = typeof(nominal(zero(T)))
    none = (typemax(N), typemin(N))
    (; system, earth, frequencies, temperature) = problem
    deepest = last(earth.layers)
    buried = findall(placement -> nominal(placement.y) < 0, system.placements)
    reach = extrema(
        (i == j ? nominal(support(system.designs[i].outer)) :
         nominal(hypot(system.placements[i].x - system.placements[j].x,
             system.placements[i].y - system.placements[j].y))
        for i in buried, j in buried); init = none)
    m = map(f -> nominal(sqrt(2π * f * μ₀ * deepest.mu_r / deepest.rho)),
        (first(frequencies), last(frequencies)))
    return (
        frequencies = ClosedInterval(nominal(first(frequencies)), nominal(last(frequencies))),
        mu_r_earth = ClosedInterval(extrema(layer -> nominal(layer.mu_r), Base.tail(earth.layers))...),
        Γ = ClosedInterval(extrema(z -> abs(nominal(z)), options.Γ)...),
        mD = isempty(buried) ? ClosedInterval(none...) :
             ClosedInterval(first(m) * first(reach), last(m) * last(reach)),
        Δθ = ClosedInterval(extrema(N[nominal(abs(temperature - region.source.material.T0))
                                      for region in system.regions
                                      if region.source.material isa
                                         Material{<:Real, :conductor}])...))
end

# The quantities of every formulation, and the highest order that the recipes keep.
function domain(problem::LineParametersProblem,
        options::FormulationOptions{CoaxialFormulation})
    presented = invoke(domain,
        Tuple{LineParametersProblem, FormulationOptions{<:Union{CoaxialFormulation, FEM}}},
        problem, options)
    return merge(presented, (order = ClosedInterval(options.order, options.order),))
end

# A formula without a declaration holds for every value.
domain(::AbstractFormulation) = (;)

"""
$(TYPEDSIGNATURES)

Warn once for each quantity of a computation outside the domain of a formula of the
formulation.

# Arguments

- `presented`: the closed interval of each quantity that the computation presents, by its key.
- `formulation`: coaxial or finite-element formulation, the formulas of its slots.

# Returns

- `presented`, unchanged.
"""
function validate(presented::NamedTuple, formulation::Union{CoaxialFormulation, FEM})
    map(formulation.formulas) do formula
        map(keys(domain(formula)), values(domain(formula))) do key, bound
            issubset(getfield(presented, key), bound) ||
                @warn("the formula :$(nameof(formula)) holds for $key in $bound; " *
                      "the computation has $key in $(getfield(presented, key))")
        end
    end
    return presented
end

# A formula's identifier, its type's parameter, in every family.
function Base.nameof(::Union{
        InternalImpedance.Formula{ID}, InsulationImpedance.Formula{ID},
        InsulationAdmittance.Formula{ID}, SemiconAdmittance.Formula{ID},
        EarthReturn.Formula{ID}, Shunt.Formula{ID}, Resistivity.Formula{ID},
        EarthProperties.Formula{ID}}) where {ID}
    return ID
end

"""
    @air(selection, rest...)

Select the formula of the air domain of an earth slot, in a line of [`@formulation`](@ref). The
tags of a chain merge into one record with a key per domain, in any order. A domain not given
takes the default of the slot.

# Arguments

- `selection`: what the slot accepts, `:id`, `:id => (key = value, ...)` or a formula.
- `rest`: the tags after this one.

# Returns

- A `NamedTuple` with the key `air` and the keys of the tags after it.
"""
macro air(selection, rest...)
    return :(merge((air = $(esc(selection)),), $(map(esc, rest)...)))
end

"""
    @earth(selection, rest...)

Select the formula of the earth domain of an earth slot, in a chain of tags as [`@air`](@ref)
describes.

# Arguments

- `selection`: what the slot accepts, `:id`, `:id => (key = value, ...)` or a formula.
- `rest`: the tags after this one.

# Returns

- A `NamedTuple` with the key `earth` and the keys of the tags after it.
"""
macro earth(selection, rest...)
    return :(merge((earth = $(esc(selection)),), $(map(esc, rest)...)))
end

"""
    @mixed(selection, rest...)

Select the formula of the mixed domain of an earth slot, in a chain of tags as [`@air`](@ref)
describes.

# Arguments

- `selection`: what the slot accepts, `:id`, `:id => (key = value, ...)` or a formula.
- `rest`: the tags after this one.

# Returns

- A `NamedTuple` with the key `mixed` and the keys of the tags after it.
"""
macro mixed(selection, rest...)
    return :(merge((mixed = $(esc(selection)),), $(map(esc, rest)...)))
end

"""
    @formulation Type begin
        lines...
    end

Build a formulation from the lines of a block. Each line is a keyword of the constructor of
`Type`, or of [`CoaxialFormulation`](@ref) without a type. A selection per domain is a chain of
the tags [`@air`](@ref), [`@earth`](@ref) and [`@mixed`](@ref).

# Arguments

- `Type`: formulation type.
- `block`: one line per keyword, a slot or an option of the formulation.

# Returns

- The formulation that the constructor builds from the keywords.

# Errors

- The constructor throws `ArgumentError` for a slot or option that the formulation lacks.
- Julia refuses a line that is not an assignment when it lowers the call.

# Examples

```jldoctest
formulation = @formulation begin
    earth_impedance = @air :gary1976 @earth :saad1996 @mixed :lucca1994
    earth_admittance = :ideal
    transposition = true
end
@assert formulation === CoaxialFormulation(;
    earth_impedance = (air = :gary1976, earth = :saad1996, mixed = :lucca1994),
    earth_admittance = :ideal, transposition = true)
# output
```
"""
macro formulation(block)
    lines = [line for line in block.args if !(line isa LineNumberNode)]
    return Expr(:call, CoaxialFormulation, Expr(:parameters, map(esc, lines)...))
end
macro formulation(type, block)
    lines = [line for line in block.args if !(line isa LineNumberNode)]
    return Expr(:call, esc(type), Expr(:parameters, map(esc, lines)...))
end
