"""
    PSCADFormulation

Store shared formula selections, requested definitions, and physical options
for PSCAD. Backend execution settings remain computation options.
"""
struct PSCADFormulation{M <: NamedTuple, O <: NamedTuple, D <: NamedTuple} <:
       AbstractFormulation
    methods::M
    options::O
    definitions::D
end

function formulation_options(::Type{PSCADFormulation}, options::NamedTuple)::FormulationOptions
    base_frequency = get(options, :base_frequency, 50.0)
    isfinite(base_frequency) && base_frequency >= 0.1 || throw(DomainError(
        base_frequency, "PSCAD base frequency must be finite and at least 0.1 Hz"))
    physical = (;
        (key => value for (key, value) in pairs(options) if key !== :base_frequency)...)
    normalized = formulation_options(LineParametersFormulation,
        merge((reduce_bundle = false, kron_reduction = false, ideal_transposition = false), physical))
    any((normalized.reduce_bundle, normalized.kron_reduction,
        normalized.ideal_transposition)) &&
        throw(ArgumentError("PSCAD currently requires unreduced, untransposed terminal matrices"))
    return merge(normalized, (; base_frequency))
end

function _pscad_formulation(internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence, options::NamedTuple)
    selections = (; internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence)
    normalized = formulation_options(PSCADFormulation, options)
    physical = Formulation(; selections...,
        options = (;
            (key => value
        for (key, value) in pairs(normalized) if key !== :base_frequency)...))
    return PSCADFormulation(physical.methods, normalized, physical.definitions)
end

"""
    Formulation(:pscad; kwargs...)

Select PSCAD through the shared formula grammar. All formula slots and options
accept Grid inputs with product or zip composition. Earth-impedance
`:default` retains its identity and selects the native direct numerical
integration setting for overhead or underground placement, or the native Lucca
setting for mixed placement. Dielectric
`:default` is explicitly lossless. An explicit Ametani selection is represented
by the equivalent capacitance and loss tangent at the export reference
frequency; PSCAD's native frequency law and loss-tangent cap of ten still apply.
Native defaults have their own indexed implementations. External selections
accept the same `(air, earth, mixed)` shorthand as LCM. Every required native
field is compiled from those cases, and conflicting demands on shared native
Z/Y controls fail before export. Custom analytical integration settings are
not native PSCAD controls.
"""
function Formulation(::Val{:pscad};
        internal_impedance = formula(:default),
        insulation_impedance = formula(:default),
        earth_impedance = formula(:default),
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_admittance = formula(:default),
        earth_properties = formula(:default),
        pipe_impedance = formula(:default),
        temperature_dependence = formula(:default),
        options = (;), combine::Symbol = :product)
    selections = (internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence)
    return parameterize(PSCADFormulation, _pscad_formulation, (selections..., options); combine)
end

function Formulation(::Val{:pscad}, problem::LineParametersProblem, requested::PSCADFormulation)
    pscad_setting(requested, problem)
    return requested
end

function Formulation(::Val{:pscad}, ::Val{:default}, ::PipeImpedance.Formula{:default}, ::Val{:coaxial})
    nothing
end
function Formulation(::Val{:pscad}, ::Val{:default}, ::PipeImpedance.Formula{:default}, ::Val{:pipe})
    throw(ArgumentError("PSCAD Cable_Coax does not support an eccentric or multicore metallic pipe enclosure"))
end

# FormulaMethod retains exactly the same tag/kind/s/t arguments. Its backend
# payload requests native settings instead of an analytical coefficient.
function earth_impedance(
        ::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad}) where {ID, Kind, S, T}
    throw(ArgumentError("PSCAD earth_impedance :$ID, kind :$Kind: formula not implemented for source in layer $S and target in layer $T"))
end
function earth_potential_coefficient(
        ::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad}) where {ID, Kind, S, T}
    throw(ArgumentError("PSCAD earth_potential_coefficient :$ID, kind :$Kind: formula not implemented for source in layer $S and target in layer $T"))
end
function internal_impedance(::Val{ID}, ::Val{Kind}, ::Val{:pscad}) where {ID, Kind}
    throw(ArgumentError("PSCAD internal_impedance :$ID: formula not implemented for kind :$Kind"))
end
function insulation_impedance(::Val{ID}, ::Val{:pscad}) where {ID}
    throw(ArgumentError("PSCAD insulation_impedance :$ID: formula not implemented"))
end

function earth_impedance(::Val{:default}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::Val{:default}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::Val{:default}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(::Val{:default}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(::Val{:Gary1976}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 0, readback = "DERISEMLYEN"),)
end
function earth_impedance(::Val{:Carson1926}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::Val{:Pollaczek1926}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::Val{:WedepohlWilcox1973}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 0, readback = "WEDEPOHL"),)
end
function earth_impedance(::Val{:Saad1996}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 3, readback = "SAAD"),)
end
function earth_impedance(
        ::Val{:Ametani2009}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 0, readback = "AMETANIL"),)
end
function earth_impedance(
        ::Val{:Ametani2009}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 0, readback = "AMETANIL"),)
end
function earth_impedance(
        ::Val{:Lucca1994}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(
        ::Val{:Lucca1994}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end

# The native default potential follows the selected ground solver. There is no
# independent switch for a Julia potential-coefficient equation.
function earth_potential_coefficient(::Val{:default}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(::Val{:default}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(
        ::Val{:default}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(
        ::Val{:default}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (;)
end
function internal_impedance(
        ::Val{:default}, ::Union{
            Val{:inner}, Val{:outer}, Val{:mutual}}, ::Val{:pscad})
    (;) # Fixed Cable_Coax conductor calculation; no native author/formula switch.
end
insulation_impedance(::Val{:default}, ::Val{:pscad}) = (;)

"""
    formulas(owner, kind, source, target)

List native external equations using the shared indexed method declarations.
`owner` is `EarthImpedance` or `EarthAdmittance`; selectors are `Val` values.
This lists individual cases, not complete models or combined author identities.
"""
function formulas(owner::Module, kind::Val, source::Val, target::Val)
    equation = owner === EarthImpedance ? earth_impedance :
               owner === EarthAdmittance ? earth_potential_coefficient :
               throw(ArgumentError(
        "indexed PSCAD formula discovery requires EarthImpedance or EarthAdmittance"))
    fallback = which(equation, Tuple{Val, Val, Val, Val, Val{:pscad}})
    return Tuple(id
    for id in owner.formulas()
    if
    which(equation, Tuple{
        Val{id}, typeof(kind), typeof(source), typeof(target), Val{:pscad}}) !== fallback)
end

"Compile every required native setting from the physical indexed interactions."
function pscad_setting(formulation::PSCADFormulation, problem::LineParametersProblem)
    validate(problem)
    model = problem.earth_props
    validate(model, Val(:pscad))
    problem.Γ === nothing || all(iszero, problem.Γ) ||
        throw(ArgumentError(
            "PSCAD does not accept an explicit nonzero longitudinal propagation constant"))
    relation = formulation.methods.earth_properties
    (relation === nothing ||
     relation === LineCableModels.Earth.FrequencyDependent.Formula(:default)) ||
        throw(ArgumentError("PSCAD does not implement the selected frequency-dependent soil relation"))
    for name in
        (:internal_impedance, :insulation_impedance, :earth_impedance, :earth_admittance)
        selection = getproperty(formulation.methods, name)
        leaves = selection isa NamedTuple ? values(selection) : (selection,)
        for selected in leaves
            isempty(selected.hooks) && isempty(selected.parameters) &&
            all(isempty, values(selected.options)) || throw(ArgumentError(
                "PSCAD cannot evaluate an analytical $name override or integration selection"))
            if name in (:earth_impedance, :earth_admittance)
                selected.equivalent_earth === nothing || throw(ArgumentError(
                    "PSCAD does not execute $name equivalent-earth reductions"))
            end
        end
    end
    for name in (:insulation_admittance, :semicon_admittance)
        selected = getproperty(formulation.methods, name)
        formula_id(selected) in (:default, :Ametani2004) || throw(ArgumentError(
            "PSCAD does not implement $name :$(formula_id(selected))"))
    end
    for design in problem.system.designs
        Formulation(Val(:pscad), formulation.methods.pipe_impedance, design)
    end
    T = eltype(problem)
    blueprints = Engine.CableBlueprint{T}[Engine.flatten(LineCableModelsCoaxial(), design, T)
                                          for design in problem.system.designs]
    input = Engine.lineinput(problem, blueprints)
    pairs = Engine.earth_pairs(first.(input.cable.assemblies), input.horz, input.vert,
        input.horz_sep, problem.earth_props)
    settings = Dict{Symbol, NamedTuple}()
    interactions = map(formulation.methods[(:earth_impedance, :earth_admittance)]) do selection
        selection isa NamedTuple && validate(selection, problem.earth_props)
        records = NamedTuple{
            (:formula, :kind, :source, :target), Tuple{Symbol, Symbol, Int, Int}}[]
        for pair in pairs
            validate(pair)
            selected = Formulation(selection, Val.(pair.layers)...)
            binding = FormulaMethod(selected, pair)
            # Explicit author choices must also exist in their analytical owner.
            formula_id(selected) === :default || validate(selected, pair)
            record = (formula = formula_id(selected),
                kind = pair.row == pair.column ? :self : :mutual,
                source = pair.layers[1], target = pair.layers[2])
            record in records && continue
            push!(records, record)
            for (field, setting) in Base.pairs(binding(Val(:pscad)))
                haskey(settings, field) && settings[field] != setting &&
                    throw(ArgumentError(
                        "PSCAD shared field $field has conflicting formula selections; requested Z/Y cannot be selected independently"))
                settings[field] = setting
            end
        end
        records
    end
    kinds = any(indices -> length(indices) > 1, input.cable.assemblies) ?
            (:inner, :outer, :mutual) : (:outer,)
    for kind in kinds
        FormulaMethod(Val(formula_id(formulation.methods.internal_impedance)),
            internal_impedance, Val(kind))(Val(:pscad))
    end
    FormulaMethod(Val(formula_id(formulation.methods.insulation_impedance)), insulation_impedance)(Val(:pscad))
    # Unused native slots are set deterministically and retained too; they do not
    # authorize any additional physical case.
    ground = (
        EarthForm2 = get(settings, :EarthForm2, (
            value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION")),
        EarthForm = get(settings, :EarthForm, (
            value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION")),
        EarthForm3 = get(settings, :EarthForm3, (value = 2, readback = "LUCCA")))
    return (ground = ground,
        frequency = (enablf = (value = 1, readback = "YES"),
            FS = (value = Float64(first(problem.frequencies)),
                readback = Float64(first(problem.frequencies))),
            FE = (value = Float64(last(problem.frequencies)),
                readback = Float64(last(problem.frequencies))),
            Numf = (value = length(problem.frequencies)-1,
                readback = length(problem.frequencies)-1)),
        configuration = (Freq = (value = Float64(formulation.options.base_frequency),
            readback = Float64(formulation.options.base_frequency)),), interactions = interactions)
end

function Base.NamedTuple(formulation::PSCADFormulation)
    identifier = function (definition)
        definition === nothing && return nothing
        definition isa NamedTuple && return map(identifier, definition)
        definition isa Symbol ? definition : formula_id(definition)
    end
    fields = keys(formulation.definitions)
    HomogeneousIDs = NamedTuple{(:air, :earth, :mixed), NTuple{3, Symbol}}
    Identifiers = NamedTuple{fields,
        NTuple{length(fields), Union{Nothing, Symbol, HomogeneousIDs}}}
    requested = Identifiers(map(identifier, formulation.definitions))
    methods = formulation.methods
    return (
        schema_version = 3,
        backend = :pscad,
        type = string(parentmodule(typeof(formulation)), ".", nameof(typeof(formulation))),
        requested,
        raw = Dict{Symbol, Any}(:selections => formulation.definitions),
        effective = Identifiers(merge(map(identifier, methods), (pipe_impedance = nothing,))),
        assumptions = (
            internal_impedance = "Fixed PSCAD Cable_Coax conductor approximation; backend-owned default, no source alias or native Bessel switch",
            insulation_impedance = "PSCAD native Cable_Coax magnetic calculation",
            earth_admittance = "Native indexed potential coefficients; coupled to the selected earth solver",
            insulation_admittance = description(methods.insulation_admittance),
            semicon_admittance = description(methods.semicon_admittance),
            dielectric_equivalence = "Reference-frequency equivalent capacitance and loss tangent; native PSCAD frequency law; loss tangent capped at 10",
            earth = "One homogeneous earth layer; no FrequencyDependent relation",
            pipe_impedance = "Cable_Coax only; shared eccentric metallic enclosure unsupported"),
        options = formulation.options)
end

function computation_details(::Type{<:PSCADFormulation}, result::LineParameters)::ComputationDetails
    return details(result)
end
