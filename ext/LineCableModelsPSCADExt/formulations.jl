import LineCableModels: formula_id

"""
    PSCADFormulation

Store shared formula selections, requested definitions, and physical options
for PSCAD. Backend execution settings remain computation options.
"""
struct PSCADFormulation{M <: NamedTuple, O <: FormulationOptions, D <: NamedTuple} <:
       AbstractFormulation
    methods::M
    options::O
    definitions::D
end

"""Identify PSCAD without server access or native setting compilation."""
description(::Type{<:PSCADFormulation}; compact::Bool=false) = "PSCAD"
description(::PSCADFormulation; compact::Bool=false) = description(PSCADFormulation;compact)
formula_id(::Type{<:PSCADFormulation}) = :pscad
formula_id(::PSCADFormulation) = :pscad
formulation_options(value::PSCADFormulation) = value.options
function Base.pairs(::Type{PSCADFormulation}; quantity=nothing)
    return pairs((; (name => owner
        for (name, owner) in pairs(LineParametersFormulation; quantity)
        if name !== :shunt_model)...))
end
Base.pairs(value::PSCADFormulation;quantity=nothing) =
    pairs(PSCADFormulation,(methods=value.methods,requested=value.definitions,options=value.options.data);quantity)
description(::Type{PSCADFormulation},slot::Val;compact::Bool=false) =
    description(LineParametersFormulation,slot;compact)
description(::Type{PSCADFormulation},slot::Union{Val{:reduce_bundle},Val{:kron_reduction},Val{:ideal_transposition}},
    value::Bool;compact::Bool=false) = description(LineParametersFormulation,slot,value;compact)
description(::Type{PSCADFormulation},::Val{:base_frequency},value::Real;compact::Bool=false) =
    "base frequency="*string(value)*" Hz"
Base.pairs(::Type{PSCADFormulation},retained::NamedTuple;quantity=nothing) =
    pairs(LineParametersFormulation,retained;quantity,owner=PSCADFormulation)

function formulation_options(::Type{PSCADFormulation}, record::FormulationOptions)::FormulationOptions
    options = record.data
    base_frequency = get(options, :base_frequency, 50.0)
    _pscad_deterministic(typeof(base_frequency))
    isfinite(base_frequency) && base_frequency >= 0.1 || throw(DomainError(
        base_frequency, "PSCAD base frequency must be finite and at least 0.1 Hz"))
    physical = (;
        (key => value for (key, value) in pairs(options) if key !== :base_frequency)...)
    normalized = formulation_options(LineParametersFormulation,
        FormulationOptions(merge((reduce_bundle = false, kron_reduction = false, ideal_transposition = false), physical)))
    any((normalized.data.reduce_bundle, normalized.data.kron_reduction,
        normalized.data.ideal_transposition)) &&
        throw(ArgumentError("PSCAD currently requires unreduced, untransposed terminal matrices"))
    return FormulationOptions(merge(normalized.data, (; base_frequency)))
end

# Resolve the backend's explicit defaults before the Engine constructs selections.
# Partial recipes retain their original branches and are validated by the owner.
_pscad_default(selected, default) = selected
_pscad_default(selected::Symbol, default) = _pscad_default(LineCableModels.formula(selected), default)
function _pscad_default(selected::LineCableModels.FormulaDefinition{:default, Order}, default) where {Order}
    Order === :default && selected.equivalent_earth === nothing || throw(ArgumentError(
        "PSCAD equations do not execute equivalent-earth reductions"))
    isempty(selected.parameters) && isempty(selected.options.data) || throw(ArgumentError(
        "PSCAD equations do not accept analytical parameters or numerical controls"))
    return default
end
function _pscad_default(selected::NamedTuple, default)
    return NamedTuple{keys(selected)}(map(keys(selected)) do name
        _pscad_default(selected[name], default isa NamedTuple ? get(default, name, selected[name]) : default)
    end)
end

function validate(selected::Union{InternalImpedance.Formula, InsulationImpedance.Formula,
        EarthImpedance.Formula, EarthAdmittance.Formula}, ::Val{:pscad})
    empty_controls = selected isa InternalImpedance.Formula ?
        all(isempty, values(selected.options.data)) : isempty(selected.options.data)
    isempty(selected.parameters) && empty_controls || throw(ArgumentError(
        "PSCAD equations do not accept analytical parameters or numerical controls"))
    if selected isa Union{EarthImpedance.Formula, EarthAdmittance.Formula}
        selected.equivalent_earth === nothing || throw(ArgumentError(
            "PSCAD equations do not execute equivalent-earth reductions"))
    end
    return selected
end

function _pscad_formulation(internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence, options::FormulationOptions)
    selections = (; internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence)
    normalized = formulation_options(PSCADFormulation, options)
    defaults = (internal_impedance = :wedepohl1973, insulation_impedance = :ametani1980,
        earth_impedance = (air = :carson1926, earth = :pollaczek1926, mixed = :lucca1994),
        earth_admittance = :ideal)
    methods = (; (name => begin
        requested = selections[name]
        if requested === nothing && name in (:earth_properties, :temperature_dependence)
            nothing
        else
            selected = Formulation(owner, haskey(defaults, name) ?
                _pscad_default(requested, defaults[name]) : requested)
            if haskey(defaults, name)
                for leaf in (selected isa NamedTuple ? values(selected) : (selected,))
                    leaf === nothing || validate(leaf, Val(:pscad))
                end
            end
            selected
        end
    end for (name, owner) in pairs(PSCADFormulation))...)
    definitions = map(selections, methods) do requested, selected
        requested isa NamedTuple ? NamedTuple{keys(selected)}(requested) : requested
    end
    return PSCADFormulation(methods, normalized, definitions)
end

"""
    Formulation(:pscad; kwargs...)

Select PSCAD through the shared formula grammar. All formula slots and options
accept Grid inputs with product or zip composition. Earth-impedance
`:default` resolves to `(air=:carson1926, earth=:pollaczek1926, mixed=:lucca1994)`.
Internal impedance selects `:wedepohl1973` for inner, outer, and transfer surfaces;
magnetic insulation impedance selects `:ametani1980`. Earth potential coefficients
select `:ideal`: electrostatic images in air, zero for buried and mixed pairs.
Every scientific selection is registered by its Engine owner. Dielectric
`:default` routes to `:lossless`. An explicit `:lossy` selection is represented
by the equivalent capacitance and loss tangent at the export reference
frequency; PSCAD's native frequency law and loss-tangent cap of ten still apply.
External selections accept the same `(air, earth, mixed)` shorthand as LCM.
Every required native field is compiled from the indexed Engine selections.
PSCAD has no separate potential-model selector. Native direct integration can
produce nonzero aerial conductance under the `:ideal` selection; raw native
matrices are preserved. Unsupported equations and analytical integration settings
fail before export.
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
        options = FormulationOptions(), combine::Symbol = :product)
    selections = (internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence)
    return parameterize(PSCADFormulation, (inputs...) -> _pscad_formulation(inputs[1:end-1]...,
        last(inputs) isa NamedTuple ? FormulationOptions(last(inputs)) : last(inputs)),
        (selections..., options); combine)
end

function Formulation(::Val{:pscad}, problem::LineParametersProblem, requested::PSCADFormulation)
    pscad_setting(requested, problem)
    return requested
end

function Formulation(::Val{:pscad}, ::PipeImpedance.Formula{:none}, ::Val{:coaxial})
    nothing
end
function Formulation(::Val{:pscad}, ::PipeImpedance.Formula{:none}, ::Val{:pipe})
    throw(ArgumentError("PSCAD Cable_Coax does not support an eccentric or multicore metallic pipe enclosure"))
end

# The Engine owns each selected formula and its kind/s/t binding. These methods
# translate that binding to PSCAD settings without evaluating the Julia equation.
function earth_impedance(
        ::EarthImpedance.Formula{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad}) where {ID, Kind, S, T}
    throw(ArgumentError("PSCAD earth_impedance :$ID, kind :$Kind: formula not implemented for source in layer $S and target in layer $T"))
end
function earth_potential_coefficient(
        ::EarthAdmittance.Formula{ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad}) where {ID, Kind, S, T}
    throw(ArgumentError("PSCAD earth_potential_coefficient :$ID, kind :$Kind: formula not implemented for source in layer $S and target in layer $T"))
end
function internal_impedance(::InternalImpedance.Formula{ID}, ::Val{Kind}, ::Val{:pscad}) where {ID, Kind}
    throw(ArgumentError("PSCAD internal_impedance :$ID: formula not implemented for kind :$Kind"))
end
function insulation_impedance(::InsulationImpedance.Formula{ID}, ::Val{:pscad}) where {ID}
    throw(ArgumentError("PSCAD insulation_impedance :$ID: formula not implemented"))
end

function earth_impedance(::EarthImpedance.Formula{:gary1976}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 0, readback = "DERISEMLYEN"),)
end
function earth_impedance(::EarthImpedance.Formula{:carson1926}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::EarthImpedance.Formula{:pollaczek1926}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::EarthImpedance.Formula{:wedepohl1973}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 0, readback = "WEDEPOHL"),)
end
function earth_impedance(::EarthImpedance.Formula{:saad1996}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 3, readback = "SAAD"),)
end
function earth_impedance(
        ::EarthImpedance.Formula{:ametani2009}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 0, readback = "AMETANIL"),)
end
function earth_impedance(
        ::EarthImpedance.Formula{:ametani2009}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 0, readback = "AMETANIL"),)
end
function earth_impedance(
        ::EarthImpedance.Formula{:lucca1994}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(
        ::EarthImpedance.Formula{:lucca1994}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end

# The documented external potential law uses ideal images in air and zero
# coefficients for buried/mixed pairs. There is no independent native selector;
# direct integration can deviate from strict ideal behavior in the aerial block.
function earth_potential_coefficient(::EarthAdmittance.Formula{:ideal}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(::EarthAdmittance.Formula{:ideal}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(
        ::EarthAdmittance.Formula{:ideal}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(
        ::EarthAdmittance.Formula{:ideal}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (;)
end
function internal_impedance(
        ::InternalImpedance.Formula{:wedepohl1973}, ::Union{
            Val{:inner}, Val{:outer}, Val{:transfer}}, ::Val{:pscad})
    (;) # PSCAD fixes the conductor surfaces to the Wedepohl-Wilcox approximation.
end
insulation_impedance(::InsulationImpedance.Formula{:ametani1980}, ::Val{:pscad}) = (;)

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
    fallback = which(equation, Tuple{owner.Formula, Val, Val, Val, Val{:pscad}})
    return Tuple(id
    for id in owner.formulas()
    if
    which(equation, Tuple{
        owner.Formula{id}, typeof(kind), typeof(source), typeof(target), Val{:pscad}}) !== fallback)
end

"""
Compile every required native setting from the physical indexed interactions.
"""
function pscad_setting(formulation::PSCADFormulation, problem::LineParametersProblem)
    _pscad_deterministic(eltype(problem), typeof(formulation.options.data.base_frequency))
    return pscad_setting(formulation, problem, _pscad_blueprints(problem.system))
end

function pscad_setting(formulation::PSCADFormulation, problem::LineParametersProblem, blueprints)
    validate(problem)
    model = problem.earth_props
    validate(model, Val(:pscad))
    relation = formulation.methods.earth_properties
    (relation === nothing ||
     relation === LineCableModels.Earth.FrequencyDependent.Formula(:default)) ||
        throw(ArgumentError("PSCAD does not implement the selected frequency-dependent soil relation"))
    for name in (:insulation_admittance, :semicon_admittance)
        selected = getproperty(formulation.methods, name)
        selected isa Union{InsulationAdmittance.Formula{:lossless},InsulationAdmittance.Formula{:lossy},
            SemiconAdmittance.Formula{:lossless},SemiconAdmittance.Formula{:lossy}} || throw(ArgumentError(
            "PSCAD does not implement $name :$(formula_id(selected))"))
    end
    for design in problem.system.designs
        Formulation(Val(:pscad), formulation.methods.pipe_impedance, design)
    end
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
            # Registration and physical validity belong to the equation owner;
            # execution availability is selected by the native binding.
            validate(selected, pair)
            record = (formula = formula_id(selected),
                kind = pair.row == pair.column ? :self : :mutual,
                source = pair.layers[1], target = pair.layers[2])
            record in records && continue
            push!(records, record)
            for (field, setting) in Base.pairs(binding(Val(:pscad)))
                haskey(settings, field) && settings[field] != setting &&
                    throw(ArgumentError(
                        "PSCAD field $field has conflicting formula selections"))
                settings[field] = setting
            end
        end
        records
    end
    kinds = any(indices -> length(indices) > 1, input.cable.assemblies) ?
            (:inner, :outer, :transfer) : (:outer,)
    for kind in kinds
        selection = formulation.methods.internal_impedance
        selected = selection isa NamedTuple ? selection[kind] : selection
        FormulaMethod(selected,
            internal_impedance, Val(kind))(Val(:pscad))
    end
    FormulaMethod(formulation.methods.insulation_impedance, insulation_impedance)(Val(:pscad))
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
        configuration = (Freq = (value = Float64(formulation.options.data.base_frequency),
            readback = Float64(formulation.options.data.base_frequency)),), interactions = interactions)
end

"""Record consumed PSCAD identifiers and fixed native assumptions with bounded field types."""
function computation_details(formulation::PSCADFormulation)::ComputationDetails
    methods = formulation.methods
    return ComputationDetails(merge((
        schema_version = 4,
        assumptions = (
            internal_impedance = "Wedepohl-Wilcox (1973) inner, outer, and transfer conductor surface approximations",
            insulation_impedance = "Ametani (1980) concentric-insulation magnetic impedance",
            earth_admittance = "Ideal-earth images in air; zero external potential coefficients for buried and mixed pairs. Native direct integration can add aerial conductance; native matrices are preserved",
            insulation_admittance = description(methods.insulation_admittance),
            semicon_admittance = description(methods.semicon_admittance),
            dielectric_equivalence = "Reference-frequency equivalent capacitance and loss tangent; native PSCAD frequency law; loss tangent capped at 10",
            earth = "One homogeneous earth layer; no FrequencyDependent relation",
            pipe_impedance = "Cable_Coax only; shared eccentric metallic enclosure unsupported"),
        ),NamedTuple(formulation)))
end

function computation_details(::Type{<:PSCADFormulation}, result::LineParameters)::ComputationDetails
    return details(result)
end

"""Expose complete requested PSCAD formula choices and native configuration options."""
function Base.NamedTuple(value::PSCADFormulation)
    record = function (selected)
        selected === nothing && return nothing
        selected isa Symbol && return NamedTuple(formula(selected))
        selected isa NamedTuple && return map(record,selected)
        return NamedTuple(selected)
    end
    Record=NamedTuple{(:backend,:requested,:methods,:options),
        Tuple{Symbol,NamedTuple,NamedTuple,NamedTuple}}
    return Record((:pscad,map(record,value.definitions),map(record,value.methods),value.options.data))
end
