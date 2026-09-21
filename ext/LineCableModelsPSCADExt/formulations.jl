import LineCableModels: formula_id

"""
$(TYPEDEF)

Select one PSCAD-owned equation or fixed component calculation. `F` identifies
the physical family; `ID` identifies its native setting. This selection carries
no analytical equation, callback, or numerical workspace.
"""
struct NativeFormula{F, ID} <: AbstractFormulation end

NativeFormula{F}(selected::NativeFormula{F}) where {F} = selected
NativeFormula{F}(identifier::Symbol) where {F} = NativeFormula{F}(Val(identifier))
NativeFormula{InternalImpedance.Formula}(::Val{:default}) =
    NativeFormula{InternalImpedance.Formula, :cable_coax}()
NativeFormula{InsulationImpedance.Formula}(::Val{:default}) =
    NativeFormula{InsulationImpedance.Formula, :cable_coax}()
NativeFormula{EarthImpedance.Formula}(::Val{:default}) =
    NativeFormula{EarthImpedance.Formula, :direct_lucca}()
NativeFormula{EarthAdmittance.Formula}(::Val{:default}) =
    NativeFormula{EarthAdmittance.Formula, :coupled}()
NativeFormula{InternalImpedance.Formula}(::Val{:cable_coax}) =
    NativeFormula{InternalImpedance.Formula, :cable_coax}()
NativeFormula{InsulationImpedance.Formula}(::Val{:cable_coax}) =
    NativeFormula{InsulationImpedance.Formula, :cable_coax}()
NativeFormula{EarthImpedance.Formula}(::Val{:direct_lucca}) =
    NativeFormula{EarthImpedance.Formula, :direct_lucca}()
NativeFormula{EarthAdmittance.Formula}(::Val{:coupled}) =
    NativeFormula{EarthAdmittance.Formula, :coupled}()
function NativeFormula{F}(::Val{ID}) where {F, ID}
    F(ID) # Admit a named scientific selection; native case support is checked at preflight.
    return NativeFormula{F, ID}()
end
function NativeFormula{F}(selection::LineCableModels.FormulaDefinition{ID, Order}) where {F, ID, Order}
    Order === :default && selection.equivalent_earth === nothing || throw(ArgumentError(
        "PSCAD native equations do not execute equivalent-earth reductions"))
    isempty(selection.parameters) && isempty(selection.options.data) ||
        throw(ArgumentError("PSCAD native equations do not accept analytical parameters or numerical controls"))
    return NativeFormula{F}(Val(ID))
end

function NativeFormula{F}(selected::F) where {F}
    empty_controls = selected isa InternalImpedance.Formula ?
        all(isempty, values(selected.options.data)) : isempty(selected.options.data)
    isempty(selected.parameters) && empty_controls || throw(ArgumentError(
        "PSCAD native equations do not accept analytical parameters or numerical controls"))
    if selected isa Union{EarthImpedance.Formula,EarthAdmittance.Formula}
        selected.equivalent_earth === nothing || throw(ArgumentError(
            "PSCAD native equations do not execute equivalent-earth reductions"))
    end
    return NativeFormula{F}(formula_id(selected))
end

formula_id(::NativeFormula{F, ID}) where {F, ID} = ID
formula_id(::Type{<:NativeFormula{F, ID}}) where {F, ID} = ID
Base.NamedTuple(selected::NativeFormula) =
    (identifier=formula_id(selected), parameters=(;), options=(;))
Base.pairs(::Type{<:NativeFormula{F}}; quantity=nothing) where {F} =
    pairs(map(_ -> NativeFormula{F}, (; pairs(F; quantity)...)))
formulation_options(::NativeFormula) = FormulationOptions()
description(selected::NativeFormula; compact::Bool=false) = description(typeof(selected); compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :gary1976}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:gary1976}; compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :carson1926}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:carson1926}; compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :pollaczek1926}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:pollaczek1926}; compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :wedepohl1973}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:wedepohl1973}; compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :saad1996}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:saad1996}; compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :ametani2009}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:ametani2009}; compact)
description(::Type{NativeFormula{EarthImpedance.Formula, :lucca1994}}; compact::Bool=false) =
    "PSCAD " * description(EarthImpedance.Formula{:lucca1994}; compact)
description(::Type{<:NativeFormula{InternalImpedance.Formula, :cable_coax}}; compact::Bool=false) =
    compact ? "Cable_Coax" : "PSCAD Cable_Coax conductor calculation"
description(::Type{<:NativeFormula{InsulationImpedance.Formula, :cable_coax}}; compact::Bool=false) =
    compact ? "Cable_Coax" : "PSCAD Cable_Coax magnetic insulation calculation"
description(::Type{<:NativeFormula{EarthImpedance.Formula, :direct_lucca}}; compact::Bool=false) =
    compact ? "Direct/Lucca" : "PSCAD direct earth integration with Lucca mixed interactions"
description(::Type{<:NativeFormula{EarthAdmittance.Formula, :coupled}}; compact::Bool=false) =
    compact ? "Coupled potential" : "PSCAD potential coefficients coupled to the native ground solver"

Formulation(selected::NativeFormula, ::Val{S}, ::Val{T}) where {S, T} = selected
function validate(selected::NativeFormula{F, ID}, pair::Engine.EarthPair) where {F, ID}
    validate(F(ID), pair)
    return selected
end
function validate(selected::Union{
        NativeFormula{EarthImpedance.Formula, :direct_lucca},
        NativeFormula{EarthAdmittance.Formula, :coupled}}, pair::Engine.EarthPair)
    validate(pair)
    return selected
end
function FormulaMethod(selected::NativeFormula{EarthImpedance.Formula}, pair::Engine.EarthPair)
    FormulaMethod(selected, earth_impedance,
        Val(pair.row == pair.column ? :self : :mutual), Val.(pair.layers)...)
end
function FormulaMethod(selected::NativeFormula{EarthAdmittance.Formula}, pair::Engine.EarthPair)
    FormulaMethod(selected, earth_potential_coefficient,
        Val(pair.row == pair.column ? :self : :mutual), Val.(pair.layers)...)
end

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
    native = (:internal_impedance, :insulation_impedance, :earth_impedance, :earth_admittance)
    return pairs((; (name => (name in native ? NativeFormula{owner} : owner)
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

function _pscad_formulation(internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence, options::FormulationOptions)
    selections = (; internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        pipe_impedance, temperature_dependence)
    normalized = formulation_options(PSCADFormulation, options)
    methods = (; (name => (selections[name] === nothing && name in (:earth_properties,:temperature_dependence) ? nothing :
        Formulation(owner, selections[name])) for (name, owner) in pairs(PSCADFormulation))...)
    definitions = map(selections, methods) do requested, selected
        requested isa NamedTuple ? NamedTuple{keys(selected)}(requested) : requested
    end
    return PSCADFormulation(methods, normalized, definitions)
end

"""
    Formulation(:pscad; kwargs...)

Select PSCAD through the shared formula grammar. All formula slots and options
accept Grid inputs with product or zip composition. Earth-impedance
`:default` resolves to the native `:direct_lucca` selection: direct numerical
integration for overhead or underground placement and Lucca for mixed placement.
Internal and insulation impedance resolve to native `:cable_coax`, and earth
potential to `:coupled`. These are not analytical formula aliases. Dielectric
`:default` routes to `:lossless`. An explicit `:lossy` selection is represented
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

# FormulaMethod retains the selected native type and kind/s/t arguments. Its backend
# payload requests native settings instead of an analytical coefficient.
function earth_impedance(
        ::NativeFormula{EarthImpedance.Formula, ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad}) where {ID, Kind, S, T}
    throw(ArgumentError("PSCAD earth_impedance :$ID, kind :$Kind: formula not implemented for source in layer $S and target in layer $T"))
end
function earth_potential_coefficient(
        ::NativeFormula{EarthAdmittance.Formula, ID}, ::Val{Kind}, ::Val{S}, ::Val{T}, ::Val{:pscad}) where {ID, Kind, S, T}
    throw(ArgumentError("PSCAD earth_potential_coefficient :$ID, kind :$Kind: formula not implemented for source in layer $S and target in layer $T"))
end
function internal_impedance(::NativeFormula{InternalImpedance.Formula, ID}, ::Val{Kind}, ::Val{:pscad}) where {ID, Kind}
    throw(ArgumentError("PSCAD internal_impedance :$ID: formula not implemented for kind :$Kind"))
end
function insulation_impedance(::NativeFormula{InsulationImpedance.Formula, ID}, ::Val{:pscad}) where {ID}
    throw(ArgumentError("PSCAD insulation_impedance :$ID: formula not implemented"))
end

function earth_impedance(::NativeFormula{EarthImpedance.Formula, :direct_lucca}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :direct_lucca}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :direct_lucca}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :direct_lucca}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :gary1976}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 0, readback = "DERISEMLYEN"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :carson1926}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (EarthForm2 = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :pollaczek1926}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 2, readback = "DIRECT_NUMERICAL_INTEGRATION"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :wedepohl1973}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 0, readback = "WEDEPOHL"),)
end
function earth_impedance(::NativeFormula{EarthImpedance.Formula, :saad1996}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (EarthForm = (value = 3, readback = "SAAD"),)
end
function earth_impedance(
        ::NativeFormula{EarthImpedance.Formula, :ametani2009}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 0, readback = "AMETANIL"),)
end
function earth_impedance(
        ::NativeFormula{EarthImpedance.Formula, :ametani2009}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 0, readback = "AMETANIL"),)
end
function earth_impedance(
        ::NativeFormula{EarthImpedance.Formula, :lucca1994}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end
function earth_impedance(
        ::NativeFormula{EarthImpedance.Formula, :lucca1994}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (EarthForm3 = (value = 2, readback = "LUCCA"),)
end

# The native default potential follows the selected ground solver. There is no
# independent switch for a Julia potential-coefficient equation.
function earth_potential_coefficient(::NativeFormula{EarthAdmittance.Formula, :coupled}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(::NativeFormula{EarthAdmittance.Formula, :coupled}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(
        ::NativeFormula{EarthAdmittance.Formula, :coupled}, ::Val{:mutual}, ::Val{1}, ::Val{2}, ::Val{:pscad})
    (;)
end
function earth_potential_coefficient(
        ::NativeFormula{EarthAdmittance.Formula, :coupled}, ::Val{:mutual}, ::Val{2}, ::Val{1}, ::Val{:pscad})
    (;)
end
function internal_impedance(
        ::NativeFormula{InternalImpedance.Formula, :cable_coax}, ::Union{
            Val{:inner}, Val{:outer}, Val{:transfer}}, ::Val{:pscad})
    (;) # Fixed Cable_Coax conductor calculation; no native author/formula switch.
end
insulation_impedance(::NativeFormula{InsulationImpedance.Formula, :cable_coax}, ::Val{:pscad}) = (;)

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
    fallback = which(equation, Tuple{NativeFormula{owner.Formula}, Val, Val, Val, Val{:pscad}})
    return Tuple(id
    for id in (owner === EarthImpedance ? :direct_lucca : :coupled, owner.formulas()...)
    if
    which(equation, Tuple{
        NativeFormula{owner.Formula, id}, typeof(kind), typeof(source), typeof(target), Val{:pscad}}) !== fallback)
end

"""
Compile every required native setting from the physical indexed interactions.
"""
function pscad_setting(formulation::PSCADFormulation, problem::LineParametersProblem)
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
            validate(selected, pair)
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
            internal_impedance = "Fixed PSCAD Cable_Coax conductor approximation; backend-owned default, no source alias or native Bessel switch",
            insulation_impedance = "PSCAD native Cable_Coax magnetic calculation",
            earth_admittance = "Native indexed potential coefficients; coupled to the selected earth solver",
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
