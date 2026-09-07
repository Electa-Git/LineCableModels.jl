"""
    PSCADFormulation

Store shared formula selections, requested definitions, and physical options
for PSCAD. Backend execution settings remain computation options.
"""
struct PSCADFormulation{M <: NamedTuple, O <: NamedTuple, D <: NamedTuple} <: AbstractFormulation
    methods::M
    options::O
    definitions::D
end

function formulation_options(::Type{PSCADFormulation}, options::NamedTuple)::FormulationOptions
    normalized = formulation_options(LineParametersFormulation,
        merge((reduce_bundle=false, kron_reduction=false, ideal_transposition=false), options))
    any((normalized.reduce_bundle, normalized.kron_reduction, normalized.ideal_transposition)) &&
        throw(ArgumentError("PSCAD Gauntlet currently requires unreduced, untransposed terminal matrices"))
    return normalized
end

function _pscad_formulation(internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        equivalent_earth, pipe_impedance, options::NamedTuple)
    selections = (; internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        equivalent_earth, pipe_impedance)
    physical = Formulation(; selections..., options=formulation_options(PSCADFormulation, options))
    return PSCADFormulation(physical.methods, physical.options, physical.definitions)
end

"""
    Formulation(:pscad; kwargs...)

Select PSCAD through the shared formula grammar. All formula slots and options
accept Grid inputs with product or zip composition. Earth-impedance
`:default` resolves from placement to direct numerical integration (overhead
or underground), or PSCAD's Lucca option for mixed placement. Dielectric
`:default` is explicitly lossless. An explicit Ametani selection is represented
by the equivalent capacitance and loss tangent at the export reference
frequency; PSCAD's native frequency law and loss-tangent cap of ten still apply.
Other fixed PSCAD calculations are recorded as backend assumptions, not as
implementations of the requested analytical kernels.
"""
function Formulation(::Val{:pscad};
        internal_impedance=formula(:default),
        insulation_impedance=formula(:default),
        earth_impedance=formula(:default),
        insulation_admittance=formula(:default),
        semicon_admittance=formula(:default),
        earth_admittance=formula(:default),
        earth_properties=formula(:default),
        equivalent_earth=formula(:default),
        pipe_impedance=formula(:default),
        options=(;), combine::Symbol=:product)
    selections = (internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance, earth_properties,
        equivalent_earth, pipe_impedance)
    return parameterize(PSCADFormulation, _pscad_formulation, (selections..., options); combine)
end

function Formulation(::Val{:pscad}, problem::LineParametersProblem, requested::PSCADFormulation)
    heights = getproperty.(problem.system.positions, :y)
    placement = all(>(0), heights) ? Val(:overhead) :
                all(<(0), heights) ? Val(:underground) : Val(:mixed)
    !problem.earth_props.vertical_layers && length(problem.earth_props.layers) == 2 ||
        throw(ArgumentError("PSCAD Cable_Coax Gauntlet supports air and one homogeneous horizontal earth half-space"))
    problem.Γ === nothing || all(iszero, problem.Γ) || throw(ArgumentError(
        "PSCAD does not accept an explicit nonzero longitudinal propagation constant"))
    requested.methods.earth_properties === nothing || throw(ArgumentError(
        "PSCAD does not implement the selected frequency-dependent soil relation"))
    for name in (:insulation_admittance, :semicon_admittance)
        selected = getproperty(requested.methods, name)
        formula_id(selected) in (:default, :Ametani2004) || throw(ArgumentError(
            "PSCAD does not implement $name :$(formula_id(selected))"))
    end
    for design in problem.system.designs
        Formulation(Val(:pscad), requested.methods.pipe_impedance, design)
    end
    earth = requested.methods.earth_impedance
    identifier = formula_id(earth)
    effective = identifier === :default ?
        (placement === Val(:mixed) ? :Lucca1994 : :DirectNumericalIntegration) : identifier
    pscad_setting(Val(effective), placement)
    isempty(EarthImpedance.routes(earth)) ||
        EarthImpedance.routes(earth) == EarthImpedance.routes(EarthImpedance.Formula(identifier)) ||
        throw(ArgumentError("PSCAD cannot apply custom analytical earth-impedance routes"))
    methods = merge(requested.methods, (earth_impedance=EarthImpedance.Formula(effective),))
    return PSCADFormulation(methods, requested.options, requested.definitions)
end

# Shared identifiers map to one native setting; placement comes from the problem.
Formulation(::Val{:pscad}, ::Val{:default}, ::PipeImpedance.Formula{:default}, ::Val{:coaxial}) = nothing
function Formulation(::Val{:pscad}, ::Val{:default}, ::PipeImpedance.Formula{:default}, ::Val{:pipe})
    throw(ArgumentError("PSCAD Cable_Coax does not support an eccentric or multicore metallic pipe enclosure"))
end

pscad_setting(::Val{:DeriSemlyen1981}, ::Val{:overhead}) =
    (field=:EarthForm2, value=0, readback="DERISEMLYEN")
pscad_setting(::Val{:DirectNumericalIntegration}, ::Val{:overhead}) =
    (field=:EarthForm2, value=2, readback="DIRECT_NUMERICAL_INTEGRATION")
pscad_setting(::Val{:WedepohlWilcox1973}, ::Val{:underground}) =
    (field=:EarthForm, value=0, readback="WEDEPOHL")
pscad_setting(::Val{:DirectNumericalIntegration}, ::Val{:underground}) =
    (field=:EarthForm, value=2, readback="DIRECT_NUMERICAL_INTEGRATION")
pscad_setting(::Val{:Saad1996}, ::Val{:underground}) =
    (field=:EarthForm, value=3, readback="SAAD")
pscad_setting(::Val{:Ametani2009}, ::Val{:mixed}) =
    (field=:EarthForm3, value=0, readback="AMETANIL")
pscad_setting(::Val{:Lucca1994}, ::Val{:mixed}) =
    (field=:EarthForm3, value=2, readback="LUCCA")
function pscad_setting(::Val{ID}, ::Val{Placement}) where {ID, Placement}
    throw(ArgumentError("PSCAD earth-impedance :$ID is not supported for $Placement placement"))
end

"""
    formulas(placement::Val)

Return shared earth-impedance identifiers with a native PSCAD setting for
`Val(:overhead)`, `Val(:underground)` or `Val(:mixed)`. Availability comes from
the backend's setting dispatch over the Engine registry, not a second author
list. No project is generated and no remote solver is invoked.
"""
function formulas(placement::Val)
    placement in (Val(:overhead), Val(:underground), Val(:mixed)) ||
        throw(ArgumentError("PSCAD formula placement must be overhead, underground or mixed"))
    identifiers = Symbol[]
    for identifier in EarthImpedance.REGISTERED
        identifier === :default && continue
        try
            pscad_setting(Val(identifier), placement)
        catch error
            error isa ArgumentError || rethrow()
            continue
        end
        push!(identifiers, identifier)
    end
    return Tuple(identifiers)
end

function pscad_setting(formulation::PSCADFormulation, problem::LineParametersProblem)
    heights = getproperty.(problem.system.positions, :y)
    placement = all(>(0), heights) ? Val(:overhead) :
                all(<(0), heights) ? Val(:underground) : Val(:mixed)
    selected = Formulation(Val(:pscad), problem, formulation)
    return pscad_setting(Val(formula_id(selected.methods.earth_impedance)), placement)
end

function formulation_record(formulation::PSCADFormulation)
    requested = map(formulation.definitions) do definition
        definition isa Symbol ? string(definition) :
        applicable(formula_id, definition) ? string(formula_id(definition)) : repr(definition)
    end
    methods = formulation.methods
    return (
        schema_version=1,
        backend=:pscad,
        type=string(parentmodule(typeof(formulation)), ".", nameof(typeof(formulation))),
        requested,
        raw=Dict{Symbol, Any}(:selections => formulation_record(LineParametersFormulation(
            methods, formulation.options, formulation.definitions))),
        effective=(
            internal_impedance=nothing, insulation_impedance=nothing,
            earth_impedance=formula_id(methods.earth_impedance),
            earth_admittance=nothing,
            insulation_admittance=formula_id(methods.insulation_admittance),
            semicon_admittance=formula_id(methods.semicon_admittance),
            earth_properties=nothing, equivalent_earth=nothing, pipe_impedance=nothing),
        assumptions=(
            internal_impedance="PSCAD native Cable_Coax conductor calculation",
            insulation_impedance="PSCAD native Cable_Coax magnetic calculation",
            earth_admittance="PSCAD native earth-admittance calculation",
            insulation_admittance=description(methods.insulation_admittance),
            semicon_admittance=description(methods.semicon_admittance),
            dielectric_equivalence="Reference-frequency equivalent capacitance and loss tangent; native PSCAD frequency law; loss tangent capped at 10",
            earth="One homogeneous earth layer; no FD relation",
            pipe_impedance="Cable_Coax only; shared eccentric metallic enclosure unsupported"),
        options=formulation.options)
end

function computation_details(::Type{<:PSCADFormulation}, result::LineParameters)::ComputationDetails
    return details(result)
end
