_engine_type_name(value) = String(nameof(typeof(value)))
_engine_unit(unit) = replace(Units.label(unit), "." => "·")

_domain_name(::Type{PhaseDomain}) = "phase domain"
_domain_name(::Type{ModalDomain}) = "sequence domain"
_domain_name(::Type{D}) where {D <: LineParamsDomain} = lowercase(String(nameof(D)))

TextDisplay.@showfields EarthPair "EarthPair" pair -> (
    row = pair.row,
    column = pair.column,
    heights = TextDisplay.engineering.(pair.heights, Ref(:meter)),
    separation = TextDisplay.engineering(pair.separation, :meter),
    layers = pair.layers
)

TextDisplay.@showfields BlueprintConductor "BlueprintConductor" row -> (
    terminal = row.terminal,
    assembly = row.assembly,
    r_in = TextDisplay.engineering(row.r_in, :meter),
    r_ex = TextDisplay.engineering(row.r_ex, :meter),
    kind = row.material.kind
)

TextDisplay.@showfields BlueprintDielectric "BlueprintDielectric" row -> (
    conductor = row.conductor,
    r_in = TextDisplay.engineering(row.r_in, :meter),
    r_ex = TextDisplay.engineering(row.r_ex, :meter),
    kind = row.material.kind
)

TextDisplay.@showfields CableBlueprint "CableBlueprint" blueprint -> (
    cable_id = blueprint.cable_id,
    conductors = length(blueprint.conductors),
    dielectrics = length(blueprint.dielectrics),
    assemblies = length(blueprint.assembly_ranges)
)

TextDisplay.@showfields ConsoleVerbosityLogger "ConsoleVerbosityLogger" logger -> (
    sink = String(nameof(typeof(logger.console))),
    levels = logger.levels
)

TextDisplay.@showfields LineCableModelsFEMError "LineCableModelsFEMError" error -> (
    category = error.category,
    object = error.object_id,
    field = error.field,
    message = error.message,
    run_directory = error.run_directory
)

function _frequency_span(values)
    isempty(values) && return "no points"
    count = length(values)
    first_value = TextDisplay.engineering(first(values), :hertz)
    last_value = TextDisplay.engineering(last(values), :hertz)
    return count == 1 ? "1 point · $first_value" :
           "$count points · $first_value … $last_value"
end

TextDisplay.name(::Type{LineCableModelsCoaxial}) = "LineCableModelsCoaxial"
Base.summary(io::IO, ::LineCableModelsCoaxial) =
    print(io, "LineCableModels coaxial backend")
Base.show(io::IO, ::LineCableModelsCoaxial) = print(io, "LineCableModelsCoaxial()")
Base.show(io::IO, ::MIME"text/plain", backend::LineCableModelsCoaxial) =
    show(io, backend)

TextDisplay.name(::Type{<:LineCableModelsFEM}) = "LineCableModelsFEM"
Base.summary(io::IO, ::LineCableModelsFEM) = print(io, "LineCableModels FEM backend")
function Base.show(io::IO, backend::LineCableModelsFEM)
    print(io, "LineCableModelsFEM(options=", length(backend.options), ")")
end
function Base.show(io::IO, ::MIME"text/plain", backend::LineCableModelsFEM)
    get(io, :compact, false) && return show(io, backend)
    return TextDisplay.fields(
        io,
        "LineCableModels FEM backend",
        (options = length(backend.options), execution = backend.execution);
        multiline = true
    )
end

TextDisplay.name(::Type{LineCableModelsFEMOptions}) = "LineCableModelsFEMOptions"
Base.summary(io::IO, ::LineCableModelsFEMOptions) = print(io, "FEM execution options")
function Base.show(io::IO, options::LineCableModelsFEMOptions)
    print(io, "LineCableModelsFEMOptions(mesh_policy=:", options.mesh_policy, ")")
end
function Base.show(io::IO, ::MIME"text/plain", options::LineCableModelsFEMOptions)
    get(io, :compact, false) && return show(io, options)
    return TextDisplay.fields(
        io,
        "FEM execution options",
        (
            ui = options.ui,
            plot_field_maps = options.plot_field_maps,
            mesh_policy = options.mesh_policy,
            mesh_path = options.mesh_path,
            keep_run_directory = options.keep_run_directory,
            getdp_executable = options.getdp_executable,
            gmsh_verbosity = options.gmsh_verbosity,
            getdp_verbosity = options.getdp_verbosity,
        );
        multiline = true
    )
end

TextDisplay.name(::Type{PhaseDomain}) = "Phase domain"
Base.summary(io::IO, ::PhaseDomain) = print(io, "Phase domain")
Base.show(io::IO, ::PhaseDomain) = print(io, "PhaseDomain()")
Base.show(io::IO, ::MIME"text/plain", value::PhaseDomain) = show(io, value)

TextDisplay.name(::Type{ModalDomain}) = "Sequence domain"
Base.summary(io::IO, ::ModalDomain) = print(io, "Sequence domain")
Base.show(io::IO, ::ModalDomain) = print(io, "ModalDomain()")
Base.show(io::IO, ::MIME"text/plain", value::ModalDomain) = show(io, value)

function Base.summary(io::IO, formulation::AbstractFormulation)
    print(io, _engine_type_name(formulation), " formulation")
end
function Base.show(io::IO, formulation::AbstractFormulation)
    print(io, _engine_type_name(formulation), "()")
end
Base.show(io::IO, ::MIME"text/plain", formulation::AbstractFormulation) =
    show(io, formulation)

TextDisplay.name(::Type{<:LineParametersFormulation}) = "LineParametersFormulation"
Base.summary(io::IO, ::LineParametersFormulation) = print(io, "Line-parameters formulation")
function Base.show(io::IO, formulation::LineParametersFormulation)
    print(io, "LineParametersFormulation(", length(formulation.methods), " methods)")
end

TextDisplay.name(::Type{<:CableConstantsFormulation}) = "CableConstantsFormulation"
Base.summary(io::IO, ::CableConstantsFormulation) =
    print(io, "Cable-constants formulation")
function Base.show(io::IO, formulation::CableConstantsFormulation)
    print(io, "CableConstantsFormulation(", length(formulation.methods), " methods)")
end
function Base.show(io::IO, ::MIME"text/plain", formulation::CableConstantsFormulation)
    get(io, :compact, false) && return show(io, formulation)
    children = Tuple((
        label = string(key, "  ", sprint(show, method; context = :compact => true)),
        noun = "methods",
    ) for (key, method) in pairs(formulation.methods))
    return TextDisplay.tree(io, "Cable-constants formulation", children; noun = "methods")
end

TextDisplay.name(::Type{<:CableConstantsProblem}) = "CableConstantsProblem"
Base.summary(io::IO, ::CableConstantsProblem) = print(io, "Cable-constants problem")
function Base.show(io::IO, problem::CableConstantsProblem)
    print(io, "CableConstantsProblem(temperature=", problem.temperature,
        ", frequency=", problem.frequency, ")")
end
function Base.show(io::IO, ::MIME"text/plain", problem::CableConstantsProblem)
    get(io, :compact, false) && return show(io, problem)
    return TextDisplay.tree(io, "CableConstantsProblem", (
        (label = "design       $(problem.design.cable_id)", noun = "fields"),
        (label = "temperature  $(TextDisplay.engineering(problem.temperature, :celsius))", noun = "fields"),
        (label = "frequency    $(TextDisplay.engineering(problem.frequency, :hertz))", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:CableConstants}) = "CableConstants"
Base.summary(io::IO, constants::CableConstants) =
    print(io, "CableConstants, ", length(constants), " assemblies")
function Base.show(io::IO, constants::CableConstants)
    print(io, "CableConstants(assemblies=", length(constants),
        ", frequency=", constants.frequency, ")")
end
function Base.show(io::IO, ::MIME"text/plain", constants::CableConstants)
    get(io, :compact, false) && return show(io, constants)
    rows = Tuple((
        label = string(
            constants.cores[index], "  R=", constants.R[index],
            " Ω/m  L=", constants.L[index], " H/m  C=", constants.C[index],
            " F/m  G=", constants.G[index], " S/m"
        ),
        noun = "assemblies",
    ) for index in eachindex(constants.cores))
    return TextDisplay.tree(io,
        "Cable constants · $(constants.frequency) Hz", rows; noun = "assemblies")
end
function Base.show(io::IO, ::MIME"text/plain", formulation::LineParametersFormulation)
    get(io, :compact, false) && return show(io, formulation)
    children = Any[(
        label = string(key, "  ", sprint(show, method; context = :compact => true)),
        noun = "methods",
    ) for (key, method) in pairs(formulation.methods)]
    isempty(formulation.options) || push!(children, (
        label = "options  $(length(formulation.options)) entries",
        noun = "methods",
    ))
    return TextDisplay.tree(io, "Line-parameters formulation", Tuple(children); noun = "methods")
end

TextDisplay.name(::Type{<:LineParametersProblem}) = "LineParametersProblem"
Base.summary(io::IO, problem::LineParametersProblem) =
    print(io, "Line-parameters problem over $(length(problem.frequencies)) frequencies")
function Base.show(io::IO, problem::LineParametersProblem)
    print(io, "LineParametersProblem(cables=", ncables(problem.system),
        ", frequencies=", length(problem.frequencies), ")")
end
function Base.show(io::IO, ::MIME"text/plain", problem::LineParametersProblem)
    get(io, :compact, false) && return show(io, problem)
    children = (
        (label = "system       $(sprint(show, problem.system; context = :compact => true))", noun = "fields"),
        (label = "temperature  $(TextDisplay.engineering(problem.temperature, :celsius))", noun = "fields"),
        (label = "earth        $(sprint(show, problem.earth_props; context = :compact => true))", noun = "fields"),
        (label = "f            $(_frequency_span(problem.frequencies))", noun = "fields"),
    )
    return TextDisplay.tree(io, "LineParametersProblem", children)
end

function _array_summary(name, value, selector)
    return string(
        name, "(", join(size(value), '×'), "; unit=",
        _engine_unit(_result_unit(value, selector)), ")"
    )
end

TextDisplay.name(::Type{<:SeriesImpedance}) = "SeriesImpedance"
Base.summary(io::IO, value::SeriesImpedance) =
    print(io, "Series impedance, ", join(size(value), '×'))
Base.show(io::IO, value::SeriesImpedance) =
    print(io, _array_summary("SeriesImpedance", value, Z))
function Base.show(io::IO, ::MIME"text/plain", value::SeriesImpedance)
    get(io, :compact, false) && return show(io, value)
    TextDisplay.tree(io, "Series impedance · $(join(size(value), '×'))", (
        (label = "unit  $(_engine_unit(_result_unit(value, Z)))", noun = "fields"),
        (label = "values summarized; use observables for extraction", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:ShuntAdmittance}) = "ShuntAdmittance"
Base.summary(io::IO, value::ShuntAdmittance) =
    print(io, "Shunt admittance, ", join(size(value), '×'))
Base.show(io::IO, value::ShuntAdmittance) =
    print(io, _array_summary("ShuntAdmittance", value, Y))
function Base.show(io::IO, ::MIME"text/plain", value::ShuntAdmittance)
    get(io, :compact, false) && return show(io, value)
    TextDisplay.tree(io, "Shunt admittance · $(join(size(value), '×'))", (
        (label = "unit  $(_engine_unit(_result_unit(value, Y)))", noun = "fields"),
        (label = "values summarized; use observables for extraction", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:LineParameters}) = "LineParameters"
function Base.summary(io::IO, parameters::LineParameters)
    print(io, "LineParameters, ", nconductors(parameters), '×', nconductors(parameters),
        " over ", nfrequencies(parameters), " frequencies")
end
function Base.show(io::IO, parameters::LineParameters)
    print(io, "LineParameters(", _domain_name(domain(parameters)), "; ",
        nconductors(parameters), '×', nconductors(parameters), '×',
        nfrequencies(parameters), ", basis=:", basis(parameters), ")")
    has_uncertainty_type(eltype(parameters)) && print(io, " ±")
end
function Base.show(io::IO, ::MIME"text/plain", parameters::LineParameters)
    get(io, :compact, false) && return show(io, parameters)
    zunit = _engine_unit(_result_unit(parameters, Z))
    yunit = _engine_unit(_result_unit(parameters, Y))
    shape = join(size(parameters.Z), '×')
    children = (
        (label = "f  $(_frequency_span(parameters.f))", noun = "fields"),
        (label = "Z  $shape · $zunit", noun = "fields"),
        (label = "Y  $shape · $yunit", noun = "fields"),
    )
    return TextDisplay.tree(io, "LineParameters · $(_domain_name(domain(parameters)))", children)
end

TextDisplay.name(::Type{<:RMSError}) = "RMSError"
Base.summary(io::IO, error::RMSError) =
    print(io, "RMS error, ", join(size(error.absolute), '×'))
function Base.show(io::IO, error::RMSError)
    print(io, "RMSError(", join(size(error.absolute), '×'), ")")
end
function Base.show(io::IO, ::MIME"text/plain", error::RMSError)
    get(io, :compact, false) && return show(io, error)
    return TextDisplay.tree(io, "RMS error · $(join(size(error.absolute), '×'))", (
        (label = "absolute  matrix summarized", noun = "fields"),
        (label = "relative  matrix summarized", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:LineParametersBenchmark}) = "LineParametersBenchmark"
Base.summary(io::IO, benchmark::LineParametersBenchmark) =
    print(io, "Line-parameters benchmark, ", join(size(benchmark.Z.absolute), '×'))
function Base.show(io::IO, benchmark::LineParametersBenchmark)
    print(io, "LineParametersBenchmark(", join(size(benchmark.Z.absolute), '×'),
        "; basis=:", basis(benchmark), ")")
end
function Base.show(io::IO, ::MIME"text/plain", benchmark::LineParametersBenchmark)
    get(io, :compact, false) && return show(io, benchmark)
    shape = join(size(benchmark.Z.absolute), '×')
    return TextDisplay.tree(io, "Line-parameters benchmark", (
        (label = "Z  $shape RMS errors", noun = "fields"),
        (label = "Y  $shape RMS errors", noun = "fields"),
        (label = "basis  :$(basis(benchmark))", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:LineParametersWorkspace}) = "LineParametersWorkspace"
Base.summary(io::IO, workspace::LineParametersWorkspace) =
    print(io, "Line-parameters workspace")
function Base.show(io::IO, workspace::LineParametersWorkspace)
    input = workspace.input
    print(io, "LineParametersWorkspace(phases=", input.n_phases,
        ", cables=", input.n_cables,
        ", frequencies=", input.n_frequencies, ")")
end
function Base.show(io::IO, ::MIME"text/plain", workspace::LineParametersWorkspace)
    get(io, :compact, false) && return show(io, workspace)
    input = workspace.input
    return TextDisplay.tree(io, "Line-parameters workspace", (
        (label = "phases       $(input.n_phases)", noun = "fields"),
        (label = "cables       $(input.n_cables)", noun = "fields"),
        (label = "frequencies  $(input.n_frequencies)", noun = "fields"),
        (label = "capture      $(workspace.capture === nothing ? "disabled" : "enabled")", noun = "fields"),
    ))
end
