const _GAUNTLET_ROOT = normpath(joinpath(@__DIR__, ".."))
const _SMOKE_ROOT = joinpath(_GAUNTLET_ROOT, "fixtures", "smoke")

"""Open a normalized Gauntlet corpus lazily."""
struct Corpus
    root::String
    cases::Dict{String, String}
    rejections::Dict{String, String}
    aliases::Dict{String, String}
end

function show(io::IO, corpus::Corpus)
    print(io, "Corpus(", length(corpus.cases), " cases, ",
        length(corpus.rejections), " rejections)")
end

function Corpus(path::AbstractString)
    root = abspath(path)
    index_path = joinpath(root, "index.toml")
    isfile(index_path) || throw(ArgumentError("Gauntlet corpus has no index.toml: $root"))
    index = TOML.parsefile(index_path)
    get(index, "schema_version", 0) == 1 || throw(ArgumentError(
        "unsupported Gauntlet corpus schema $(get(index, "schema_version", nothing))",
    ))
    cases = Dict{String, String}(
        String(id) => String(relative) for (id, relative) in get(index, "cases", Dict())
    )
    rejections = Dict{String, String}(
        String(id) => String(relative)
    for (id, relative) in get(index, "rejections", Dict())
    )
    aliases_path = joinpath(root, "aliases.toml")
    aliases = isfile(aliases_path) ?
              Dict{String, String}(
        String(id) => String(target)
    for (id, target) in get(TOML.parsefile(aliases_path), "aliases", Dict())
    ) : Dict{String, String}()
    return Corpus(root, cases, rejections, aliases)
end

function Corpus(name::Symbol)
    name === :smoke ? Corpus(_SMOKE_ROOT) :
    throw(ArgumentError(
        "Corpus(:$name) is not locally bound; use Corpus(path) for a staged full corpus",
    ))
end

length(corpus::Corpus) = length(corpus.cases) + length(corpus.rejections)
keys(corpus::Corpus) = sort!(collect(union(keys(corpus.cases), keys(corpus.rejections))))

function iterate(corpus::Corpus, state = 1)
    ids = keys(corpus)
    state > length(ids) && return nothing
    return corpus[ids[state]], state + 1
end

function _canonical(corpus::Corpus, id::AbstractString)
    candidate = String(id)
    seen = Set{String}()
    while haskey(corpus.aliases, candidate)
        candidate in seen && throw(ArgumentError("cyclic Gauntlet alias at $candidate"))
        push!(seen, candidate)
        candidate = corpus.aliases[candidate]
    end
    return candidate
end

_family(name::AbstractString) = _family(Val(Symbol(name)))

function _fidelity(name::AbstractString)
    symbol = Symbol(name)
    symbol === :exact && return Exact()
    symbol === :approximate && return Approximate()
    symbol === :rejected && return Rejected()
    throw(ArgumentError("unknown fidelity '$name'"))
end

_reduction(name::AbstractString) = _case_reduction(name)

function _provenance(values::AbstractDict)
    return Provenance(
        Symbol(values["source"]),
        String(values["version"]),
        String(values["campaign"]),
        String(values["definition"]),
        String(values["source_sha256"]),
        String(values["specification_sha256"]),
        String(values["case_sha256"]),
        get(values, "elapsed_seconds", nothing),
        Dict{String, String}(
            String(key) => String(value)
        for (key, value) in get(values, "artifact_sha256", Dict())
        ),
        String(get(values, "cohort", values["case_sha256"]))
    )
end

function _assumptions(values)
    return Assumption[Assumption(Symbol(value["subject"]), String(value["detail"]))
                      for value in values]
end

function _dict_block(values::AbstractDict)
    result = PSCADIO.Block(String(values["name"]); value = get(values, "value", nothing))
    merge!(result.fields, Dict{String, Any}(
        String(key) => value for (key, value) in values["fields"]
    ))
    append!(result.children, _dict_block.(values["children"]))
    return result
end

function _line_parameters(values::Nothing, domain_type = LineCableModels.PhaseDomain)
    return nothing
end

function _line_parameters(values::AbstractDict, domain_type = LineCableModels.PhaseDomain)
    return EN.LineParameters(
        domain_type,
        Array{ComplexF64, 3}(values["Z"]),
        Array{ComplexF64, 3}(values["Y"]),
        Float64.(values["frequency"]);
        basis = Symbol(values["basis"])
    )
end

function _fit(values::Nothing)
    return nothing
end

function _fit(values::AbstractDict)
    columns = RationalColumn{Float64}[RationalColumn(
                                          ComplexF64.(column["poles"]),
                                          Matrix{ComplexF64}(column["residues"])
                                      ) for column in values["columns"]]
    groups = DelayGroup{Float64}[DelayGroup(
                                     Float64(group["delay"]),
                                     ComplexF64.(group["poles"]),
                                     Array{ComplexF64, 3}(group["residues"])
                                 ) for group in values["groups"]]
    range = Tuple(Float64.(values["frequency_range"]))
    return Fit(columns, Matrix{Float64}(values["constant"]), groups, range)
end

function _modes(values::Nothing)
    return nothing
end

function _optional_array(values, key, ::Type{T}) where {T}
    value = get(values, key, nothing)
    return value === nothing ? nothing : T(value)
end

function _modes(values::AbstractDict)
    return Modes(
        Float64.(values["frequency"]);
        transform = _optional_array(values, "transform", Array{ComplexF64, 3}),
        propagation = _optional_array(values, "propagation", Matrix{ComplexF64}),
        characteristic = _optional_array(values, "characteristic", Array{ComplexF64, 3}),
        fitted_characteristic = _optional_array(
            values, "fitted_characteristic", Array{ComplexF64, 3}
        ),
        phase_propagation = _optional_array(
            values, "phase_propagation", Array{ComplexF64, 3}
        ),
        modal_propagation = _optional_array(
            values, "modal_propagation", Matrix{ComplexF64}
        ),
        fitted_phase_propagation = _optional_array(
            values, "fitted_phase_propagation", Array{ComplexF64, 3}
        )
    )
end

function _terminal(values::Nothing)
    return nothing
end

function _terminal(values::AbstractDict)
    return Terminal(
        Float64.(values["frequency"]);
        open = _optional_array(values, "open", Matrix{ComplexF64}),
        short = _optional_array(values, "short", Matrix{ComplexF64}),
        empty = Symbol.(get(values, "empty", String[]))
    )
end

function _reference(values::AbstractDict, fallback_ports)
    ports = isempty(get(values, "ports", [])) ? fallback_ports :
            Port[Port(
                     String(port["id"]),
                     Int(port["cable"]),
                     String(port["component"]),
                     Int(port["phase"])
                 ) for port in values["ports"]]
    return Reference(
        phase = _line_parameters(get(values, "phase", nothing)),
        sequence = _line_parameters(
            get(values, "sequence", nothing), LineCableModels.ModalDomain
        ),
        sequence_transform = _optional_array(
            values, "sequence_transform", Matrix{ComplexF64}
        ),
        modes = _modes(get(values, "modes", nothing)),
        fit = _fit(get(values, "fit", nothing)),
        terminal = _terminal(get(values, "terminal", nothing)),
        ports = ports
    )
end

function _checks(reference::Reference)
    values = Check[]
    reference.phase === nothing || append!(values, (MatrixCheck{:Z}(), MatrixCheck{:Y}()))
    reference.sequence === nothing || append!(
        values, (ModalCheck{:Z}(), ModalCheck{:Y}())
    )
    if reference.modes !== nothing
        reference.modes.transform === nothing || append!(values, (
            ModalCheck{:transform}(), ModalCheck{:reconstruction}()
        ))
        reference.modes.propagation === nothing || push!(
            values, ModalCheck{:propagation}()
        )
        reference.modes.characteristic === nothing || push!(
            values, ModalCheck{:characteristic}()
        )
        reference.modes.phase_propagation === nothing || push!(
            values, ModalCheck{:phase_propagation}()
        )
        reference.modes.modal_propagation === nothing || push!(
            values, ModalCheck{:modal_propagation}()
        )
    end
    reference.fit === nothing || append!(values, (FitCheck{:Yc}(), FitCheck{:H}()))
    reference.terminal === nothing || push!(values, PhysicalCheck{:terminal}())
    append!(values,
        (
            PhysicalCheck{:symmetry}(),
            PhysicalCheck{:reciprocity}(),
            PhysicalCheck{:positive_real}()
        ))
    return tuple(values...)
end

function _load_success(path::AbstractString)
    payload = JLD2.load(path, "record")
    id = String(payload["id"])
    family = _family(String(payload["family"]))
    fidelity = _fidelity(String(payload["fidelity"]))
    reduction = _reduction(String(payload["reduction"]))
    input = PSCADIO.ParsedInput(Symbol(payload["input_kind"]), _dict_block(payload["input"]))
    problem, fallback_ports = materialize(family, input, reduction, id)
    reference_value = _reference(payload["reference"], fallback_ports)
    assumptions_value = _assumptions(payload["assumptions"])
    formulation = _formulation(input, reduction)
    return Case(
        id,
        family,
        fidelity,
        reduction,
        problem,
        formulation,
        reference_value,
        _checks(reference_value),
        assumptions_value,
        _provenance(payload["provenance"])
    )
end

function _load_rejection(path::AbstractString)
    payload = TOML.parsefile(path)
    reference_value = Reference()
    return Case(
        String(payload["id"]),
        _family(String(payload["family"])),
        Rejected(),
        _reduction(String(payload["reduction"])),
        nothing,
        nothing,
        reference_value,
        (),
        _assumptions(get(payload, "assumptions", [])),
        _provenance(payload["provenance"])
    )
end

function getindex(corpus::Corpus, requested::AbstractString)
    id = _canonical(corpus, requested)
    if haskey(corpus.cases, id)
        return _load_success(joinpath(corpus.root, corpus.cases[id]))
    elseif haskey(corpus.rejections, id)
        return _load_rejection(joinpath(corpus.root, corpus.rejections[id]))
    end
    throw(KeyError(requested))
end

function Suite(
        name::Symbol;
        corpus = Corpus(name),
        ids = keys(corpus),
        checks = nothing,
        policy::Policy = ExactOnly(),
        tolerances = Dict{String, Tolerance{Float64}}(),
        performance = String[]
)
    selected = String.(ids)
    all(id -> id in keys(corpus) || haskey(corpus.aliases, id), selected) ||
        throw(KeyError("suite includes a case absent from its corpus"))
    performance_ids = performance === true ? copy(selected) :
                      performance === false ? String[] : String.(performance)
    all(id -> id in selected, performance_ids) || throw(KeyError(
        "performance selection must be a subset of the suite",
    ))
    converted_tolerances = Dict{String, Tolerance{Float64}}(
        String(key) => Tolerance(Float64(value.rtol), Float64(value.atol))
    for (key, value) in pairs(tolerances)
    )
    return Suite(
        name,
        corpus,
        selected,
        checks === nothing ? nothing : tuple(checks...),
        policy,
        converted_tolerances,
        performance_ids
    )
end
