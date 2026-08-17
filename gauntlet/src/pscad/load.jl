_family(name::AbstractString) = _family(Val(Symbol(name)))

struct PSCADCase{T}
    payload::T
end

struct PSCADRejection{T}
    payload::T
end

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
        Symbol(values["datasource"]),
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

_line_parameters(::Nothing, domain_type = LineCableModels.PhaseDomain) = nothing

function _line_parameters(values::AbstractDict, domain_type = LineCableModels.PhaseDomain)
    return EN.LineParameters(
        domain_type,
        Array{ComplexF64, 3}(values["Z"]),
        Array{ComplexF64, 3}(values["Y"]),
        Float64.(values["frequency"]);
        basis = Symbol(values["basis"])
    )
end

_fit(::Nothing) = nothing

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

_modes(::Nothing) = nothing

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

_terminal(::Nothing) = nothing

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

function load(source::PSCAD, record::CaseRecord)
    return load(source, PSCADCase(JLD2.load(record.path, "record")))
end

function load(source::PSCAD, record::RejectionRecord)
    return load(source, PSCADRejection(TOML.parsefile(record.path)))
end

function load(source::PSCAD, path::AbstractString)
    extension = lowercase(splitext(path)[2])
    extension == ".jld2" && return load(source, CaseRecord(abspath(path)))
    extension == ".toml" && return load(source, RejectionRecord(abspath(path)))
    throw(ArgumentError("PSCAD dataset records must use .jld2 or .toml"))
end

function _pscad_provenance(values)
    result = _provenance(values)
    result.datasource === :pscad || throw(ArgumentError(
        "PSCAD record declares datasource :$(result.datasource)",
    ))
    return result
end

function load(source::PSCAD, record::PSCADCase)
    payload = record.payload
    fidelity = _fidelity(String(payload["fidelity"]))
    id = String(payload["id"])
    family = _family(String(payload["family"]))
    reduction = _reduction(String(payload["reduction"]))
    input = PSCADIO.ParsedInput(Symbol(payload["input_kind"]), _dict_block(payload["input"]))
    problem, fallback_ports = materialize(family, input, reduction, id)
    reference_value = decode(source, record, fallback_ports)
    return Case(
        id,
        family,
        fidelity,
        reduction,
        problem,
        _formulation(input, reduction),
        reference_value,
        _checks(reference_value),
        _assumptions(payload["assumptions"]),
        _pscad_provenance(payload["provenance"])
    )
end

function load(source::PSCAD, record::PSCADRejection)
    payload = record.payload
    reference_value = decode(source, record)
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
        _pscad_provenance(payload["provenance"])
    )
end

function decode(::PSCAD, record::PSCADCase, fallback_ports = Port[])
    _reference(record.payload["reference"], fallback_ports)
end

decode(::PSCAD, ::PSCADRejection, fallback_ports = Port[]) = Reference()
