module NumericalReferences

using JLD2
using LineCableModels
using LineCableModels.Engine: compare
using Pkg.Artifacts
using SHA
using TOML
using Test: @inferred

# This reader deliberately does not load a case catalogue or any backend adapter.
# The reviewed file supplies the physical declarations; the current package
# supplies the implementation under test.
function read_reference(path::AbstractString, expected_sha256::AbstractString)
    occursin(r"^[0-9a-f]{64}$", expected_sha256) || throw(ArgumentError(
        "numerical reference needs a reviewed SHA-256 digest: $path"))
    isfile(path) || throw(ArgumentError("numerical reference file is missing: $path"))
    bytes2hex(open(sha256, path)) == expected_sha256 || throw(ArgumentError(
        "numerical reference checksum mismatch: $path"))
    document = JLD2.load(path)
    required = ("kind", "status", "backend", "problem", "formulation", "Z", "Y",
        "frequencies", "port_order", "basis", "domain")
    missing = filter(name -> !haskey(document, name), required)
    isempty(missing) || throw(ArgumentError(
        "numerical reference $path lacks $(join(missing, ", ")); historical reports " *
        "remain readable, but replay requires stored problem and formulation declarations"))
    document["kind"] === :gauntlet_calculation && document["status"] === :complete ||
        throw(ArgumentError("numerical reference must contain completed phase matrices: $path"))
    document["backend"] === :coaxial || throw(ArgumentError(
        "CI replays owned coaxial calculations, not external solvers; select a reviewed " *
        "coaxial result validated through the backend comparisons: $path"))
    document["basis"] === :pul && document["domain"] === :PhaseDomain ||
        throw(ArgumentError("numerical reference requires phase-domain Z/Y per metre: $path"))
    problem = LineCableModels.ImportExport.deserialize_value(document["problem"])
    problem isa LineParametersProblem || throw(ArgumentError(
        "numerical reference must declare one scalar LineParametersProblem: $path"))
    declared = document["formulation"]
    declared isa NamedTuple && hasproperty(declared, :definitions) &&
        declared.definitions isa NamedTuple &&
        hasproperty(declared, :options) && declared.options isa NamedTuple || throw(ArgumentError(
        "numerical reference requires author selectors and explicit formulation options: $path"))
    formulation = Formulation(; declared.definitions..., options=declared.options)
    formulation isa LineParametersFormulation || throw(ArgumentError(
        "numerical reference must declare one scalar formulation, not a formulation grid: $path"))
    all(isfinite, document["Z"]) && all(isfinite, document["Y"]) || throw(ArgumentError(
        "numerical reference contains nonfinite matrix entries: $path"))
    parameters = LineParameters(PhaseDomain, document["Z"], document["Y"], document["frequencies"];
        basis=:pul)
    parameters.f == problem.frequencies || throw(ArgumentError(
        "numerical reference frequency samples differ from its stored problem: $path"))
    length(document["port_order"]) == size(parameters.Z, 1) && allunique(document["port_order"]) ||
        throw(ArgumentError("numerical reference terminal labels do not match its matrix axes: $path"))
    return (; problem, formulation, parameters)
end

function compare_reference(reference)
    actual = @inferred compute(reference.problem, reference.formulation)
    # Scientific comparison's negligible-signal policy is not CI acceptance.
    # Only the explicitly reviewed manifest tolerances may accept a difference.
    return compare(reference.parameters, actual; atol=0.0)
end

function check(manifest::AbstractString=joinpath(@__DIR__, "approved.toml");
        artifacts_toml::AbstractString=joinpath(dirname(manifest), "Artifacts.toml"))
    document = TOML.parsefile(manifest)
    get(document, "schema_version", nothing) == 1 || throw(ArgumentError(
        "unsupported numerical-reference manifest: $manifest"))
    entries = get(document, "references", [])
    entries isa AbstractVector || throw(ArgumentError(
        "numerical references must be an array of reviewed entries in $manifest"))
    isempty(entries) && throw(ArgumentError(
        "no numerical references have been approved in $manifest; review and pin them " *
        "explicitly as described in test/numerical/README.md. No Gauntlet run will be launched."))
    required = ("id", "artifact", "file", "sha256", "review", "Z_atol", "Y_atol", "rtol")
    for entry in entries
        entry isa AbstractDict && all(name -> haskey(entry, name), required) || throw(ArgumentError(
            "numerical-reference entries require $(join(required, ", "))"))
        for name in ("id", "artifact", "file", "sha256", "review")
            entry[name] isa AbstractString && !isempty(strip(entry[name])) || throw(ArgumentError(
                "numerical reference requires a nonempty $name record"))
        end
        for name in ("Z_atol", "Y_atol", "rtol")
            value = entry[name]
            value isa Real && !(value isa Bool) && isfinite(value) && value >= 0 ||
                throw(ArgumentError("$(entry["id"]) requires a finite nonnegative $name"))
        end
    end
    identifiers = [entry["id"] for entry in entries]
    allunique(identifiers) || throw(ArgumentError("duplicate numerical-reference IDs in $manifest"))
    rows = NamedTuple[]
    for entry in entries
        hash = artifact_hash(entry["artifact"], artifacts_toml)
        hash === nothing && throw(ArgumentError(
            "reviewed artifact $(entry["artifact"]) is not bound in $artifacts_toml"))
        ensure_artifact_installed(entry["artifact"], artifacts_toml)
        root = realpath(artifact_path(hash))
        relative = normpath(entry["file"])
        (isabspath(relative) || relative == ".." || startswith(relative, ".." * Base.Filesystem.path_separator)) &&
            throw(ArgumentError("numerical reference file must be inside its pinned artifact"))
        path = joinpath(root, relative)
        isfile(path) || throw(ArgumentError("reviewed numerical reference is missing: $path"))
        startswith(realpath(path), root * Base.Filesystem.path_separator) || throw(ArgumentError(
            "numerical reference file resolves outside its pinned artifact"))
        reference = read_reference(path, entry["sha256"])
        comparison = compare_reference(reference)
        for (quantity, tolerance) in ((:Z, entry["Z_atol"]), (:Y, entry["Y_atol"]))
            errors = getproperty(comparison, quantity)
            for index in CartesianIndices(errors.absolute)
                absolute, relative = errors.absolute[index], errors.relative[index]
                passed = absolute <= tolerance || relative <= entry["rtol"]
                push!(rows, (id=entry["id"], quantity, row=index[1], column=index[2],
                    passed, absolute, relative))
            end
        end
    end
    return rows
end

end
