const REPOSITORY_ROOT = pkgdir(LineCableModels)

_selection_value(value::Union{Nothing, Bool, Number, AbstractString, Symbol}) = value
function _selection_value(value::Type)
    sprint(show, value; context = (:module=>nothing, :compact=>false))
end
function _selection_value(value::NamedTuple{Names}) where {Names}
    record = map(_selection_value, value)
    # Keep declared record fields when their portable values still fit. Narrowing
    # a NamedTuple field to each point's options breaks homogeneous result spaces.
    types = map(fieldtypes(typeof(value)), values(record)) do declared, item
        item isa declared ? declared : typeof(item)
    end
    return NamedTuple{Names, Tuple{types...}}(record)
end
_selection_value(value::Tuple) = map(_selection_value, value)
_selection_value(value::AbstractArray) = map(_selection_value, value)
function _selection_value(value::AbstractDict)
    Dict(_selection_value(k)=>_selection_value(v) for (k, v) in value)
end
_selection_value(value::Missing) = nothing
function _selection_value(value)
    ismutabletype(typeof(value)) && throw(ArgumentError(
        "calculation records cannot encode mutable $(typeof(value)); provide explicit serializable inputs"))
    names = fieldnames(typeof(value))
    return (type = _selection_value(typeof(value)),
        fields = NamedTuple{names}(
            map(name -> _selection_value(getfield(value, name)), names)))
end

"Return a digest of the materialised numerical declarations in one problem."
function numerical_input_sha256(problem::Engine.LineParametersProblem)
    # JSON transport retains the declaration, not the resolved physical tree.
    # Both are inputs to reuse: a placement change must not resurrect results
    # computed from the same declaration by an older resolver.
    resolved = map(problem.system.designs) do design
        _selection_value((geometry = design.geometry,
            terminal_order = design.terminal_order, terminal_map = design.terminal_map))
    end
    return semantic_sha256((schema_version = 2,
        declaration = ImportExport.serialize_value(problem), resolved))
end

function repository_revision()
    command = `git -C $REPOSITORY_ROOT rev-parse HEAD`
    success(pipeline(command; stdout = devnull, stderr = devnull)) ||
        return (commit = nothing, dirty = nothing)
    commit = readchomp(command)
    dirty = !isempty(readchomp(`git -C $REPOSITORY_ROOT status --porcelain`))
    return (; commit, dirty)
end

function formulation_record(formulation::Engine.LineParametersFormulation)
    return _selection_value(NamedTuple(formulation))
end

function _selection_value(value::Union{Engine.LineParametersFormulation,Engine.LineCableModelsFEM,PSCAD.PSCADFormulation})
    return _selection_value(NamedTuple(value))
end

"Capture the runtime sources and environment used by a campaign."
function implementation_record()
    paths = String[]
    for folder in ("src", "ext", "gauntlet")
        for (directory, directories, names) in walkdir(joinpath(REPOSITORY_ROOT, folder))
            filter!(
                name -> !startswith(name, ".") && !(folder == "gauntlet" &&
                          directory == joinpath(REPOSITORY_ROOT, folder) &&
                          name in ("cases", "benchmarks")),
                directories)
            for name in names
                splitext(name)[2] in (".jl", ".pro", ".py", ".ps1", ".toml") || continue
                name in ("local.jl", "Artifacts.toml") && continue
                path = joinpath(directory, name)
                any(part -> startswith(part, "."), splitpath(relpath(path, REPOSITORY_ROOT))) &&
                    continue
                push!(paths, path)
            end
        end
    end
    for name in ("Project.toml", "Manifest.toml")
        path = joinpath(REPOSITORY_ROOT, name)
        isfile(path) && push!(paths, path)
    end
    project=Base.active_project()
    if project !== nothing
        for path in (project,joinpath(dirname(project),"Manifest.toml"),
                joinpath(dirname(project),"Manifest-v$(VERSION.major).$(VERSION.minor).toml"))
            isfile(path) && push!(paths,path)
        end
    end
    return [(path = relpath(path, REPOSITORY_ROOT),
                sha256 = bytes2hex(open(sha256, path)),
                source = read(path)) for path in sort!(unique(paths))]
end

function implementation_record(formulation; external_sources = ())
    source=implementation_record()
    for relative in external_sources
        path=joinpath(REPOSITORY_ROOT, relative)
        any(row -> row.path == relative, source) && continue
        push!(source,
            (path = String(relative),
                sha256 = bytes2hex(open(sha256, path)), source = read(path)))
    end
    selection=_selection_value(formulation)
    return (; selection, selection_sha256 = semantic_sha256(selection), sources = source)
end

function calculation_record(calculation::BenchmarkCalculation)
    options=Base.structdiff(calculation.options, (;
        on_result = get(calculation.options, :on_result, nothing)))
    problem=calculation.problem
    numerical=problem isa LineParametersProblem ? numerical_input_sha256(problem) :
              semantic_sha256(_selection_value(problem))
    return (id = calculation.id, input_sha256 = numerical,
        formulation = _selection_value(calculation.formulation), options = _selection_value(options))
end
