const _GAUNTLET_ROOT = normpath(joinpath(@__DIR__, ".."))
const _SMOKE_ROOT = joinpath(_GAUNTLET_ROOT, "fixtures", "smoke")

abstract type DatasetRecord end

struct CaseRecord <: DatasetRecord
    path::String
end

struct RejectionRecord <: DatasetRecord
    path::String
end

"""Open a normalized Gauntlet dataset lazily."""
struct Dataset{D <: Datasource}
    datasource::D
    root::String
    cases::Dict{String, String}
    rejections::Dict{String, String}
    aliases::Dict{String, String}
end

function show(io::IO, dataset::Dataset)
    print(io, "Dataset(", nameof(typeof(dataset.datasource)), ", ",
        length(dataset.cases), " cases, ", length(dataset.rejections), " rejections)")
end

function Dataset(path::AbstractString)
    root = abspath(path)
    index_path = joinpath(root, "index.toml")
    isfile(index_path) || throw(ArgumentError("Gauntlet dataset has no index.toml: $root"))
    index = TOML.parsefile(index_path)
    get(index, "schema_version", 0) == 1 || throw(ArgumentError(
        "unsupported Gauntlet dataset schema $(get(index, "schema_version", nothing))",
    ))
    haskey(index, "datasource") || throw(ArgumentError(
        "Gauntlet dataset index does not declare a datasource",
    ))
    source = datasource(Symbol(index["datasource"]))
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
    return Dataset(source, root, cases, rejections, aliases)
end

function Dataset(name::Symbol)
    name === :smoke ? Dataset(_SMOKE_ROOT) :
    throw(ArgumentError(
        "Dataset(:$name) is not locally bound; use Dataset(path) for a staged full dataset",
    ))
end

length(dataset::Dataset) = length(dataset.cases) + length(dataset.rejections)
function keys(dataset::Dataset)
    sort!(collect(union(keys(dataset.cases), keys(dataset.rejections))))
end

function iterate(dataset::Dataset, state = 1)
    ids = keys(dataset)
    state > length(ids) && return nothing
    return dataset[ids[state]], state + 1
end

function _canonical(dataset::Dataset, id::AbstractString)
    candidate = String(id)
    seen = Set{String}()
    while haskey(dataset.aliases, candidate)
        candidate in seen && throw(ArgumentError("cyclic Gauntlet alias at $candidate"))
        push!(seen, candidate)
        candidate = dataset.aliases[candidate]
    end
    return candidate
end

function getindex(dataset::Dataset, requested::AbstractString)
    id = _canonical(dataset, requested)
    if haskey(dataset.cases, id)
        return load(
            dataset.datasource,
            CaseRecord(joinpath(dataset.root, dataset.cases[id]))
        )
    elseif haskey(dataset.rejections, id)
        return load(
            dataset.datasource,
            RejectionRecord(joinpath(dataset.root, dataset.rejections[id]))
        )
    end
    throw(KeyError(requested))
end

function Suite(
        name::Symbol;
        dataset = Dataset(name),
        ids = keys(dataset),
        checks = nothing,
        policy::Policy = ExactOnly(),
        tolerances = Dict{String, Tolerance{Float64}}(),
        performance = String[]
)
    selected = String.(ids)
    all(id -> id in keys(dataset) || haskey(dataset.aliases, id), selected) ||
        throw(KeyError("suite includes a case absent from its dataset"))
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
        dataset,
        selected,
        checks === nothing ? nothing : tuple(checks...),
        policy,
        converted_tolerances,
        performance_ids
    )
end
