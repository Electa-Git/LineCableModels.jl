# The test taxonomy, defined once. `test/support/runner.jl`, the architecture guards
# and the tools under `test/tools/` include this file. It uses only Base.
#
# Each item outside `quality` has exactly one owner tag, the latest layer, in
# load order, whose code the item exercises. A change in one layer can affect only
# that layer and the layers loaded after it.

# Owner tags in load order: the core layers, then the extensions.
const OWNERS = (:units, :commons, :materials, :earth, :datamodel, :engine,
    :parametric, :uq, :report, :importexport, :pscad,
    :measurements, :distributions, :xlsx, :fem, :makie)

# The core modules in load order, each with its owner tag. Modules without a tag of
# their own take the tag of the preceding owner.
const MODULE_OWNERS = (Units = :units, Commons = :commons, TextDisplay = :commons,
    PlotBuilder = :commons, Materials = :materials, Earth = :earth,
    DataModel = :datamodel, Engine = :engine,
    ParametricBuilder = :parametric, UQ = :uq, ReportBuilder = :report,
    ImportExport = :importexport, PSCAD = :pscad)

const EXTENSION_OWNERS = (LineCableModelsMeasurementsExt = :measurements,
    LineCableModelsDistributionsExt = :distributions, LineCableModelsXLSXExt = :xlsx,
    LineCableModelsGmshExt = :fem, LineCableModelsMakieExt = :makie,
    LineCableModelsCairoMakieExt = :makie, LineCableModelsGLMakieExt = :makie,
    LineCableModelsWGLMakieExt = :makie)

# Each item has at least one kind tag.
const KINDS = (:unit, :integration, :extension, :visual, :quality, :aqua)

# A path relative to `directory`, with `/` separators on every platform.
relative(path, directory) = join(splitpath(relpath(path, directory)), "/")

# The files under `directory` that belong to the checkout: the files git tracks and the
# new files it would add, never ignored ones such as local captures under
# `test/fixtures/`. A local run and CI then read the same files. Outside a git work tree,
# every file under `directory`.
function repository_files(directory)
    isdir(directory) || return String[]
    listing = IOBuffer()
    command = `git -C $directory ls-files --cached --others --exclude-standard -z`
    if success(pipeline(command; stdout = listing, stderr = devnull))
        names = split(String(take!(listing)), '\0'; keepempty = false)
        return sort!(filter(isfile, [joinpath(directory, name) for name in names]))
    end
    return sort!([joinpath(path, name) for (path, _, names) in walkdir(directory) for name in names])
end

rank(owner::Symbol) = something(findfirst(==(owner), OWNERS), 0)
# The first owner tag among `tags`, or nothing.
owner_tag(tags) = (i = findfirst(in(OWNERS), tags); i === nothing ? nothing : tags[i])
latest(owners) = isempty(owners) ? nothing : OWNERS[maximum(rank, owners)]
earliest(owners) = isempty(owners) ? nothing : OWNERS[minimum(rank, owners)]

# The literal `include` paths of a source file, in order.
function included(path)
    found = String[]
    visit(x) = if Meta.isexpr(x, :call) && x.args[1] === :include && length(x.args) == 2 &&
                  x.args[2] isa String
        push!(found, normpath(joinpath(dirname(path), x.args[2])))
    elseif x isa Expr
        foreach(visit, x.args)
    end
    visit(Meta.parseall(read(path, String); filename = path))
    return found
end

# The module a file declares at top level, also behind a docstring.
function declared_module(path)
    for x in Meta.parseall(read(path, String); filename = path).args
        Meta.isexpr(x, :macrocall) && (x = x.args[end])
        Meta.isexpr(x, :module) && return x.args[2]
    end
    return nothing
end

# Each Julia file the package loads, relative to `repository`, with the owner of its
# load position. `src/LineCableModels.jl` includes the core modules in load order. A file
# it includes directly takes the owner of the last module included before it, or the
# first owner. Extensions take their own owners.
function loaded_owners(repository)
    found = Dict{String, Symbol}()
    function load!(path, owner)
        file = relative(path, repository)
        haskey(found, file) && error("$file is included more than once")
        found[file] = owner
        foreach(p -> load!(p, owner), included(path))
    end
    entry = joinpath(repository, "src", "LineCableModels.jl")
    found[relative(entry, repository)] = first(OWNERS)
    owner = first(OWNERS)
    for path in included(entry)
        name = declared_module(path)
        name !== nothing && haskey(MODULE_OWNERS, name) && (owner = MODULE_OWNERS[name])
        load!(path, owner)
    end
    for (name, extension) in pairs(EXTENSION_OWNERS)
        for path in (joinpath(repository, "ext", "$name.jl"),
            joinpath(repository, "ext", "$name", "$name.jl"))
            isfile(path) && load!(path, extension)
        end
    end
    return found
end

# The earliest owner loaded from `directory`, or nothing.
directory_owner(owners, directory) =
    earliest([owner for (file, owner) in owners if startswith(file, directory * "/")])

# The owner of a path under `src/` or `ext/`: the owner of its load position, or, for a
# file the package does not load, the earliest owner loaded from its nearest enclosing
# directory. Nothing for other paths.
function path_owner(owners, path)
    haskey(owners, path) && return owners[path]
    startswith(path, "src/") || startswith(path, "ext/") || return nothing
    directory = dirname(path)
    while !isempty(directory)
        owner = directory_owner(owners, directory)
        owner === nothing || return owner
        directory = dirname(directory)
    end
    return first(OWNERS)
end
