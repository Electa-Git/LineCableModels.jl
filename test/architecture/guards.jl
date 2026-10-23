# Each guard takes its inputs as arguments and returns its violations, keyed, with a count.
# `architecture.jl` applies it to the package and, as negative controls, to the probe
# packages under `probes/`.
module Guards

using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports,
                       check_all_explicit_imports_via_owners,
                       check_all_explicit_imports_are_public,
                       check_all_qualified_accesses_via_owners,
                       check_all_qualified_accesses_are_public,
                       check_no_self_qualified_accesses

# The module layout, in load order. Engine's children take Engine's position.
const ORDER = (:Commons, :Materials, :Earth, :DataModel, :Engine)
# Base and the standard libraries whose exported names share a user's scope with the names
# that the root module exports.
const DEPENDENCIES = (Base,)

# The core modules, the loaded extension modules and the repository
# directory that holds `src/` and `ext/`.
struct PackageTree
    root::Module
    core::Set{Module}
    extensions::Set{Module}
    modules::Set{Module}
    directory::String
end

function submodules(m::Module, found = Module[])
    push!(found, m)
    for name in names(m; all = true)
        isdefined(m, name) || continue
        value = getfield(m, name)
        value isa Module && value !== m && parentmodule(value) === m &&
            submodules(value, found)
    end
    return found
end

function PackageTree(root::Module, extensions, directory::AbstractString)
    core = Set(submodules(root))
    loaded = Set{Module}(Iterators.flatten(submodules(e) for e in extensions))
    return PackageTree(root, core, loaded, union(core, loaded), normpath(directory))
end

function is_descendant(m::Module, ancestor::Module)
    while true
        m === ancestor && return true
        parentmodule(m) === m && return false
        m = parentmodule(m)
    end
end

# For the root module, a descendant does not count as its own.
function within(m, owner::Module, root::Module)
    m isa Module && (m === owner || (owner !== root && is_descendant(m, owner)))
end

function module_name(tree::PackageTree, m::Module)
    m === tree.root && return string(nameof(m))
    path, prefix = fullname(m), fullname(tree.root)
    length(path) > length(prefix) && path[1:length(prefix)] == prefix &&
        return join(path[(length(prefix) + 1):end], ".")
    return join(path, ".")
end

inside(path, directory) = startswith(path, joinpath(directory, ""))

# Methods generated in `boot.jl` and similar report a relative file.
function source_path(tree::PackageTree, file)
    path = string(file)
    isabspath(path) || return nothing
    path = normpath(path)
    return inside(path, tree.directory) ? relative(path, tree.directory) : nothing
end
function source_name(tree::PackageTree, m::Method)
    something(source_path(tree, m.file), string(m.file))
end

function package_methods(tree::PackageTree)
    found = Method[]
    Base.visit(Core.methodtable) do m
        m.module in tree.modules && push!(found, m)
    end
    return found
end

function source_files(directory::AbstractString)
    files = String[]
    for base in ("src", "ext")
        append!(files, filter(endswith(".jl"), repository_files(joinpath(directory, base))))
    end
    return sort!(files)
end

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
    return sort!([joinpath(path, name) for (path, _, names) in walkdir(directory)
                  for name in names])
end

count!(found, key) = (found[key] = get(found, key, 0) + 1; found)

type_owner(T) = (T = Base.unwrap_unionall(T); T isa DataType ? parentmodule(T) : nothing)
value_owner(value::Module) = value
value_owner(value::Function) = parentmodule(typeof(value))
value_owner(value::Type) = type_owner(value)
value_owner(_) = nothing

# A constructor belongs to the owner of the constructed type.
function function_owner(m::Method)
    F = Base.unwrap_unionall(m.sig).parameters[1]
    F isa TypeVar && (F = F.ub)
    F = Base.unwrap_unionall(F)
    F isa DataType || return nothing
    if F.name === Type.body.name
        T = F.parameters[1]
        return type_owner(T isa TypeVar ? T.ub : T)
    end
    return parentmodule(F)
end

function mentioned!(found::Base.IdSet{DataType}, @nospecialize(type))
    if type isa TypeVar
        mentioned!(found, type.lb)
        mentioned!(found, type.ub)
    elseif type isa UnionAll
        mentioned!(found, type.var)
        mentioned!(found, type.body)
    elseif type isa Union
        mentioned!(found, type.a)
        mentioned!(found, type.b)
    elseif type isa Core.TypeofVararg
        isdefined(type, :T) && mentioned!(found, type.T)
    elseif type isa DataType && type ∉ found
        push!(found, type)
        for p in type.parameters
            mentioned!(found, p)
        end
    end
    return found
end

# Ownership. A core method that extends a function owned by another package
# module mentions a type owned by its own module or a descendant.
function ownership(methods, tree::PackageTree)
    found = Dict{String, Int}()
    for m in methods
        M = m.module
        M in tree.core || continue
        startswith(string(m.name), "#") && continue
        F = function_owner(m)
        (F isa Module && F in tree.modules) || continue
        within(F, M, tree.root) && continue
        types = Base.IdSet{DataType}()
        foreach(p -> mentioned!(types, p), Base.unwrap_unionall(m.sig).parameters[2:end])
        any(T -> within(parentmodule(T), M, tree.root), types) && continue
        count!(found,
            string(module_name(tree, M), " | ", module_name(tree, F), ".",
                m.name, " | ", source_name(tree, m)))
    end
    return found
end

# Each directory holding a module's `<ModuleName>.jl` is that module's home.
function home_directories(tree::PackageTree)
    files = source_files(tree.directory)
    homes = Dict{String, Vector{Module}}()
    for M in tree.modules
        name = string(nameof(M))
        declaration = Regex("(?m)^\\s*(?:bare)?module\\s+" * name * "\\b")
        candidates = filter(files) do file
            basename(file) == name * ".jl" && occursin(declaration, read(file, String))
        end
        # A module declared inline in another file has no home of its own.
        isempty(candidates) && continue
        length(candidates) == 1 ||
            error("Ambiguous home directory for $M: $(join(candidates, ", "))")
        push!(get!(homes, dirname(only(candidates)), Module[]), M)
    end
    return homes
end

function nearest_home(homes, file)
    nearest = nothing
    for directory in keys(homes)
        inside(file, directory) || continue
        (nearest === nothing || length(directory) > length(nearest)) &&
            (nearest = directory)
    end
    return nearest
end

# Placement. No core method is defined under `ext/`. The nearest home directory
# around a method's file belongs to its module or one of its ancestors.
function placement(methods, tree::PackageTree, homes)
    found = Dict{String, Int}()
    extensions = joinpath(tree.directory, "ext")
    for m in methods
        M = m.module
        path = source_path(tree, m.file)
        path === nothing && continue
        file = normpath(string(m.file))
        misplaced = M in tree.core && inside(file, extensions)
        nearest = nearest_home(homes, file)
        misplaced |= nearest === nothing ||
                     !any(owner -> is_descendant(M, owner), homes[nearest])
        misplaced && count!(found, string(path, " | ", module_name(tree, M)))
    end
    return found
end

function layer(tree::PackageTree, order, m::Module)
    m === tree.root && return nothing
    while parentmodule(m) !== tree.root
        parentmodule(m) === m && return nothing
        m = parentmodule(m)
    end
    index = findfirst(==(nameof(m)), order)
    index === nothing && error("Module $m is missing from the declared module order")
    return index
end

function lowered_code(m::Method)
    isdefined(m, :source) && m.source !== nothing &&
        !isdefined(m, :generator) ? Base.uncompressed_ir(m).code : nothing
end

# Every value that a lowered body names through a global reference or a
# `getproperty(module, :name)` chain of SSA values.
function referenced_values(code)
    resolved = Dict{Int, Any}()
    found = Any[]
    resolve(x::GlobalRef) = isdefined(x.mod, x.name) ? getglobal(x.mod, x.name) : nothing
    resolve(x::Core.SSAValue) = get(resolved, x.id, nothing)
    resolve(_) = nothing
    function references!(x)
        if x isa GlobalRef
            value = resolve(x)
            value === nothing || push!(found, value)
        elseif x isa Expr
            foreach(references!, x.args)
        end
    end
    for (i, statement) in enumerate(code)
        references!(statement)
        value = statement isa GlobalRef ? resolve(statement) : nothing
        if statement isa Expr && statement.head === :call && length(statement.args) == 3 &&
           resolve(statement.args[1]) === getproperty &&
           statement.args[3] isa QuoteNode
            parent = resolve(statement.args[2])
            name = statement.args[3].value
            parent isa Module && name isa Symbol && isdefined(parent, name) &&
                (value = getglobal(parent, name))
            value === nothing || push!(found, value)
        end
        value === nothing || (resolved[i] = value)
    end
    return found
end

# Direction. No method body of a core submodule references a value owned by a
# module in a later position of `order`.
function direction(methods, tree::PackageTree, order)
    found = Dict{String, Int}()
    for m in methods
        M = m.module
        M in tree.core || continue
        L = layer(tree, order, M)
        L === nothing && continue
        code = lowered_code(m)
        code === nothing && continue
        targets = Set{Module}()
        for value in referenced_values(code)
            owner = value_owner(value)
            (owner isa Module && owner in tree.core) || continue
            target = layer(tree, order, owner)
            target !== nothing && target > L && push!(targets, owner)
        end
        for owner in targets
            count!(found,
                string(module_name(tree, M), " -> ", module_name(tree, owner),
                    " | ", source_name(tree, m)))
        end
    end
    return found
end

# Names. Distinct functions or types, owned by different package modules,
# that share a name. The value is the number of modules that define the name.
function shared_names(tree::PackageTree)
    owners = Dict{Symbol, Vector{Module}}()
    for M in tree.modules, name in names(M; all = true)

        (startswith(string(name), "#") || name in (:eval, :include) ||
         name === nameof(M)) && continue
        isdefined(M, name) || continue
        value = getfield(M, name)
        (value isa Function || value isa Type) && value_owner(value) === M || continue
        push!(get!(owners, name, Module[]), M)
    end
    return Dict{String, Int}(string(name) => length(modules)
    for (name, modules) in owners if length(modules) > 1)
end

# Shadowing. A package-owned name that the root module exports, and so `using`
# brings into scope, and that a dependency exports for a different object. A name
# that only a submodule exports stays in that namespace.
function shadowing(tree::PackageTree, dependencies)
    found = Dict{String, Int}()
    for name in names(tree.root)
        Base.isexported(tree.root, name) && isdefined(tree.root, name) || continue
        value = getfield(tree.root, name)
        value isa Function || value isa Type || continue
        owner = value_owner(value)
        (owner isa Module && owner in tree.modules) || continue
        for dependency in dependencies
            Base.isexported(dependency, name) && isdefined(dependency, name) &&
            getfield(dependency, name) !== value || continue
            found[string(module_name(tree, owner), ".", name, " -> ",
                nameof(dependency))] = 1
        end
    end
    return found
end

# The file that defines an extension module: `ext/<name>.jl` or `ext/<name>/<name>.jl`.
function extension_file(tree::PackageTree, extension::Module)
    only(filter(isfile,
        [joinpath(tree.directory, "ext", string(nameof(extension), ".jl")),
            joinpath(tree.directory, "ext", string(nameof(extension)),
                string(nameof(extension), ".jl"))]))
end

# Explicit imports. Each ExplicitImports check runs on the root module, with its submodules,
# and on each extension. Imports and accesses inside the package are checked like any
# other, and a name used from a package module is public. Each reported name counts once,
# keyed by check, module and name. A check reports the first module that fails it.
function explicit_imports(tree::PackageTree, extensions)
    found = Dict{String, Int}()
    targets = [tree.root => joinpath(tree.directory, "src", string(nameof(tree.root), ".jl"))]
    append!(targets, [E => extension_file(tree, E) for E in extensions])
    package = (tree.root,)
    for (M, file) in targets
        checks = (
            "no_implicit_imports" => check_no_implicit_imports(
                M, file; skip = (
                    Base, Core), throw = false),
            "no_stale_explicit_imports" => check_no_stale_explicit_imports(M, file; throw = false),
            "all_explicit_imports_via_owners" => check_all_explicit_imports_via_owners(
                M, file;
                allow_internal_imports = false, throw = false),
            "all_explicit_imports_are_public" => check_all_explicit_imports_are_public(
                M, file; from = package,
                allow_internal_imports = false, throw = false),
            "all_qualified_accesses_via_owners" => check_all_qualified_accesses_via_owners(
                M, file;
                allow_internal_accesses = false, throw = false),
            "all_qualified_accesses_are_public" => check_all_qualified_accesses_are_public(
                M, file; from = package,
                allow_internal_accesses = false, throw = false),
            "no_self_qualified_accesses" => check_no_self_qualified_accesses(M, file; throw = false))
        for (check, exception) in checks
            exception === nothing && continue
            # Every check exception holds its module, then its rows.
            for row in getfield(exception, 2)
                found[string(check, " | ", nameof(exception.mod), " | ", row.name)] = 1
            end
        end
    end
    return found
end

end
