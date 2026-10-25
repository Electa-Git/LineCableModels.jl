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
using TOML: TOML
using JuliaSyntax: JuliaSyntax, SyntaxNode, @K_str, kind, children, numchildren, is_leaf,
                   is_prefix_call

# The module layout, in load order. Engine's children take Engine's position.
const ORDER = (:Commons, :Materials, :Earth, :DataModel, :Engine)
# Base and the standard libraries whose exported names share a user's scope with the names
# that the root module exports.
const DEPENDENCIES = (Base,)
# The home of Commons, the one module that defines `*Buffers` types.
const COMMONS_DIRECTORY = "src/commons"
const RESERVED_VERB = r"^_*(validate|check|require|assert|verify|ensure)_"
const COMMONS = :Commons
const COMMONS_TESTS = "test/commons"
# Each reserving name maps to its home and to the definition names it reserves outside that
# home. The physical constants and the matrix reductions enter it with their definitions.
const VOCABULARY = Dict{Symbol, Pair{String, Regex}}(
    :vacuum_permittivity => "src/engine/constants" => r"(?i)^_*(\w+_)?(vacuum_permittivity|eps(ilon)?_?[0₀]|[εϵ]_?[0₀])$",
    :vacuum_permeability => "src/engine/constants" => r"(?i)^_*(\w+_)?(vacuum_permeability|mu_?[0₀]|[μµ]_?[0₀])$")
const SHINGLE = 8
const CLONE_TOKENS = 30
const CLONE_SHARE = 0.75

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

# Names. Distinct functions or types, owned by different package modules, that share a name.
# Types share a name when each descends from one abstract type that Commons owns: the
# namespace specializes one concept. The value is the number of modules that define the name.
function shared_names(tree::PackageTree)
    owners = Dict{Symbol, Vector{Module}}()
    values = Dict{Symbol, Vector{Any}}()
    for M in tree.modules, name in names(M; all = true)

        (startswith(string(name), "#") || name in (:eval, :include) ||
         name === nameof(M)) && continue
        isdefined(M, name) || continue
        value = getfield(M, name)
        (value isa Function || value isa Type) && value_owner(value) === M || continue
        push!(get!(owners, name, Module[]), M)
        push!(get!(values, name, Any[]), value)
    end
    C = commons_module(tree)
    return Dict{String, Int}(string(name) => length(modules)
    for (name, modules) in owners
    if length(modules) > 1 && !family(values[name], tree, C))
end

# Whether all `values` are types that descend from one abstract type that Commons owns.
function family(values, tree::PackageTree, C)
    C === nothing && return false
    all(value -> value isa Type, values) || return false
    function ancestors(value)
        found = Any[]
        S = supertype(Base.unwrap_unionall(value))
        while S !== Any
            within(S.name.module, C, tree.root) && push!(found, S.name.wrapper)
            S = supertype(S)
        end
        return found
    end
    return !isempty(intersect(map(ancestors, values)...))
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
            # A check exception stores its module, then its rows.
            for row in getfield(exception, 2)
                found[string(check, " | ", nameof(exception.mod), " | ", row.name)] = 1
            end
        end
    end
    return found
end

# Source spelling only. No binding is resolved.
function syntax_name(node)
    kind(node) in (K"Identifier", K"MacroName") && return (node.val,)
    if kind(node) == K"." && numchildren(node) == 2
        member = node[2]
        kind(member) == K"quote" && numchildren(member) == 1 && (member = member[1])
        left, right = syntax_name(node[1]), syntax_name(member)
        !isempty(left) && length(right) == 1 && return (left..., right...)
    end
    return ()
end

# Calls `f` on every node. A `false` result skips the children of that node.
function visit(f, node)
    f(node) === false && return
    is_leaf(node) || foreach(child -> visit(f, child), children(node))
    return
end

parse_source(file) = JuliaSyntax.parseall(SyntaxNode, read(file, String); filename = file)
last_child(node) = children(node)[end]

# The call of a method definition, without `where` and return annotations.
function signature_call(node)
    kind(node) == K"function" && numchildren(node) == 2 || return nothing
    signature = node[1]
    while kind(signature) in (K"where", K"::") && numchildren(signature) == 2
        signature = signature[1]
    end
    return kind(signature) == K"call" && is_prefix_call(signature) ? signature : nothing
end

function definition_name(node)
    kind(node) == K"function" || return nothing
    call = signature_call(node)
    head = call !== nothing ? call[1] : numchildren(node) == 1 ? node[1] : nothing
    head === nothing && return nothing
    name = syntax_name(head)
    return isempty(name) ? nothing : last(name)
end

function first_argument(call)
    arguments = [a for a in children(call)[2:end] if kind(a) != K"parameters"]
    isempty(arguments) && return nothing
    argument = arguments[1]
    kind(argument) == K"=" && numchildren(argument) == 2 && (argument = argument[1])
    kind(argument) == K"::" && numchildren(argument) == 2 && (argument = argument[1])
    return kind(argument) == K"Identifier" ? argument.val : nothing
end

const THROWS = ((:throw,), (:rethrow,), (:error,), (:Base, :throw), (:Core, :throw),
    (:Base, :rethrow), (:Base, :error))

# Whether an expression in tail position yields `name` or throws.
function returns_subject(node, name)
    k = kind(node)
    k == K"Identifier" && return node.val === name
    k in (K"block", K"let", K"parens") &&
        return numchildren(node) > 0 && returns_subject(last_child(node), name)
    k == K"return" && return numchildren(node) == 1 && returns_subject(node[1], name)
    k in (K"if", K"elseif", K"?") && return numchildren(node) == 3 &&
           returns_subject(node[2], name) && returns_subject(node[3], name)
    if k == K"try"
        clauses = Dict(kind(c) => c for c in children(node)[2:end])
        value = haskey(clauses, K"else") ? last_child(clauses[K"else"]) : node[1]
        return returns_subject(value, name) && (!haskey(clauses, K"catch") ||
                returns_subject(last_child(clauses[K"catch"]), name))
    end
    return k == K"call" && is_prefix_call(node) && syntax_name(node[1]) in THROWS
end

# Every `return` of this method, outside nested functions and quotations.
function early_returns_subject(node, name)
    is_leaf(node) && return true
    for child in children(node)
        kind(child) in (K"function", K"->", K"do", K"macro", K"quote") && continue
        kind(child) == K"return" && !returns_subject(child, name) && return false
        early_returns_subject(child, name) || return false
    end
    return true
end

# The occurrences of `name` that a body returns: its tail values and the values of
# its `return` statements outside nested functions.
function returned!(found, node, name)
    k = kind(node)
    if k == K"Identifier"
        node.val === name && push!(found, node)
    elseif k in (K"block", K"let", K"parens")
        numchildren(node) > 0 && returned!(found, last_child(node), name)
    elseif k == K"return"
        numchildren(node) == 1 && returned!(found, node[1], name)
    elseif k in (K"if", K"elseif", K"?")
        foreach(child -> returned!(found, child, name), children(node)[2:end])
    elseif k == K"try"
        foreach(child -> returned!(found, child, name), children(node))
    end
    return found
end

function early_returned!(found, node, name)
    is_leaf(node) && return found
    for child in children(node)
        kind(child) in (K"function", K"->", K"do", K"macro", K"quote") && continue
        kind(child) == K"return" && returned!(found, child, name)
        early_returned!(found, child, name)
    end
    return found
end

# Whether a body uses `name` other than by returning it. A body that only returns
# `name` is an identity method: its consumer does not restrict the subject.
function uses_subject(body, name)
    kind(body) == K"Identifier" && return true
    kind(body) == K"block" && numchildren(body) == 1 &&
        kind(body[1]) in (K"Identifier", K"return") && return true
    found = early_returned!(returned!(Base.IdSet{SyntaxNode}(), body, name), body, name)
    used = false
    visit(body) do node
        kind(node) == K"Identifier" && node.val === name && node ∉ found && (used = true)
        return !used
    end
    return used
end

function validate_returns_subject(node)
    name = first_argument(signature_call(node))
    name === nothing && return false
    body = node[2]
    return returns_subject(body, name) && early_returns_subject(body, name) &&
           uses_subject(body, name)
end

# validate returns its subject. Every `validate` method names its first positional argument,
# returns it on every path, or throws, and uses it in its work unless it only returns it.
function validate_returns(files, directory)
    found = Dict{String, Int}()
    for file in files
        visit(parse_source(file)) do node
            signature_call(node) === nothing && return true
            definition_name(node) === :validate || return true
            validate_returns_subject(node) || count!(found, relative(file, directory))
            return true
        end
    end
    return found
end

# Reserved verbs. No function name spells an input check with a reserved verb prefix.
function reserved_verbs(files, directory)
    found = Dict{String, Int}()
    for file in files
        visit(parse_source(file)) do node
            name = definition_name(node)
            name isa Symbol && occursin(RESERVED_VERB, string(name)) &&
                count!(found, string(relative(file, directory), " | ", name))
            return true
        end
    end
    return found
end

function symbol_literal(node)
    kind(node) == K"quote" && numchildren(node) == 1 &&
        kind(node[1]) == K"Identifier"
end
function symbol_operand(node)
    symbol_literal(node) ||
        (kind(node) in (K"tuple", K"vect") && numchildren(node) > 0 &&
         all(symbol_literal, children(node)))
end

const NOTHING_TYPES = ((:Nothing,), (:Base, :Nothing), (:Core, :Nothing))

# Whether a `::` node types a positional or keyword argument, with or without a default.
function typed_argument(node)
    parent = node.parent
    parent !== nothing && kind(parent) == K"=" && parent[1] === node &&
        (node = parent; parent = node.parent)
    parent === nothing && return false
    kind(parent) == K"->" && return parent[1] === node
    return kind(parent) in (K"call", K"parameters", K"tuple")
end

const NOT_ASSIGNED = (K"call", K"dotcall", K"parameters", K"tuple", K"vect", K"braces",
    K"curly", K"ref")

function assigned_names(target)
    kind(target) == K"Identifier" && return [target.val]
    kind(target) == K"::" && numchildren(target) == 2 && return assigned_names(target[1])
    kind(target) == K"tuple" && return reduce(vcat, map(assigned_names, children(target));
        init = Symbol[])
    return Symbol[]
end

const STORAGE_NAMES = (:work, :scratch, :storage, :cache)
const STORAGE_TYPE = r"(Buffer|Storage|Scratch|Cache)$"
const STORAGE_FUNCTION = r"_(workspace|storage|buffers)!?$"

# The name of a struct, abstract or primitive type definition.
function declared_name(node)
    head = node[1]
    while kind(head) in (K"<:", K"curly") && numchildren(head) > 0
        head = head[1]
    end
    return kind(head) == K"Identifier" ? head.val : nothing
end

# The names that a parameter, an assignment target or a loop variable binds.
function bound_names(node)
    k = kind(node)
    k == K"Identifier" && return [node.val]
    (k == K"::" && numchildren(node) == 2 || k in (K"=", K"...")) &&
        return bound_names(node[1])
    k in (K"tuple", K"parameters") &&
        return reduce(vcat, map(bound_names, children(node)); init = Symbol[])
    return Symbol[]
end

# The field that one entry of a struct body declares, or nothing.
function field_name(node)
    k = kind(node)
    k == K"Identifier" && return node.val
    k == K"doc" && return field_name(last_child(node))
    k in (K"::", K"const", K"=") && numchildren(node) > 0 && return field_name(node[1])
    return nothing
end

# The keys of a named tuple literal.
function tuple_keys(node)
    keys = Symbol[]
    for child in children(node)
        if kind(child) == K"=" && kind(child[1]) == K"Identifier"
            push!(keys, child[1].val)
        elseif kind(child) == K"parameters"
            for entry in children(child)
                kind(entry) == K"Identifier" && push!(keys, entry.val)
                kind(entry) == K"=" && kind(entry[1]) == K"Identifier" &&
                    push!(keys, entry[1].val)
            end
        end
    end
    return keys
end

function in_struct_body(node)
    node.parent !== nothing && kind(node.parent) == K"block" &&
        node.parent.parent !== nothing && kind(node.parent.parent) == K"struct"
end

# Storage vocabulary. No function name ends in `_workspace` or `_storage`, and only
# `initialize_buffers` ends in `_buffers`. Only Commons defines a `*Buffers` struct, and no
# type name ends in `Buffer`, `Storage`, `Scratch` or `Cache`. No field, named tuple key,
# parameter, local or loop variable takes the name `work`, `scratch`, `storage` or `cache`.
function storage_vocabulary(files, directory)
    found = Dict{String, Int}()
    commons = joinpath(directory, COMMONS_DIRECTORY)
    for file in files
        path = relative(file, directory)
        entry(name,
            n = 1) = n > 0 &&
                     (found[string(path, " | ", name)] = get(found, string(path, " | ", name), 0) +
                                                         n)
        names!(names) = foreach(name -> name in STORAGE_NAMES && entry(name), names)
        visit(parse_source(file)) do node
            k = kind(node)
            if k == K"struct"
                name = declared_name(node)
                name === nothing && return true
                text = string(name)
                fields = filter(!isnothing, map(field_name, children(node[2])))
                names!(fields)
                endswith(text, "Buffers") && !inside(file, commons) && entry(name)
                occursin(STORAGE_TYPE, text) && entry(name)
            elseif k in (K"abstract", K"primitive")
                name = declared_name(node)
                name !== nothing && occursin(STORAGE_TYPE, string(name)) && entry(name)
            elseif k == K"function"
                name = definition_name(node)
                name isa Symbol && name !== :initialize_buffers &&
                    occursin(STORAGE_FUNCTION, string(name)) && entry(name)
                call = signature_call(node)
                call === nothing || names!(reduce(
                    vcat, map(bound_names,
                        children(call)[2:end]); init = Symbol[]))
            elseif k in (K"->", K"do") && numchildren(node) > 0
                names!(bound_names(node[1]))
            elseif k == K"in" && node.parent !== nothing &&
                   kind(node.parent) == K"iteration"
                names!(bound_names(node[1]))
            elseif k == K"tuple" &&
                   !(node.parent !== nothing && kind(node.parent) in (K"->", K"do"))
                names!(tuple_keys(node))
            elseif k == K"=" && numchildren(node) == 2 && !in_struct_body(node) &&
                   !(node.parent !== nothing &&
                     kind(node.parent) in (NOT_ASSIGNED..., K"const"))
                names!(assigned_names(node[1]))
            elseif k == K"local"
                for child in children(node)
                    kind(child) == K"Identifier" && names!([child.val])
                end
            end
            return true
        end
    end
    return found
end

# Repeated access names. Per file, the chains whose base binding has the name of their
# first field, as in `x.x` or `x.:x`. The count leaves out a field repeated inside a
# field, as in `a.x.x`.
function repeated_access_names(files, directory)
    found = Dict{String, Int}()
    for file in files
        visit(parse_source(file)) do node
            kind(node) == K"." && numchildren(node) == 2 &&
            kind(node[1]) == K"Identifier" || return true
            member = node[2]
            kind(member) == K"quote" && numchildren(member) == 1 && (member = member[1])
            kind(member) == K"Identifier" && member.val === node[1].val &&
                count!(found, relative(file, directory))
            return true
        end
    end
    return found
end

const APPLICABLE = ((:applicable,), (:Base, :applicable), (:Core, :applicable))
const EVAL = ((Symbol("@eval"),), (:Base, Symbol("@eval")), (:Core, Symbol("@eval")))
const ISA = ((:isa,), (:Base, :isa), (:Core, :isa))
# Comparisons, negated ones included, so negating a switch does not hide it.
const COMPARISONS = ((:(===),), (:(==),), (:in,), (:∈,), (:(!==),), (:(!=),), (:∉,))
const EQUALITIES = ((:(===),), (:(==),), (:(!==),), (:(!=),))
# The parents of a `<:` that declares a bound or a supertype instead of testing one.
const DECLARATIONS = (K"curly", K"braces", K"where", K"struct", K"abstract", K"primitive")
# The calls whose predicate selects from a collection.
const SELECTIONS = (:filter, :filter!, :findall, :findfirst, :findlast)

# The operator and the operands of a call with two arguments, infix or prefix.
function binary(node)
    kind(node) == K"call" && numchildren(node) == 3 || return nothing
    return is_prefix_call(node) ? (node[1], node[2], node[3]) : (node[2], node[1], node[3])
end

# Whether `node` calls a function whose last name is in `names`, qualified or not.
function calls(node, names)
    kind(node) == K"call" && is_prefix_call(node) && numchildren(node) > 0 || return false
    name = syntax_name(node[1])
    return !isempty(name) && last(name) in names
end

function symbol_comparison(node)
    operands = binary(node)
    operands === nothing && return false
    operator, left, right = operands
    return syntax_name(operator) in COMPARISONS &&
           (symbol_operand(left) || symbol_operand(right))
end

# A type test: `isa`, a `<:` used as a value, or an equality with a `typeof` call.
function type_test(node)
    kind(node) == K"<:" && return kind(node.parent) ∉ DECLARATIONS
    operands = binary(node)
    operands === nothing && return false
    operator, left, right = operands
    syntax_name(operator) in ISA && return true
    return syntax_name(operator) in EQUALITIES &&
           (calls(left, (:typeof,)) || calls(right, (:typeof,)))
end

# Whether a type test is admitted: inside a `catch` clause, in the condition of a
# comprehension or generator, or in the predicate of a selection call, given as an
# anonymous function or a `do` block.
function admitted(node)
    child, parent = node, node.parent
    while parent !== nothing
        k = kind(parent)
        k == K"catch" && return true
        k == K"filter" && child === last_child(parent) && return true
        k in (K"->", K"do") && calls(parent.parent, SELECTIONS) && return true
        child, parent = parent, parent.parent
    end
    return false
end

# Symbol switches and probes. Per file and form: calls to `applicable`, `@eval` calls,
# comparisons against a Symbol literal, and type tests that are not admitted.
function symbol_switches(files, directory)
    found = Dict{String, Int}()
    for file in files
        path = relative(file, directory)
        visit(parse_source(file)) do node
            calls(node, (:applicable,)) && syntax_name(node[1]) in APPLICABLE &&
                count!(found, string(path, " | applicable"))
            kind(node) == K"macrocall" && numchildren(node) > 0 &&
                syntax_name(node[1]) in EVAL && count!(found, string(path, " | eval"))
            symbol_comparison(node) && count!(found, string(path, " | symbol"))
            type_test(node) && !admitted(node) &&
                count!(found, string(path, " | type test"))
            return true
        end
    end
    return found
end

const MISSING_TYPES = ((:Missing,), (:Base, :Missing), (:Core, :Missing))
const ABSENT_TYPES = (NOTHING_TYPES..., MISSING_TYPES...)
const ABSENT_VALUES = ((:nothing,), (:missing,))
# The Julia functions that return `nothing` by their protocol.
const PROTOCOL = (:iterate, :findfirst, :findnext, :findlast, :findprev, :match, :tryparse)

# The nearest enclosing function, anonymous function or `do` block, or the file.
function scope(node)
    while node.parent !== nothing
        node = node.parent
        kind(node) in (K"function", K"->", K"do") && return node
    end
    return node
end

# The right-hand sides of the assignments to `name` within `scope`.
function assignments(scope, name)
    found = SyntaxNode[]
    visit(scope) do node
        kind(node) == K"=" && numchildren(node) == 2 && kind(node[1]) == K"Identifier" &&
            node[1].val === name && push!(found, node[2])
        return true
    end
    return found
end

# Whether `node` is the result of a protocol call: the call, an assignment from it, or a
# local that every assignment in its scope takes from such a call.
function protocol_result(node)
    calls(node, PROTOCOL) && return true
    kind(node) == K"=" && numchildren(node) == 2 && return calls(node[2], PROTOCOL)
    kind(node) == K"Identifier" || return false
    assigned = assignments(scope(node), node.val)
    return !isempty(assigned) && all(value -> calls(value, PROTOCOL), assigned)
end

# The call whose callee is `node`, qualified or not, or nothing.
function callee_call(node)
    parent = node.parent
    while parent !== nothing && kind(parent) in (K".", K"quote")
        node, parent = parent, parent.parent
    end
    parent !== nothing && kind(parent) == K"call" && is_prefix_call(parent) &&
    parent[1] === node || return nothing
    return parent
end

# Decisions by absence. Per file and form: comparisons with `nothing` or `missing`;
# references to `isnothing`, `ismissing` and `something`, called or passed; arguments typed
# `Nothing` or `Missing`; and `Nothing` or `Missing` as a type parameter. A comparison with
# `nothing`, an `isnothing` or a `something` on the result of a protocol call is exempt.
function absences(files, directory)
    found = Dict{String, Int}()
    for file in files
        path = relative(file, directory)
        visit(parse_source(file)) do node
            k = kind(node)
            operands = binary(node)
            if operands !== nothing && syntax_name(operands[1]) in EQUALITIES
                _, left, right = operands
                for (value, other) in ((left, right), (right, left))
                    name = syntax_name(value)
                    name in ABSENT_VALUES || continue
                    name == (:nothing,) && protocol_result(other) && continue
                    count!(found, string(path, " | comparison"))
                end
            elseif k == K"Identifier" && node.val in (:isnothing, :ismissing, :something)
                call = callee_call(node)
                node.val !== :ismissing && call !== nothing && numchildren(call) > 1 &&
                    protocol_result(call[2]) && return true
                count!(found, string(path, " | ", node.val))
            elseif k == K"::" && syntax_name(last_child(node)) in ABSENT_TYPES &&
                   typed_argument(node)
                count!(found, string(path, " | argument"))
            elseif k == K"curly"
                for parameter in children(node)[2:end]
                    syntax_name(parameter) in ABSENT_TYPES &&
                        count!(found, string(path, " | parameter"))
                end
            end
            return true
        end
    end
    return found
end

const UNHYGIENIC = (:tempdir, :mktempdir, :tempname, :homedir)

# File hygiene. Per file and name, the references to `tempdir`, `mktempdir`, `tempname`
# and `homedir`, called or passed.
function file_hygiene(files, directory)
    found = Dict{String, Int}()
    for file in files
        visit(parse_source(file)) do node
            kind(node) == K"Identifier" && node.val in UNHYGIENIC &&
                count!(found, string(relative(file, directory), " | ", node.val))
            return true
        end
    end
    return found
end

function upper(T)
    T isa TypeVar ? upper(T.ub) :
    T isa Core.TypeofVararg ? (isdefined(T, :T) ? upper(T.T) : Any) : T
end

type_name(T) = (T = Base.unwrap_unionall(T); T isa DataType ? nameof(T) : nothing)

# The functions, types and constants that a module defines, by name. A const binding of a
# function that another module owns is a name the module defines.
function definitions(M::Module)
    found = Dict{Symbol, Any}()
    for name in names(M; all = true)
        (startswith(string(name), "#") || name in (:eval, :include) ||
         name === nameof(M)) && continue
        isdefined(M, name) || continue
        value = getfield(M, name)
        value isa Module && continue
        owner = value isa Function ? value_owner(value) :
                value isa Type && type_name(value) === name ? type_owner(value) :
                isconst(M, name) ? Base.binding_module(M, name) : nothing
        alias = value isa Function && isconst(M, name) && Base.binding_module(M, name) === M
        (owner === M || alias) && (found[name] = value)
    end
    return found
end

# The function or type whose method `m` is. A keyword sorter belongs to its
# function, and a method of a nested function to the type of that function.
function method_subject(m::Method)
    signature = Base.unwrap_unionall(m.sig).parameters
    F = Base.unwrap_unionall(signature[1])
    F isa DataType || return nothing
    F === typeof(Core.kwcall) && length(signature) >= 3 &&
        (F = Base.unwrap_unionall(signature[3]); F isa DataType || return nothing)
    if F.name === Type.body.name
        T = F.parameters[1]
        T = Base.unwrap_unionall(T isa TypeVar ? T.ub : T)
        return T isa DataType ? T.name.wrapper : nothing
    end
    return isdefined(F, :instance) ? F.instance : F.name.wrapper
end

# The root module, a top-level core submodule or a top-level extension module.
function top_owner(tree::PackageTree, m::Module)
    while m !== tree.root && parentmodule(m) !== m && parentmodule(m) !== tree.root &&
              parentmodule(m) in tree.modules
        m = parentmodule(m)
    end
    return m
end

# The texts of the checkout's `.jl` files under `directory`, by path, except those under
# `excluded` subdirectories.
function test_sources(directory; excluded = ())
    texts = Dict{String, String}()
    for file in repository_files(directory)
        endswith(file, ".jl") &&
            !any(x -> inside(file, joinpath(directory, x)), excluded) &&
            (texts[file] = read(file, String))
    end
    return texts
end

# The paths of the texts that name `name`, sorted.
function naming(texts, name)
    pattern = Regex("(?<![\\w!@])\\Q" * string(name) * "\\E(?![\\w!])")
    return sort!([path for (path, text) in texts if occursin(pattern, text)])
end

function commons_module(tree::PackageTree)
    (C = isdefined(tree.root, COMMONS) ?
         getfield(tree.root, COMMONS) : nothing;
        C isa Module ? C : nothing)
end

# The grammar and the user API of Commons. The grammar is the Commons functions with a method
# defined outside Commons and the Commons abstract types with a subtype defined outside
# Commons. The user API is the Commons names that the root module exports or declares public.
function commons_roles(methods, tree::PackageTree, C::Module, defined)
    functions = IdDict{Any, Symbol}(value => name
    for (name, value) in defined
    if value isa Function)
    abstracts = IdDict{Any, Symbol}(value => name
    for (name, value) in defined
    if value isa Type && isabstracttype(Base.unwrap_unionall(value)))
    grammar = Set{Symbol}()
    for m in methods
        within(m.module, C, tree.root) && continue
        name = get(functions, method_subject(m), nothing)
        name === nothing || push!(grammar, name)
    end
    for M in tree.modules
        within(M, C, tree.root) && continue
        for value in values(definitions(M))
            value isa Type || continue
            S = supertype(Base.unwrap_unionall(value))
            while S !== Any
                name = get(abstracts, S.name.wrapper, nothing)
                name === nothing || push!(grammar, name)
                S = supertype(S)
            end
        end
    end
    api = Set(name
    for (name, value) in defined
    if Base.ispublic(tree.root, name) &&
       isdefined(tree.root, name) && getfield(tree.root, name) === value)
    return grammar, api
end

# The owners outside Commons whose methods use each Commons name, directly or through
# the `public` Commons definitions. Macros are used where source code calls them.
function commons_owners(methods, tree::PackageTree, C::Module, defined, public, homes)
    name_of = IdDict{Any, Symbol}(value => name for (name, value) in defined)
    # Owners of every Commons value, nested functions and keyword bodies included.
    owners = IdDict{Any, Set{Module}}(value => Set{Module}() for value in values(defined))
    # The Commons values that each Commons method subject references.
    references = IdDict{Any, Base.IdSet{Any}}()
    unnamed(value) = !haskey(name_of, value) && (value isa Function &&
                      within(parentmodule(typeof(value)), C, tree.root) ||
                      value isa Type && value <: Function &&
                      within(parentmodule(value), C, tree.root))
    function use!(M::Module, @nospecialize(subject), @nospecialize(value))
        haskey(name_of, value) || unnamed(value) || return
        if within(M, C, tree.root)
            subject === nothing || subject === value ||
                push!(get!(Base.IdSet{Any}, references, subject), value)
        elseif haskey(name_of, value)
            push!(owners[value], top_owner(tree, M))
        end
    end
    # Loops instead of anonymous functions: an anonymous function compiles once
    # per type it receives.
    for m in methods
        subject = method_subject(m)
        use!(m.module, subject, subject)
        types = Base.IdSet{DataType}()
        for p in Base.unwrap_unionall(m.sig).parameters[2:end]
            mentioned!(types, p)
        end
        for T in types
            use!(m.module, subject, T.name.wrapper)
        end
        code = lowered_code(m)
        code === nothing && continue
        for value in referenced_values(code)
            use!(m.module, subject, value)
        end
    end
    # Supertypes and field types of every package type.
    for M in tree.modules, value in values(definitions(M))

        T = value isa Type ? Base.unwrap_unionall(value) : nothing
        T isa DataType || continue
        related = Base.IdSet{DataType}()
        S = supertype(T)
        while S !== Any
            push!(related, S)
            S = supertype(S)
        end
        if isstructtype(T)
            for F in Base.datatype_fieldtypes(T)
                mentioned!(related, F)
            end
        end
        for S in related
            use!(M, value, S.name.wrapper)
        end
    end
    macros = Set(name for name in keys(defined) if startswith(string(name), "@"))
    for file in source_files(tree.directory)
        isempty(macros) && break
        home = nearest_home(homes, file)
        home === nothing && continue
        M = first(homes[home])
        visit(parse_source(file)) do node
            kind(node) == K"macrocall" && numchildren(node) > 0 || return true
            name = syntax_name(node[1])
            !isempty(name) && last(name) in macros && use!(M, nothing, defined[last(name)])
            return true
        end
    end
    # Public definitions pass their owners to what they use, also through
    # their nested functions and keyword bodies. Private definitions pass none.
    changed = true
    while changed
        changed = false
        for (subject, targets) in references
            haskey(name_of, subject) && name_of[subject] ∉ public && continue
            from = get!(Set{Module}, owners, subject)
            for target in targets
                to = get!(Set{Module}, owners, target)
                issubset(from, to) && continue
                union!(to, from)
                changed = true
            end
        end
    end
    return Dict(name => owners[value] for (name, value) in defined)
end

# Commons admission. Each Commons name is public, documented and named in a Commons test.
# Unless it is grammar or user API, methods of at least two owners outside Commons use it,
# directly or through public Commons definitions. An extension is not an owner.
function commons_admission(methods, tree::PackageTree, homes)
    found = Dict{String, Int}()
    C = commons_module(tree)
    C === nothing && return found
    defined = definitions(C)
    grammar, api = commons_roles(methods, tree, C, defined)
    public = Set(name for name in keys(defined) if Base.ispublic(C, name))
    users = commons_owners(methods, tree, C, defined, public, homes)
    tests = test_sources(joinpath(tree.directory, COMMONS_TESTS))
    for name in keys(defined)
        owners = count(in(tree.core), users[name])
        for (criterion, met) in (("public", name in public),
            ("docstring", Base.Docs.hasdoc(C, name)),
            ("owners", name in grammar || name in api || owners >= 2),
            ("tests", !isempty(naming(tests, name))))
            met || (found[string(name, " | ", criterion)] = 1)
        end
    end
    return found
end

function reserving(vocabulary, name)
    for (owner, (_, pattern)) in vocabulary
        occursin(pattern, string(name)) && return owner
    end
    return nothing
end

# Reserved vocabulary. Outside the home of a reserving name, no function, constant or
# assigned local takes a name that it reserves. A local assigned from a call to the
# reserving definition is exempt.
function reserved_vocabulary(files, directory, vocabulary)
    found = Dict{String, Int}()
    for file in files
        path = relative(file, directory)
        home(owner) = inside(file, joinpath(directory, first(vocabulary[owner])))
        visit(parse_source(file)) do node
            if kind(node) == K"function"
                name = definition_name(node)
                owner = name isa Symbol ? reserving(vocabulary, name) : nothing
                owner === nothing || home(owner) ||
                    count!(found, string(path, " | ", name))
            elseif kind(node) == K"=" && numchildren(node) == 2 &&
                   !(node.parent !== nothing && kind(node.parent) in NOT_ASSIGNED)
                for name in assigned_names(node[1])
                    owner = reserving(vocabulary, name)
                    owner === nothing || home(owner) || calls(node[2], (owner,)) ||
                        count!(found, string(path, " | ", name))
                end
            end
            return true
        end
    end
    return found
end

function any_node(predicate, node)
    predicate(node) ||
        (!is_leaf(node) && any(child -> any_node(predicate, child), children(node)))
end

# The innermost statement around a node.
function statement(node)
    while node.parent !== nothing && !(kind(node.parent) in (K"block", K"toplevel"))
        node = node.parent
    end
    return node
end

is_pi(node) = kind(node) == K"Identifier" && node.val in (:π, :pi)
function is_number(node)
    kind(node) in (K"Integer", K"Float", K"Float32", K"HexInt", K"OctInt",
        K"BinInt")
end

# `base ^ -7` with a literal 10 in `base`.
function minus_seventh_power(node)
    kind(node) == K"call" && numchildren(node) == 3 || return false
    operator, base,
    exponent = is_prefix_call(node) ?
               (node[1], node[2], node[3]) : (node[2], node[1], node[3])
    return syntax_name(operator) == (:^,) && kind(exponent) == K"Integer" &&
           exponent.val == -7 && any_node(n -> kind(n) == K"Integer" && n.val == 10, base)
end

# Whether `node` is a numeric fingerprint of a vacuum constant: the digits of ε₀ or c₀, or
# 10⁻⁷ in a statement that names π.
function fingerprint(node)
    if is_number(node)
        digits = filter(isdigit, JuliaSyntax.sourcetext(node))
        (occursin("8854187", digits) || occursin("299792458", digits)) && return true
        value = node.val
        value isa AbstractFloat && value == oftype(value, 1e-7) || return false
    else
        minus_seventh_power(node) || return false
    end
    return any_node(is_pi, statement(node))
end

# Literal fingerprints. The numeric fingerprints of the vacuum constants appear in one file,
# the home of the physical constants. When several files hold them, each file counts them.
function fingerprints(files, directory)
    found = Dict{String, Int}()
    for file in files
        path = relative(file, directory)
        visit(parse_source(file)) do node
            fingerprint(node) && count!(found, path)
            return true
        end
    end
    return length(found) > 1 ? found : Dict{String, Int}()
end

# A statement `x || throw(...)` or `x && throw(...)`, chains included.
function throw_guard(node)
    kind(node) in (K"||", K"&&") || return false
    while kind(node) in (K"||", K"&&") && numchildren(node) > 0
        node = last_child(node)
    end
    return kind(node) == K"call" && is_prefix_call(node) && syntax_name(node[1]) in THROWS
end

# A body's tokens, without throw guards, with each identifier replaced by
# its role.
function body_tokens(body)
    bytes = Vector{UInt8}(JuliaSyntax.sourcetext(body))
    offset = first(JuliaSyntax.byte_range(body)) - 1
    visit(body) do node
        statement = node === body ||
                    (node.parent !== nothing &&
                     kind(node.parent) in (K"block", K"toplevel"))
        statement && throw_guard(node) || return true
        bytes[JuliaSyntax.byte_range(node) .- offset] .= UInt8(' ')
        return false
    end
    text = String(bytes)
    raw = [(kind(t), JuliaSyntax.untokenize(t, text))
           for t in JuliaSyntax.tokenize(text)
           if !(kind(t) in (K"Whitespace", K"NewlineWs", K"Comment"))]
    tokens = String[]
    for (i, (k, token)) in enumerate(raw)
        if k == K"Identifier" && !Base.isoperator(Symbol(token))
            previous = i > 1 ? raw[i - 1][2] : ""
            next = i < length(raw) ? raw[i + 1][2] : ""
            token = previous == "." ? "FIELD" : next == "(" ? "CALL" : "NAME"
        end
        push!(tokens, token)
    end
    return tokens
end

function shingles(tokens)
    Set(join(view(tokens, i:(i + SHINGLE - 1)), '\x1f')
    for i in 1:(length(tokens) - SHINGLE + 1))
end

function function_bodies(trees)
    found = Tuple{String, Any, Any}[]
    for (path, tree) in trees
        visit(tree) do node
            kind(node) == K"function" && numchildren(node) == 2 &&
                push!(found, (path, node, node[2]))
            return true
        end
    end
    return found
end

# Clones of Commons bodies. No function body contains `CLONE_SHARE` of the token shingles of a
# Commons public function body of at least `CLONE_TOKENS` tokens. Hooks are
# not references.
function clones(files, methods, tree::PackageTree; share = CLONE_SHARE)
    found = Dict{String, Int}()
    C = commons_module(tree)
    C === nothing && return found
    hooks, _ = commons_roles(methods, tree, C, definitions(C))
    public = Set(name
    for name in names(C; all = true)
    if Base.ispublic(C, name) && name ∉ hooks)
    trees = [(relative(file, tree.directory), parse_source(file)) for file in files]
    bodies = [(path, node, shingles(body_tokens(body)))
              for (path, node, body) in function_bodies(trees)]
    commons = COMMONS_DIRECTORY * "/"
    for (path, reference, body) in function_bodies(trees)
        startswith(path, commons) && definition_name(reference) in public || continue
        tokens = body_tokens(body)
        length(tokens) >= CLONE_TOKENS || continue
        expected = shingles(tokens)
        for (candidate_path, candidate, candidate_shingles) in bodies
            candidate === reference && continue
            count(in(candidate_shingles), expected) >= share * length(expected) &&
                count!(found,
                    string(candidate_path, " | ",
                        something(definition_name(candidate), "anonymous"), " ~ ",
                        definition_name(reference)))
        end
    end
    return found
end

function definition_node(trees, file, line, name)
    tree = get!(() -> parse_source(file), trees, file)
    found = nothing
    visit(tree) do node
        kind(node) == K"function" && JuliaSyntax.source_line(node) == line &&
            definition_name(node) === name && (found = node)
        return found === nothing
    end
    return found
end

function argument_name(node)
    kind(node) == K"Identifier" && return node.val
    kind(node) in (K"::", K"=") && numchildren(node) == 2 && return argument_name(node[1])
    kind(node) == K"..." && numchildren(node) == 1 &&
        (name = argument_name(node[1]); return name === nothing ? nothing : (name, :...))
    return nothing
end

function forwarded_name(node)
    kind(node) == K"=" && numchildren(node) == 2 && kind(node[2]) == K"Identifier" &&
        argument_name(node[1]) === node[2].val && return node[2].val
    kind(node) == K"..." && numchildren(node) == 1 &&
        (name = forwarded_name(node[1]); return name === nothing ? nothing : (name, :...))
    return kind(node) == K"Identifier" ? node.val : nothing
end

# Positional and keyword argument names of a call, by `name`.
function call_arguments(call, name)
    positional, keywords = Any[], Any[]
    for argument in children(call)[2:end]
        if kind(argument) == K"parameters"
            append!(keywords, map(name, children(argument)))
        else
            push!(positional, name(argument))
        end
    end
    return positional, keywords
end

# A body that is one call passing the definition's own arguments, unchanged,
# in order and keywords included.
function forwards(node)
    signature = signature_call(node)
    signature === nothing && return false
    body = node[2]
    kind(body) == K"block" && numchildren(body) == 1 && (body = body[1])
    kind(body) == K"return" && numchildren(body) == 1 && (body = body[1])
    kind(body) == K"call" && is_prefix_call(body) || return false
    own = call_arguments(signature, argument_name)
    passed = call_arguments(body, forwarded_name)
    return !isempty(own[1]) && !any(isnothing, own[1]) && !any(isnothing, own[2]) &&
           own == passed
end

# The register. The functions (macros included), types and constants each package module
# defines, against its section of the register in `names`: a table of public names and one of
# private names, each name with its meaning. A name that Julia or Pkg imposes is left out.
function register_differences(tree::PackageTree, names)
    found = Dict{String, Int}()
    register = get(names, "register", Dict{String, Any}())
    mandated = names["shape"]["mandated"]
    label(path) = path == "" ? module_name(tree, tree.root) : path
    entry(path, name, problem) = (found[join((label(path), name, problem), " | ")] = 1)
    defined = Dict{String, Module}()
    for M in tree.modules
        path = module_path(tree, M)
        defined[path] = M
        section = get(register, path, Dict{String, Any}())
        public, private = get(section, "public", Dict()), get(section, "private", Dict())
        for name in keys(definitions(M))
            key = string(name)
            key in mandated && continue
            listed = haskey(public, key) ? "public" :
                     haskey(private, key) ? "private" :
                     nothing
            listed === nothing && (entry(path, key, "unregistered"); continue)
            visibility = Base.ispublic(M, name) ? "public" : "private"
            visibility == listed ||
                entry(path, key, string(visibility, ", registered ", listed))
        end
    end
    for (path, section) in register
        M = get(defined, path, nothing)
        M === nothing && (found[string(label(path), " | stale module")] = 1; continue)
        names_defined = Set(string(name) for name in keys(definitions(M)))
        public, private = get(section, "public", Dict()), get(section, "private", Dict())
        for (visibility, table) in ("public" => public, "private" => private),
            (key, meaning) in table

            key in names_defined || entry(path, key, "stale")
            visibility == "public" && haskey(private, key) && entry(path, key, "twice")
            meaning isa AbstractString && !isempty(strip(meaning)) ||
                entry(path, key, "meaning")
        end
    end
    return found
end

# The kinds whose children are a body: a definition inside one is local.
const BODIES = (K"function", K"->", K"do", K"let", K"macro")
# The statement lists in which an assignment binds a local.
const STATEMENTS = (K"block", K"local", K"global")

# Named local functions. A method defined inside a function body, an anonymous function, a
# `do` block, a `let` or a macro, and a local bound to an anonymous function, by file and
# name. Quoted code is not read.
function local_functions(files, directory)
    found = Dict{String, Int}()
    for file in files
        local_definitions!(
            found, relative(file, directory), parse_source(file), K"toplevel",
            false)
    end
    return found
end

function local_definitions!(found, path, node, parent, inside)
    k = kind(node)
    k == K"quote" && return found
    if inside && k == K"function" && signature_call(node) !== nothing
        name = definition_name(node)
        count!(found, string(path, " | ", name === nothing ? "call" : name))
    elseif inside && k == K"=" && parent in STATEMENTS && numchildren(node) == 2 &&
           anonymous(node[2])
        target = kind(node[1]) == K"::" ? node[1][1] : node[1]
        kind(target) == K"Identifier" && count!(found, string(path, " | ", target.val))
    end
    is_leaf(node) && return found
    for child in children(node)
        local_definitions!(found, path, child, k, inside || k in BODIES)
    end
    return found
end

function anonymous(node)
    kind(node) == K"->" ||
        (kind(node) == K"function" && signature_call(node) === nothing)
end

# Forwarders. The functions the package owns whose every method is one call that passes the
# method's own arguments, unchanged: a second name for that call. By module and name.
function forwarders(tree::PackageTree)
    found = Dict{String, Int}()
    trees = Dict{String, Any}()
    for M in tree.modules, (name, value) in definitions(M)

        value isa Function && !startswith(string(name), "@") || continue
        nodes = [definition_node(trees, normpath(string(m.file)), m.line, name)
                 for m in Base.methods(value) if source_path(tree, m.file) !== nothing]
        !isempty(nodes) && all(node -> node !== nothing && forwards(node), nodes) &&
            count!(found,
                string(
                    module_path(tree, M) == "" ? module_name(tree, M) :
                    module_path(tree, M),
                    " | ",
                    name))
    end
    return found
end

# Root module. The number of functions, types and constants the root module defines.
function root_definitions(tree::PackageTree)
    n = length(definitions(tree.root))
    return n == 0 ? Dict{String, Int}() : Dict(module_name(tree, tree.root) => n)
end

boxes(x::GlobalRef) = x.mod === Core && x.name === :Box
boxes(x::Expr) = any(boxes, x.args)
boxes(_) = false

# A method of `Core.kwcall`, which sorts the keywords of a call into its keyword body.
sorts(m::Method) = Base.unwrap_unionall(m.sig).parameters[1] === typeof(Core.kwcall)

# The method that each keyword body among `methods` belongs to, by the body: the method
# with keywords and the positional arguments, whose last call is its body, not the keyword
# sorter.
function keyword_bodies(methods)
    found = IdDict{Any, Method}()
    for m in methods
        (sorts(m) || isempty(Base.kwarg_decl(m))) && continue
        body = Base.bodyfunction(m)
        body === nothing || (found[body] = m)
    end
    return found
end

# The method that `m` is the keyword body of, or `m`.
keyword_method(m::Method, bodies) = get(bodies, method_subject(m), m)

# Whether a method constructs an instance of a type.
function constructs(m::Method)
    F = Base.unwrap_unionall(Base.unwrap_unionall(m.sig).parameters[1])
    return F isa DataType && F.name === Type.body.name
end

# Boxed captures. The methods whose lowered code holds a `Core.Box`, the cell of a captured
# variable that is assigned after its capture, keyword bodies and closures included. A
# keyword body counts under its function.
function boxed_captures(methods, tree::PackageTree)
    found = Dict{String, Int}()
    bodies = keyword_bodies(methods)
    for m in methods
        code = lowered_code(m)
        code === nothing && continue
        any(boxes, code) &&
            count!(found,
                string(module_name(tree, m.module), " | ",
                    keyword_method(m, bodies).name, " | ", source_name(tree, m)))
    end
    return found
end

# A `Val` type, or a type variable bounded by one.
tagged(p::DataType) = p.name === Val.body.name
tagged(p::TypeVar) = tagged(Base.unwrap_unionall(p.ub))
tagged(_) = false

# The number of `Val(x)` calls in lowered code whose argument is not a literal.
function computed_tags(code)
    resolved = Dict{Int, Any}()
    resolve(x::GlobalRef) = isdefined(x.mod, x.name) ? getglobal(x.mod, x.name) : nothing
    resolve(x::Core.SSAValue) = get(resolved, x.id, nothing)
    resolve(_) = nothing
    n = 0
    for (i, statement) in enumerate(code)
        statement isa GlobalRef && (resolved[i] = resolve(statement))
        call = statement isa Expr && statement.head === :(=) ? statement.args[2] : statement
        call isa Expr && call.head === :call && length(call.args) == 2 &&
            resolve(call.args[1]) === Val &&
            call.args[2] isa
            Union{Core.SlotNumber, Core.SSAValue, Core.Argument, GlobalRef} &&
            (n += 1)
    end
    return n
end

# `Val` tags. A method whose first positional argument is a `Val` is a constructor or a
# method of a function that the register lists under the `lifts` of the module that owns the
# function, wherever the method is defined; so is a method that calls `Val(x)` with an
# argument that is not a literal, a keyword body counting as its method; the keyword sorters
# only pass the call on. A listed lift that is no function of its module is stale.
function val_tags(methods, tree::PackageTree, names)
    found = Dict{String, Int}()
    register = get(names, "register", Dict{String, Any}())
    label(M) = (
        path = module_path(tree, M); path == "" ? module_name(tree, tree.root) : path)
    lifts(M) = get(get(register, module_path(tree, M), Dict{String, Any}()), "lifts",
        String[])
    bodies = keyword_bodies(methods)
    for m in methods
        sorts(m) && continue
        owner = keyword_method(m, bodies)
        constructs(owner) && continue
        home = value_owner(method_subject(owner))
        string(owner.name) in lifts(home in tree.modules ? home : owner.module) && continue
        key(what) = string(label(owner.module), " | ", owner.name, " | ", what)
        signature = Base.unwrap_unionall(m.sig).parameters
        owner === m && length(signature) >= 2 && tagged(signature[2]) &&
            count!(found, key("Val first"))
        code = lowered_code(m)
        n = code === nothing ? 0 : computed_tags(code)
        n > 0 && (found[key("Val(x)")] = get(found, key("Val(x)"), 0) + n)
    end
    for M in tree.modules, name in lifts(M)

        get(definitions(M), Symbol(name), nothing) isa Function ||
            count!(found, string(label(M), " | ", name, " | stale lift"))
    end
    return found
end

# The identifier `ID` of a formula type `Formula{ID}` of the type name `name`, or `nothing`.
function formula_identifier(T, name)
    T = Base.unwrap_unionall(T)
    return T isa DataType && T.name === name && !isempty(T.parameters) &&
           T.parameters[1] isa Symbol ? T.parameters[1] : nothing
end

# The formula registry. Per module that defines a type `Formula`, a constant `FORMULAS` or a
# `formulas` table in its section of the register: the identifiers with a method of Commons'
# `defaults` on their formula type, `FORMULAS` and the table hold the same identifiers, and
# each identifier has a method of a package function outside Commons whose first positional
# argument is its formula. An empty `FORMULAS` and no table agree.
function formula_registry(methods, tree::PackageTree, names)
    found = Dict{String, Int}()
    register = get(names, "register", Dict{String, Any}())
    C = commons_module(tree)
    defaults = C !== nothing && isdefined(C, :defaults) ? getfield(C, :defaults) : nothing
    for M in tree.modules
        defined = definitions(M)
        path = module_path(tree, M)
        table = get(get(register, path, Dict{String, Any}()), "formulas", Dict{
            String, Any}())
        F = get(defined, :Formula, nothing)
        haskey(defined, :FORMULAS) || F isa Type || !isempty(table) || continue
        name = F isa Type ? Base.unwrap_unionall(F).name : nothing
        declared, physics = Set{Symbol}(), Set{Symbol}()
        for m in methods
            sorts(m) && continue
            signature = Base.unwrap_unionall(m.sig).parameters
            length(signature) >= 2 || continue
            subject = method_subject(m)
            if subject === defaults
                T = Base.unwrap_unionall(signature[2])
                T isa DataType && T.name === Type.body.name || continue
                id = formula_identifier(T.parameters[1], name)
                id === nothing || push!(declared, id)
            elseif subject isa Function && value_owner(subject) in tree.modules &&
                   !(C !== nothing && within(value_owner(subject), C, tree.root))
                id = formula_identifier(signature[2], name)
                id === nothing || push!(physics, id)
            end
        end
        listed = Set{Symbol}(get(defined, :FORMULAS, ()))
        tabled = Set{Symbol}(Symbol(key) for key in keys(table))
        label = path == "" ? module_name(tree, tree.root) : path
        for id in union(declared, listed, tabled),
            (what, holders) in (
                ("defaults", declared), ("FORMULAS", listed), ("table", tabled),
                ("physics", physics))

            id in holders || count!(found, string(label, " | ", id, " | ", what))
        end
    end
    return found
end

# A field type that admits values of more than one concrete layout: an abstract type, a
# `Union`, a `UnionAll`, a type parameter bounded by `Any`, or a type with such a type among
# its parameters (`Vector{Any}`, `Dict{Symbol, Any}`, `Vector{AbstractCablePart}`).
loose(F::TypeVar) = F.ub === Any
loose(F::DataType) = isabstracttype(F) || any(loose_parameter, F.parameters)
loose(_) = true

# A parameter of a field type: a type variable is fixed per instance, and a value such as an
# integer or a symbol is no type.
loose_parameter(p::Type) = loose(p)
loose_parameter(p::Core.TypeofVararg) = isdefined(p, :T) && loose_parameter(p.T)
loose_parameter(_) = false

# Typed values. Per package struct, each field of a loose type, and per constructor method of
# the struct, the positional values typed `Any`; Julia's converting default constructor counts.
function typed_values(tree::PackageTree)
    found = Dict{String, Int}()
    for M in tree.modules, (name, value) in definitions(M)

        value isa Type || continue
        T = Base.unwrap_unionall(value)
        T isa DataType && isstructtype(T) || continue
        key(what) = string(module_name(tree, M), " | ", name, " | ", what)
        for (field, F) in zip(fieldnames(T), T.types)
            loose(F) && count!(found, key(field))
        end
        for m in Base.methods(value)
            m.module in tree.modules || continue
            F = Base.unwrap_unionall(m.sig).parameters[1]
            F isa DataType && F.name === Type.body.name && F.parameters[1] === value ||
                continue
            untyped = count(p -> upper(p) === Any, Base.unwrap_unionall(m.sig).parameters[2:end])
            untyped == 0 && continue
            what = key(string("constructor | ", source_name(tree, m)))
            found[what] = get(found, what, 0) + untyped
        end
    end
    return found
end

# The path of a module below the root, or from its top extension module for an extension.
function module_path(tree::PackageTree, M::Module)
    M === tree.root && return ""
    M in tree.extensions || return module_name(tree, M)
    parts = [nameof(M)]
    # A package extension is a top-level module: its parent is itself.
    while parentmodule(M) !== M && parentmodule(M) in tree.extensions
        M = parentmodule(M)
        pushfirst!(parts, nameof(M))
    end
    return join(parts, ".")
end

# The types that a signature parameter admits: the members of a `Union`, the bound of a type
# variable.
function admitted_types(T)
    T isa Union ? [admitted_types(T.a); admitted_types(T.b)] :
    T isa TypeVar ? admitted_types(T.ub) :
    T isa Core.TypeofVararg ? Any[] : Any[Base.unwrap_unionall(T)]
end

# Names. A name's chunks are its lowercase pieces between underscores and case transitions,
# without leading underscores, a leading `@` or a trailing `!`. Digits stay with their piece.
function chunks(name)
    text = replace(string(name), r"^[_@]+" => "", r"!$" => "")
    found = String[]
    for part in split(text, '_'; keepempty = false)
        for piece in
            eachmatch(r"[A-Z]+(?=[A-Z][a-z])|[A-Z]?[a-z0-9]+|[A-Z]+[0-9]*|[^\W_]+", part)
            push!(found, lowercase(piece.match))
        end
    end
    return found
end

# The role of the terms in each list of a section of the names file. A preposition is
# sanctioned only as a whole field or keyword name.
const ROLES = Dict("verbs" => :verb, "nouns" => :noun, "qualifiers" => :qualifier,
    "symbols" => :symbol, "prepositions" => :preposition)

# The sanctioned terms that a name in `path` may use, with their roles: the global section and
# the sections of the module path's prefixes, root side first.
function vocabulary(names, path)
    found = Dict{String, Symbol}()
    sections = Any[names["global"]]
    parts = split(path, '.'; keepempty = false)
    for i in eachindex(parts)
        section = get(get(names, "modules", Dict()), join(parts[1:i], "."), nothing)
        section === nothing || push!(sections, section)
    end
    for section in sections, (role, entries) in section

        for term in (entries isa AbstractDict ? keys(entries) : entries)
            found[term] = ROLES[role]
        end
    end
    for term in names["shape"]["predicates"]
        found[term] = :predicate
    end
    found[names["shape"]["marker"]] = :marker
    return found
end

# The splits of `name` into terms with the fewest terms. A term matches whole chunks, several
# chunks joined by underscores (a compound) or part of a chunk (a squashed word).
function splits(name, terms)
    text = collect(join(chunks(name), "_"))
    candidates = [collect(term) for term in keys(terms)]
    memo = Dict{Int, Vector{Vector{String}}}()
    function from(i)
        i > length(text) && return [String[]]
        text[i] == '_' && return from(i + 1)
        haskey(memo, i) && return memo[i]
        best = Vector{String}[]
        for term in candidates
            j = i + length(term) - 1
            j <= length(text) && text[i:j] == term || continue
            for rest in from(j + 1)
                option = [String(term); rest]
                if isempty(best) || length(option) < length(first(best))
                    best = [option]
                elseif length(option) == length(first(best))
                    push!(best, option)
                end
            end
        end
        return memo[i] = best
    end
    return from(1)
end

# Every name the package owns: modules below the root, types, functions (macros included),
# constants, struct fields and method keywords, each with its kind, its module path, whether
# the root exports it and its value. A keyword belongs to the module that defines its
# method, which may extend a function of another module. The names Julia or Pkg impose are
# left out; an extension module is no module's submodule, so its own name is never read.
function owned_names(tree::PackageTree, mandated)
    found = NamedTuple{
        (:kind, :path, :name, :exported, :value), Tuple{Symbol, String, Symbol, Bool, Any}}[]
    path(M) = module_path(tree, M)
    exported(M,
        name,
        value) = Base.isexported(tree.root, name) &&
                 isdefined(tree.root, name) && getfield(tree.root, name) === value
    for M in tree.modules
        for name in names(M; all = true)
            (startswith(string(name), "#") || name in (:eval, :include) ||
             string(name) in mandated || name === nameof(M)) && continue
            isdefined(M, name) || continue
            value = getfield(M, name)
            value isa Module && parentmodule(value) === M && value !== M &&
                push!(found,
                    (kind = :module, path = path(M), name,
                        exported = exported(M, name, value), value))
        end
        for (name, value) in definitions(M)
            string(name) in mandated && continue
            kind = value isa Function ? :function : value isa Type ? :type : :constant
            push!(found, (;
                kind, path = path(M), name, exported = exported(M, name, value), value))
            if value isa Type
                T = Base.unwrap_unionall(value)
                T isa DataType && isstructtype(T) &&
                    for field in fieldnames(T)
                        field isa Symbol && push!(found,
                            (kind = :field, path = path(M), name = field,
                                exported = false, value = nothing))
                    end
            end
            if value isa Function || value isa Type
                for m in Base.methods(value), keyword in Base.kwarg_decl(m)

                    m.module in tree.modules && !endswith(string(keyword), "...") &&
                        push!(found,
                            (kind = :keyword,
                                path = path(parentmodule(m)), name = keyword,
                                exported = false, value = nothing))
                end
            end
        end
    end
    return unique(item -> (item.kind, item.path, item.name), found)
end

# The form of a name, from the roles of its terms: `nothing` when it conforms.
function form(kind, roles, abstract)
    roles = copy(roles)
    kind === :type && abstract && !isempty(roles) && roles[1] === :marker &&
        popfirst!(roles)
    (isempty(roles) || :marker in roles) && return "form"
    :preposition in roles &&
        return kind in (:field, :keyword) && roles == [:preposition] ? nothing : "form"
    if kind === :function
        roles[1] === :predicate &&
            return length(roles) == 2 && !(roles[2] in (:verb, :predicate)) ? nothing :
                   "form"
        roles[1] === :verb &&
            return length(roles) == 1 || roles[2] in (:noun, :symbol) ? nothing : "form"
        (:verb in roles || :predicate in roles) && return "form"
        return length(roles) == 1 || roles[end] in (:noun, :symbol) ? nothing : "form"
    end
    (:verb in roles || :predicate in roles) && return "form"
    kind in (:field, :keyword) && roles == [:qualifier] && return nothing
    return roles[end] in (:noun, :symbol) ? nothing : "form"
end

# The chunks of a module path, without the root's chunks and an extension's `ext`, as the
# context of the names it holds.
function context(tree::PackageTree, path)
    root = Set(chunks(nameof(tree.root)))
    return Set(c for part in split(path, '.'; keepempty = false)
    for c in chunks(part)
    if !(c in root) && c != "ext")
end

# Names. The checks of the naming rule over every name the package owns, with the sanctioned
# terms of `names`: terms, shape, subject, context, verb homes and Base names. Each violation
# counts once, keyed by check, kind, module path and name.
function naming_rule(methods, tree::PackageTree, names)
    found = Dict{String, Int}()
    shape = names["shape"]
    rejected = names["rejected"]
    verbs = keys(names["global"]["verbs"])
    C = commons_module(tree)
    base = Set(string(n) for M in (Base, Core) for n in Base.names(M))
    split_of = Dict{Tuple{String, Symbol}, Vector{String}}()
    entry(check,
        item,
        detail = "") = (found[join(
        filter(!isempty, [check, string(item.kind),
            item.path, string(item.name), detail]), " | ")] = 1)
    items = owned_names(tree, shape["mandated"])
    for item in items
        terms = vocabulary(names, item.path)
        options = splits(item.name, terms)
        # The rejected words are looked up among the terms of the name's one split, or among its
        # chunks when it has none, so that a sanctioned compound counts as one term.
        for word in (length(options) == 1 ? only(options) : chunks(item.name))
            haskey(rejected, word) &&
                entry("terms", item, string("rejected ", word, " -> ", rejected[word]))
        end
        if length(options) != 1
            entry("terms", item, isempty(options) ? "unknown" : "ambiguous")
            continue
        end
        split_of[(item.path, item.name)] = only(options)
        roles = [terms[t] for t in only(options)]
        abstract = item.value isa Type && isabstracttype(Base.unwrap_unionall(item.value))
        cap = shape[string(item.kind)] + (abstract && roles[1] === :marker ? 1 : 0)
        length(roles) > cap && entry("shape", item, "length")
        length(roles) <= cap && form(item.kind, roles, abstract) !== nothing &&
            entry("shape", item, "form")
        if item.kind in (:function, :type, :constant, :module) && !item.exported
            clash = intersect(Set(c for t in only(options) for c in split(t, '_')), context(tree, item.path))
            isempty(clash) || entry("context", item, join(sort!(collect(clash)), ","))
        end
        if item.kind === :function
            roles[1] === :verb && only(options)[1] in verbs &&
                (C === nothing || !within(parentmodule(item.value), C, tree.root)) &&
                entry("verbs", item, only(options)[1])
            string(item.name) in base && entry("base", item)
        end
    end
    # Compare the terms of each function with the terms of the package types its methods take.
    for item in items
        item.kind === :function && haskey(split_of, (item.path, item.name)) || continue
        own = Set(split_of[(item.path, item.name)])
        for m in Base.methods(item.value),
            p in Base.unwrap_unionall(m.sig).parameters[2:end]

            for T in admitted_types(p)
                T isa DataType && T.name.module in tree.modules || continue
                type_terms = get(split_of, (module_path(tree, T.name.module), T.name.name), String[])
                clash = intersect(own, type_terms)
                isempty(clash) || entry("subject", item, join(sort!(collect(clash)), ","))
            end
        end
    end
    return merge!(found, naming_file(names, split_of))
end

# The file of sanctioned terms. Each term is used by a name in its section's scope and is
# sanctioned once along a module path. No term is a preposition or a rejected word. A global
# noun, qualifier or symbol is used by root names or by names of two top-level modules.
function naming_file(names, split_of)
    found = Dict{String, Int}()
    entry(section, term, problem) = (found[join(("file", section, term, problem), " | ")] = 1)
    users = Dict{String, Set{String}}()
    for ((path, _), terms) in split_of, term in terms

        push!(get!(Set{String}, users, term), path)
    end
    sections = Pair{String, Any}["global" => names["global"]]
    for (path, section) in get(names, "modules", Dict())
        push!(sections, path => section)
    end
    sanctioned(section) = Set(term for (_, entries) in section
    for term in (entries isa AbstractDict ? keys(entries) : entries))
    for (path, section) in sections
        seen = Set{String}()
        above = path == "global" ? Set{String}() :
                union(sanctioned(names["global"]),
            (sanctioned(names["modules"][prefix])
            for prefix in prefixes(path)
            if prefix != path && haskey(names["modules"], prefix))...)
        for (role, entries) in section,
            term in (entries isa AbstractDict ? keys(entries) : entries)

            (term in seen || term in above) && entry(path, term, "sanctioned twice")
            push!(seen, term)
            role != "prepositions" && term in names["shape"]["prepositions"] &&
                entry(path, term, "preposition")
            haskey(names["rejected"], term) && entry(path, term, "rejected")
            paths = get(users, term, Set{String}())
            scope = path == "global" ? paths :
                    Set(p for p in paths if p == path || startswith(p, path * "."))
            isempty(scope) && entry(path, term, "unused")
            path == "global" && role != "verbs" && !isempty(paths) && !("" in paths) &&
                length(Set(first(split(p, '.')) for p in paths)) < 2 &&
                entry(path, term, "one module")
        end
    end
    return found
end

# The module path and its prefixes, root side first.
function prefixes(path)
    (parts = split(path, '.'); [join(parts[1:i], ".") for i in eachindex(parts)])
end

end
