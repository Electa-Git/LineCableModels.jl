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

const WORKSPACE_FIELDS = (:input, :plan, :buffers)
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

# Storage vocabulary. A `*Workspace` struct has exactly the fields `input`, `plan` and
# `buffers`, and optionally `trace`. Its entry counts the missing and the other fields.
# No function name ends in `_workspace` or `_storage`, and only `initialize_buffers`
# ends in `_buffers`. Only Commons defines a `*Buffers` struct, and no type name ends
# in `Buffer`, `Storage`, `Scratch` or `Cache`. No field, named tuple key, parameter,
# local or loop variable takes the name `work`, `scratch`, `storage` or `cache`.
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
                if endswith(text, "Workspace")
                    entry(name,
                        count(!in(fields), WORKSPACE_FIELDS) +
                        count(!in((WORKSPACE_FIELDS..., :trace)), fields))
                end
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

end
