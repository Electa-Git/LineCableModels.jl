# Each guard takes its inputs as arguments. The controls below apply it to probe
# packages and probe sources. `test/tools/architecture_inventory.jl` evaluates
# this module to print the live inventory.
@testmodule ArchitectureGuards begin
    using LineCableModels
    using JuliaSyntax: JuliaSyntax, SyntaxNode, @K_str, kind, children, numchildren,
        is_leaf, is_prefix_call
    using TOML: TOML
    import DataFrames, LinearAlgebra, Statistics, Random, Logging
    # These packages load the extensions checked by `explicit_imports.jl`.
    import Measurements, Distributions, Gmsh, Calculus, XLSX, CairoMakie

    # A module may reference only earlier modules, its ancestors and its
    # descendants. The root module and package extensions are exempt. The load
    # order is defined once, in the test taxonomy.
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "support", "taxonomy.jl"))
    const ORDER = keys(MODULE_OWNERS)
    const EXTENSIONS = (:LineCableModelsMeasurementsExt, :LineCableModelsDistributionsExt,
        :LineCableModelsGmshExt, :LineCableModelsXLSXExt, :LineCableModelsMakieExt,
        :LineCableModelsCairoMakieExt)
    const DEPENDENCIES = (Base, DataFrames, LinearAlgebra, Statistics, Random,
        Base.require(LineCableModels, :Dates), Logging)
    # Per-family formula selections are the only name shared across modules.
    const SHARED_NAMES = (:Formula,)
    const RESERVED_VERB = r"^_*(validate|check|require|assert|verify|ensure)_"
    # Tables of the structural guards (A) and of the Commons and helper guards (C).
    const A_TABLES = ("ownership", "placement", "direction", "names", "shadowing",
        "validate", "reserved_verbs", "switches")
    const C_TABLES = ("commons", "vocabulary", "fingerprints", "clones", "helpers", "root")
    const TABLES = (A_TABLES..., C_TABLES...)
    const SWITCHES = ("applicable", "eval", "kind")
    const BASELINE = joinpath(@__DIR__, "architecture_baseline.toml")

    const COMMONS = :Commons
    const COMMONS_DIRECTORY = "src/commons"
    const CONSTANTS_FILE = "src/commons/consts.jl"
    const COMMONS_TESTS = "test/unit/commons"
    # Directories under `test/` whose files name definitions without testing them.
    const NOT_TESTS = ("quality", "tools")
    # Each Commons public name and the definition names it reserves outside Commons.
    const VOCABULARY = Dict{Symbol, Regex}(
        :vacuum_permittivity =>
            r"(?i)^_*(\w+_)?(vacuum_permittivity|eps(ilon)?_?[0₀]|[εϵ]_?[0₀])$",
        :vacuum_permeability =>
            r"(?i)^_*(\w+_)?(vacuum_permeability|mu_?[0₀]|[μµ]_?[0₀])$",
        :ideal_transposition! => r"transpos",
        :kron_reduce! => r"kron",
        :merge_bundles! => r"merge_bundle",
        :bundle_operations => r"^_*bundle_(operations|pairs)$",
        :reorder_indices => r"^_*reorder(_|$)",
        :ReductionPlan => r"(?i)^_*reduction_?(plan|map)$",
        :ReductionBuffers => r"(?i)^_*reduction_?buffers$",
        :reduce_line_matrices! => r"^_*reduce_(line|primitive)_matrices")
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
    within(m, owner::Module, root::Module) =
        m isa Module && (m === owner || (owner !== root && is_descendant(m, owner)))

    function module_name(tree::PackageTree, m::Module)
        m === tree.root && return string(nameof(m))
        path, prefix = fullname(m), fullname(tree.root)
        length(path) > length(prefix) && path[1:length(prefix)] == prefix &&
            return join(path[length(prefix)+1:end], ".")
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
    source_name(tree::PackageTree, m::Method) =
        something(source_path(tree, m.file), string(m.file))

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

    function mentioned!(found::Base.IdSet{DataType}, @nospecialize(T))
        if T isa TypeVar
            mentioned!(found, T.lb)
            mentioned!(found, T.ub)
        elseif T isa UnionAll
            mentioned!(found, T.var)
            mentioned!(found, T.body)
        elseif T isa Union
            mentioned!(found, T.a)
            mentioned!(found, T.b)
        elseif T isa Core.TypeofVararg
            isdefined(T, :T) && mentioned!(found, T.T)
        elseif T isa DataType && T ∉ found
            push!(found, T)
            for p in T.parameters
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
            count!(found, string(module_name(tree, M), " | ", module_name(tree, F), ".",
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

    lowered_code(m::Method) = isdefined(m, :source) && m.source !== nothing &&
        !isdefined(m, :generator) ? Base.uncompressed_ir(m).code : nothing

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
                count!(found, string(module_name(tree, M), " -> ", module_name(tree, owner),
                    " | ", source_name(tree, m)))
            end
        end
        return found
    end

    # Names. Distinct functions or types, owned by different package modules,
    # that share a name. The value is the number of modules that define the name.
    function shared_names(tree::PackageTree; exceptions = SHARED_NAMES)
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
            for (name, modules) in owners if length(modules) > 1 && name ∉ exceptions)
    end

    # Shadowing. A public package-owned name that a dependency exports for a
    # different object.
    function shadowing(tree::PackageTree, dependencies)
        found = Dict{String, Int}()
        for M in tree.modules, name in names(M)
            isdefined(M, name) || continue
            value = getfield(M, name)
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

    is_required_block(node) = kind(node) == K"macrocall" && numchildren(node) > 0 &&
        syntax_name(node[1]) in ((Symbol("@required"),),
            (:RequiredInterfaces, Symbol("@required")))

    # validate returns its subject. Every `validate` method names its first positional argument,
    # returns it on every path, or throws, and uses it in its work unless it only returns it.
    function validate_returns(files, directory)
        found = Dict{String, Int}()
        for file in files
            visit(parse_source(file)) do node
                is_required_block(node) && return false
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

    symbol_literal(node) = kind(node) == K"quote" && numchildren(node) == 1 &&
        kind(node[1]) == K"Identifier"
    symbol_operand(node) = symbol_literal(node) ||
        (kind(node) in (K"tuple", K"vect") && numchildren(node) > 0 &&
         all(symbol_literal, children(node)))

    function kind_field(node)
        kind(node) == K"." && numchildren(node) == 2 || return false
        member = node[2]
        kind(member) == K"quote" && numchildren(member) == 1 && (member = member[1])
        return kind(member) == K"Identifier" && member.val === :kind
    end

    # Negated comparisons count too, so negating a switch does not hide it.
    const KIND_COMPARISONS = ((:(===),), (:(==),), (:in,), (:∈,),
        (:(!==),), (:(!=),), (:∉,))

    function kind_switch(node)
        kind(node) == K"call" && numchildren(node) == 3 || return false
        operator, left, right = is_prefix_call(node) ?
            (node[1], node[2], node[3]) : (node[2], node[1], node[3])
        return syntax_name(operator) in KIND_COMPARISONS &&
            kind_field(left) && symbol_operand(right)
    end

    # Symbol switches and probes. Per file, calls to `applicable`, `@eval` calls and comparisons of a
    # `kind` field against symbols, including negated comparisons.
    function switches(files, directory)
        found = Dict{String, Dict{String, Int}}()
        for file in files
            counts = Dict(name => 0 for name in SWITCHES)
            visit(parse_source(file)) do node
                k = kind(node)
                if k == K"call" && is_prefix_call(node) && syntax_name(node[1]) in
                        ((:applicable,), (:Base, :applicable), (:Core, :applicable))
                    counts["applicable"] += 1
                elseif k == K"macrocall" && numchildren(node) > 0 && syntax_name(node[1]) in
                        ((Symbol("@eval"),), (:Base, Symbol("@eval")), (:Core, Symbol("@eval")))
                    counts["eval"] += 1
                end
                kind_switch(node) && (counts["kind"] += 1)
                return true
            end
            any(>(0), values(counts)) && (found[relative(file, directory)] = counts)
        end
        return found
    end

    type_name(T) = (T = Base.unwrap_unionall(T); T isa DataType ? nameof(T) : nothing)

    # The functions, types and constants that a module defines, by name.
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
            owner === M && (found[name] = value)
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

    # The texts of the checkout's `.jl` files under `directory`, except those under
    # `excluded` subdirectories.
    function test_sources(directory; excluded = ())
        texts = String[]
        for file in repository_files(directory)
            endswith(file, ".jl") &&
                !any(x -> inside(file, joinpath(directory, x)), excluded) &&
                push!(texts, read(file, String))
        end
        return texts
    end

    named(texts, name) = (pattern = Regex("(?<![\\w!@])\\Q" * string(name) * "\\E(?![\\w!])");
        any(text -> occursin(pattern, text), texts))

    commons_module(tree::PackageTree) = (C = isdefined(tree.root, COMMONS) ?
        getfield(tree.root, COMMONS) : nothing; C isa Module ? C : nothing)

    # Hooks are Commons functions with a method defined outside Commons. The
    # user API is the Commons names that the root module exports or declares public.
    function commons_roles(methods, tree::PackageTree, C::Module, defined)
        functions = IdDict{Any, Symbol}(value => name for (name, value) in defined
            if value isa Function)
        hooks = Set{Symbol}()
        for m in methods
            within(m.module, C, tree.root) && continue
            name = get(functions, method_subject(m), nothing)
            name === nothing || push!(hooks, name)
        end
        api = Set(name for (name, value) in defined if Base.ispublic(tree.root, name) &&
            isdefined(tree.root, name) && getfield(tree.root, name) === value)
        return hooks, api
    end

    # Commons admission. Each Commons name is public, documented and named in a Commons unit
    # test. Unless it is user API, methods of at least two owners outside
    # Commons use it, directly or through public Commons definitions. Unless it
    # is a hook or user API, it has a vocabulary entry. Macros are used where
    # source code calls them.
    function commons_admission(methods, tree::PackageTree, homes, vocabulary)
        found = Dict{String, Int}()
        C = commons_module(tree)
        C === nothing && return found
        defined = definitions(C)
        hooks, api = commons_roles(methods, tree, C, defined)
        name_of = IdDict{Any, Symbol}(value => name for (name, value) in defined)
        # Owners of every Commons value, nested functions and keyword bodies included.
        owners = IdDict{Any, Set{Module}}(value => Set{Module}() for value in values(defined))
        # The Commons values that each Commons method subject references.
        references = IdDict{Any, Base.IdSet{Any}}()
        unnamed(value) = !haskey(name_of, value) && (value isa Function &&
            within(parentmodule(typeof(value)), C, tree.root) ||
            value isa Type && value <: Function && within(parentmodule(value), C, tree.root))
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
        public = Set(name for name in keys(defined) if Base.ispublic(C, name))
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
        users = Dict(name => owners[value] for (name, value) in defined)
        tests = test_sources(joinpath(tree.directory, COMMONS_TESTS))
        for name in keys(defined)
            for (criterion, met) in (("public", name in public),
                    ("docstring", Base.Docs.hasdoc(C, name)),
                    ("owners", name in api || length(users[name]) >= 2),
                    ("tests", named(tests, name)),
                    ("vocabulary", name in hooks || name in api || haskey(vocabulary, name)))
                met || (found[string(name, " | ", criterion)] = 1)
            end
        end
        return found
    end

    function reserving(vocabulary, name)
        for (owner, pattern) in vocabulary
            occursin(pattern, string(name)) && return owner
        end
        return nothing
    end

    function calls(node, owner)
        kind(node) == K"call" && is_prefix_call(node) || return false
        name = syntax_name(node[1])
        return name == (owner,) || (length(name) >= 2 && name[end-1:end] == (COMMONS, owner))
    end

    # `=` under these kinds binds a keyword, a default or a field.
    const NOT_ASSIGNED = (K"call", K"dotcall", K"parameters", K"tuple", K"vect", K"braces",
        K"curly", K"ref")

    function assigned_names(target)
        kind(target) == K"Identifier" && return [target.val]
        kind(target) == K"::" && numchildren(target) == 2 && return assigned_names(target[1])
        kind(target) == K"tuple" && return reduce(vcat, map(assigned_names, children(target));
            init = Symbol[])
        return Symbol[]
    end

    # Reserved vocabulary. Outside Commons, no function, constant or assigned local takes a
    # name reserved by a Commons definition. A local assigned from a call to
    # the reserving definition is exempt.
    function reserved_vocabulary(files, directory, vocabulary)
        found = Dict{String, Int}()
        commons = joinpath(directory, COMMONS_DIRECTORY)
        for file in files
            inside(file, commons) && continue
            path = relative(file, directory)
            visit(parse_source(file)) do node
                if kind(node) == K"function"
                    name = definition_name(node)
                    name isa Symbol && reserving(vocabulary, name) !== nothing &&
                        count!(found, string(path, " | ", name))
                elseif kind(node) == K"=" && numchildren(node) == 2 &&
                        !(node.parent !== nothing && kind(node.parent) in NOT_ASSIGNED)
                    for name in assigned_names(node[1])
                        owner = reserving(vocabulary, name)
                        owner === nothing || calls(node[2], owner) ||
                            count!(found, string(path, " | ", name))
                    end
                end
                return true
            end
        end
        return found
    end

    any_node(predicate, node) = predicate(node) ||
        (!is_leaf(node) && any(child -> any_node(predicate, child), children(node)))

    # The innermost statement around a node.
    function statement(node)
        while node.parent !== nothing && !(kind(node.parent) in (K"block", K"toplevel"))
            node = node.parent
        end
        return node
    end

    is_pi(node) = kind(node) == K"Identifier" && node.val in (:π, :pi)
    is_number(node) = kind(node) in (K"Integer", K"Float", K"Float32", K"HexInt", K"OctInt",
        K"BinInt")

    # `base ^ -7` with a literal 10 in `base`.
    function minus_seventh_power(node)
        kind(node) == K"call" && numchildren(node) == 3 || return false
        operator, base, exponent = is_prefix_call(node) ?
            (node[1], node[2], node[3]) : (node[2], node[1], node[3])
        return syntax_name(operator) == (:^,) && kind(exponent) == K"Integer" &&
            exponent.val == -7 && any_node(n -> kind(n) == K"Integer" && n.val == 10, base)
    end

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

    # Literal fingerprints. Numeric fingerprints of the Commons constants appear only in the
    # constants file: digits of ε₀ or c₀, and 10⁻⁷ in a statement with π.
    function fingerprints(files, directory)
        found = Dict{String, Int}()
        for file in files
            path = relative(file, directory)
            path == CONSTANTS_FILE && continue
            visit(parse_source(file)) do node
                fingerprint(node) && count!(found, path)
                return true
            end
        end
        return found
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
                (node.parent !== nothing && kind(node.parent) in (K"block", K"toplevel"))
            statement && throw_guard(node) || return true
            bytes[JuliaSyntax.byte_range(node) .- offset] .= UInt8(' ')
            return false
        end
        text = String(bytes)
        raw = [(kind(t), JuliaSyntax.untokenize(t, text)) for t in JuliaSyntax.tokenize(text)
            if !(kind(t) in (K"Whitespace", K"NewlineWs", K"Comment"))]
        tokens = String[]
        for (i, (k, token)) in enumerate(raw)
            if k == K"Identifier" && !Base.isoperator(Symbol(token))
                previous = i > 1 ? raw[i-1][2] : ""
                next = i < length(raw) ? raw[i+1][2] : ""
                token = previous == "." ? "FIELD" : next == "(" ? "CALL" : "NAME"
            end
            push!(tokens, token)
        end
        return tokens
    end

    shingles(tokens) = Set(join(view(tokens, i:i+SHINGLE-1), '\x1f')
        for i in 1:length(tokens)-SHINGLE+1)

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
        public = Set(name for name in names(C; all = true)
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
                    count!(found, string(candidate_path, " | ",
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

    function body_statements(node)
        body = node[2]
        return kind(body) == K"block" ? numchildren(body) : 1
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

    # Tiny private helpers. Private module-level functions with one method and at most three body
    # statements, referenced by exactly one method and named in no test file, are
    # tiny helpers. Exact forwarders count whatever their references.
    function tiny_helpers(methods, tree::PackageTree, tests)
        callers = IdDict{Any, Set{Method}}()
        for m in methods
            code = lowered_code(m)
            code === nothing && continue
            for value in referenced_values(code)
                value isa Function && push!(get!(Set{Method}, callers, value), m)
            end
        end
        trees = Dict{String, Any}()
        found = Dict{String, Int}()
        for M in tree.modules, (name, value) in definitions(M)
            value isa Function && !startswith(string(name), "@") &&
                !Base.ispublic(M, name) && length(Base.methods(value)) == 1 || continue
            m = only(Base.methods(value))
            path = source_path(tree, m.file)
            path === nothing && continue
            node = definition_node(trees, normpath(string(m.file)), m.line, name)
            (node !== nothing && body_statements(node) <= 3) || continue
            named(tests, name) && continue
            references = length(setdiff(get(callers, value, Set{Method}()), (m,)))
            (references == 1 || forwards(node)) && count!(found, path)
        end
        return found
    end

    # Root freeze. The number of functions, types and constants the root module defines.
    function root_definitions(tree::PackageTree)
        n = length(definitions(tree.root))
        return n == 0 ? Dict{String, Int}() : Dict(module_name(tree, tree.root) => n)
    end

    function inventory(tree::PackageTree; order = ORDER, dependencies = DEPENDENCIES,
            vocabulary = VOCABULARY, tables = TABLES)
        methods = package_methods(tree)
        files = source_files(tree.directory)
        homes = home_directories(tree)
        guards = Dict{String, Function}(
            "ownership" => () -> ownership(methods, tree),
            "placement" => () -> placement(methods, tree, homes),
            "direction" => () -> direction(methods, tree, order),
            "names" => () -> shared_names(tree),
            "shadowing" => () -> shadowing(tree, dependencies),
            "validate" => () -> validate_returns(files, tree.directory),
            "reserved_verbs" => () -> reserved_verbs(files, tree.directory),
            "switches" => () -> switches(files, tree.directory),
            "commons" => () -> commons_admission(methods, tree, homes, vocabulary),
            "vocabulary" => () -> reserved_vocabulary(files, tree.directory, vocabulary),
            "fingerprints" => () -> fingerprints(files, tree.directory),
            "clones" => () -> clones(files, methods, tree),
            "helpers" => () -> tiny_helpers(methods, tree,
                test_sources(joinpath(tree.directory, "test"); excluded = NOT_TESTS)),
            "root" => () -> root_definitions(tree))
        return Dict{String, Any}(table => guards[table]() for table in tables)
    end

    function live_tree()
        extensions = map(EXTENSIONS) do name
            extension = Base.get_extension(LineCableModels, name)
            extension === nothing && error("Extension $name is not loaded")
            extension
        end
        return PackageTree(LineCableModels, extensions, pkgdir(LineCableModels))
    end

    # Writes a probe package into a temporary directory, loads it and returns
    # its inventory.
    function probe_inventory(files, name; extension = nothing, options...)
        mktempdir() do directory
            for (path, text) in files
                mkpath(dirname(joinpath(directory, path)))
                write(joinpath(directory, path), text)
            end
            root = Base.include(Main, joinpath(directory, "src", name * ".jl"))
            extensions = extension === nothing ? Module[] :
                [Base.include(Main, joinpath(directory, "ext", extension * ".jl"))]
            tree = Base.invokelatest(PackageTree, root, extensions, directory)
            Base.invokelatest(inventory, tree; options...)
        end
    end

    # Baseline tables that no guard owns.
    unknown_tables(document) = sort!([table for table in keys(document) if table ∉ TABLES])

    const LIVE = Dict{String, Any}()
    live() = isempty(LIVE) ? merge!(LIVE, inventory(live_tree())) : LIVE
    baseline(path = BASELINE) = TOML.parsefile(path)

    # Each counter of a file entry is compared on its own.
    function flatten(entries)
        flat = Dict{String, Int}()
        for (key, value) in entries
            if value isa AbstractDict
                for (counter, n) in value
                    n == 0 || (flat[string(key, " | ", counter)] = n)
                end
            else
                flat[key] = value
            end
        end
        return flat
    end

    # `added` lists unlisted keys and raised counts. `stale` lists listed keys
    # that are gone or whose count dropped.
    function compare(live, listed)
        actual, expected = flatten(live), flatten(listed)
        added = sort!([string(key, ": ", n, " (listed ", get(expected, key, 0), ")")
            for (key, n) in actual if n > get(expected, key, 0)])
        stale = sort!([string(key, ": ", get(actual, key, 0), " (listed ", n, ")")
            for (key, n) in expected if get(actual, key, 0) < n])
        return (; added, stale)
    end

    check(table) = compare(live()[table], get(baseline(), table, Dict{String, Any}()))

    toml_string(s) = "\"" * escape_string(s, "\"") * "\""

    function render(inventory)
        io = IOBuffer()
        println(io, "# Architecture violations present when the guards in")
        println(io, "# `test/quality/architecture.jl` were introduced. Delete or lower entries only.")
        println(io, "# Print the live inventory with `test/tools/architecture_inventory.jl`.")
        for table in TABLES
            haskey(inventory, table) || continue
            println(io, "\n[", table, "]")
            entries = inventory[table]
            for key in sort!(collect(keys(entries)))
                value = entries[key]
                if value isa AbstractDict
                    counts = join((string(c, " = ", get(value, c, 0)) for c in SWITCHES), ", ")
                    println(io, toml_string(key), " = { ", counts, " }")
                else
                    println(io, toml_string(key), " = ", value)
                end
            end
        end
        return String(take!(io))
    end
end

@testitem "Quality / architecture / ownership" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("ownership")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / placement" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("placement")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / direction" tags=[:quality] setup=[ArchitectureGuards] begin
    A = ArchitectureGuards
    # Each top-level submodule has a position in the declared order.
    top = Set(nameof(m) for m in A.live_tree().core
        if m !== LineCableModels && parentmodule(m) === LineCableModels)
    @test top == Set(A.ORDER)
    result = A.check("direction")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / names" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("names")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / shadowing" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("shadowing")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / validate returns its subject" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("validate")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / reserved verbs" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("reserved_verbs")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / symbol switches and probes" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("switches")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / negative controls of the architecture guards" tags=[:quality] setup=[ArchitectureGuards] begin
    A = ArchitectureGuards

    # The planted package has one or two violations for each guard. `Late` is
    # declared first in the source and has the second position in the order.
    planted_files = Dict(
        "src/ArchitectureProbePlanted.jl" => """
            module ArchitectureProbePlanted
            include("late/Late.jl")
            include("early/Early.jl")
            include("consumer/Consumer.jl")
            end
            """,
        "src/late/Late.jl" => """
            module Late
            marker() = 1
            end
            """,
        "src/early/Early.jl" => """
            module Early
            import ..Late
            struct Item end
            struct Formula end
            function extend end
            function quantity end
            upward() = Late.marker()
            end
            """,
        "src/consumer/Consumer.jl" => """
            module Consumer
            import ..Early
            export filter
            struct Local end
            struct Formula end
            function quantity end
            function filter end
            filter(::Local) = 1
            Early.extend(::Int) = 1
            Early.Item(::Int) = Early.Item()
            include("../early/stray.jl")
            include("../../ext/core.jl")
            end
            """,
        "src/early/stray.jl" => "stray() = 1\n",
        "ext/core.jl" => "core_method() = 1\n",
        "src/sources.jl" => raw"""
            validate(x::Int) = nothing
            validate(::Float64) = 1
            function validate(x::String)
                if isempty(x)
                    throw(ArgumentError("empty"))
                end
            end
            function validate(x::Symbol)
                x === :none && return
                return x
            end
            Owner.validate(x::Char) = normalize(x)
            function validate(x::Vector, context)
                validate(context)
                return x
            end
            _check_input(x) = x
            function ensure_ready end
            Owner.require_kind(x) = x
            __validate_shape(x) = x
            probe(f, x) = applicable(f, x) || Base.applicable(f, x)
            @eval generated() = 1
            Base.@eval generated_too() = 2
            switch(x) = x.kind === :a || x.kind == :b || x.kind in (:c, :d) || in(x.kind, [:e])
            negated(x) = x.kind !== :a && x.kind != :b && x.kind ∉ (:c, :d)
            """)
    planted_expected = Dict{String, Any}(
        "ownership" => Dict(
            "Consumer | Early.extend | src/consumer/Consumer.jl" => 1,
            "Consumer | Early.Item | src/consumer/Consumer.jl" => 1),
        "placement" => Dict("src/early/stray.jl | Consumer" => 1, "ext/core.jl | Consumer" => 1),
        "direction" => Dict("Early -> Late | src/early/Early.jl" => 1),
        "names" => Dict("quantity" => 2),
        "shadowing" => Dict("Consumer.filter -> Base" => 1),
        "validate" => Dict("src/sources.jl" => 6),
        "reserved_verbs" => Dict(
            "src/sources.jl | _check_input" => 1, "src/sources.jl | ensure_ready" => 1,
            "src/sources.jl | require_kind" => 1, "src/sources.jl | __validate_shape" => 1),
        "switches" => Dict("src/sources.jl" => Dict("applicable" => 2, "eval" => 2, "kind" => 7)))

    # A clean package with the same modules and one extension method.
    clean_files = Dict(
        "src/ArchitectureProbeClean.jl" => """
            module ArchitectureProbeClean
            include("early/Early.jl")
            include("late/Late.jl")
            end
            """,
        "src/early/Early.jl" => """
            module Early
            struct Item end
            struct Formula end
            function extend end
            function quantity end
            extend(::Item) = 1
            end
            """,
        "src/late/Late.jl" => """
            module Late
            import ..Early
            import Base: sum
            export sum
            struct Local end
            struct Formula end
            sum(::Local) = 0
            Early.extend(::Local) = 2
            Early.quantity(::Local) = 3
            Early.Item(::Local) = Early.Item()
            downward() = Early.extend(Early.Item())
            include("helpers.jl")
            end
            """,
        "src/late/helpers.jl" => "helper() = 1\n",
        "ext/ArchitectureProbeCleanExt.jl" => """
            module ArchitectureProbeCleanExt
            import ..ArchitectureProbeClean: Early
            Early.extend(::Int) = 4
            end
            """,
        "src/sources.jl" => raw"""
            validate(x::Int) = x
            function Owner.validate(x::Real, context)
                x > 0 || throw(ArgumentError("nonpositive"))
                foreach(context) do y
                    return y
                end
                return x
            end
            validate(x::Symbol) = x === :none ? error("none") : x
            function validate(x::String)
                if isempty(x)
                    throw(ArgumentError("empty"))
                elseif length(x) > 9
                    return x
                else
                    x
                end
            end
            function validate(x::Char)
                try
                    isvalid(x) || throw(ArgumentError("invalid"))
                    x
                catch
                    rethrow()
                end
            end
            validate(x::Float32, ::Type{Int}) = x
            function validate end
            @required Abstract begin
                validate(::Abstract)
            end
            validate!(x) = nothing
            validated(x) = nothing
            checker(x) = x
            check(x) = x
            text = "_check_input(x) = x, applicable(f, x), @eval f() = 1, x.kind === :a"
            # _check_input(x) = x, applicable(f, x), @eval f() = 1, x.kind === :a
            call() = _check_input(1)
            different(x, y) = x.kind === y || x.kind !== y || x.mode !== :a || kind === :a
            generate(x) = eval(x)
            """)

    planted = A.probe_inventory(planted_files, "ArchitectureProbePlanted";
        order = (:Early, :Late, :Consumer), dependencies = (Base,), tables = A.A_TABLES)
    clean = A.probe_inventory(clean_files, "ArchitectureProbeClean";
        extension = "ArchitectureProbeCleanExt", order = (:Early, :Late),
        dependencies = (Base,), tables = A.A_TABLES)
    for table in A.A_TABLES
        @testset "$table" begin
            @test planted[table] == planted_expected[table]
            @test isempty(clean[table])
        end
    end

    @testset "baseline comparison" begin
        listed = Dict{String, Any}("a" => 2, "b" => 1,
            "f" => Dict("applicable" => 1, "eval" => 0, "kind" => 2))
        @test A.compare(listed, listed) == (added = String[], stale = String[])
        raised = A.compare(Dict{String, Any}("a" => 3, "b" => 1, "c" => 1,
            "f" => Dict("applicable" => 1, "eval" => 1, "kind" => 2)), listed)
        @test raised.added == ["a: 3 (listed 2)", "c: 1 (listed 0)", "f | eval: 1 (listed 0)"]
        @test raised.stale == String[]
        lowered = A.compare(Dict{String, Any}("a" => 1,
            "f" => Dict("applicable" => 1, "eval" => 0, "kind" => 1)), listed)
        @test lowered.added == String[]
        @test lowered.stale == ["a: 1 (listed 2)", "b: 0 (listed 1)", "f | kind: 1 (listed 2)"]
        # The rendered baseline parses back to the inventory it was printed from.
        @test A.TOML.parse(A.render(planted)) == planted
        @test A.TOML.parse(A.render(clean)) == clean
    end
end

@testitem "Quality / architecture / Commons admission" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("commons")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / reserved vocabulary" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("vocabulary")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / literal fingerprints" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("fingerprints")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / clones of Commons bodies" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("clones")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / tiny private helpers" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("helpers")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / root freeze" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("root")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / baseline tables belong to guards" tags=[:quality] setup=[ArchitectureGuards] begin
    A = ArchitectureGuards
    @test A.unknown_tables(A.baseline()) == String[]
    @test A.unknown_tables(Dict("ownership" => Dict(), "retired" => Dict())) == ["retired"]
end

@testitem "Quality / architecture / negative controls of the Commons guards" tags=[:quality] setup=[ArchitectureGuards] begin
    A = ArchitectureGuards

    vocabulary = Dict{Symbol, Regex}(:shared => r"^_*shared$", :lonely => r"^_*lonely$",
        :undocumented => r"^_*undocumented$", :untested => r"^_*untested$",
        :ideal_transposition! => r"transpos", :inner_only => r"^_*inner_only$",
        :kernel => r"^_*kernel$", :scale_of => r"^_*scale_of$", :front => r"^_*front$",
        :vacuum_permittivity => A.VOCABULARY[:vacuum_permittivity],
        :vacuum_permeability => A.VOCABULARY[:vacuum_permeability])
    matrixops = """
        "Average a square matrix over its cyclic diagonals."
        function ideal_transposition!(matrix::AbstractMatrix)
            n = checksquare(matrix)
            coefficients = similar(diag(matrix))
            @inbounds for offset in 0:(n - 1)
                total = zero(eltype(matrix))
                for row in 1:n
                    total += matrix[row, 1 + mod(row - 1 + offset, n)]
                end
                coefficients[offset + 1] = total / n
            end
            @inbounds for row in 1:n, column in 1:n
                matrix[row, column] = coefficients[mod1(column - row + 1, n)]
            end
            return matrix
        end
        """

    # Commons has one definition failing each Commons admission criterion and `hidden`
    # failing all. `Alpha` and `Beta` plant violations of the reserved vocabulary, the
    # literal fingerprints, the clones and the tiny helpers. The root breaks the root freeze.
    planted_files = Dict(
        "src/CommonsProbePlanted.jl" => """
            module CommonsProbePlanted
            include("commons/Commons.jl")
            using .Commons: welcome
            export welcome
            include("alpha/Alpha.jl")
            include("beta/Beta.jl")
            root_helper() = 1
            const ROOT_LIMIT = 2
            end
            """,
        "src/commons/Commons.jl" => """
            module Commons
            export shared, lonely, undocumented, untested, unreserved, ideal_transposition!
            export welcome, extended, inner_only
            "Reached only through the closure of a private wrapper."
            inner_only(x) = x + 7
            hidden_wrapper(values) = map(value -> inner_only(value), values)
            "User API that no test names."
            welcome(x) = x
            "A hook extended by one owner."
            function extended end
            "Used by two owners."
            shared(x) = x + 1
            "Used by one owner."
            lonely(x) = x + 2
            undocumented(x) = x + 3
            "Named in no Commons test."
            untested(x) = x + 4
            "Reserves no vocabulary."
            unreserved(x) = x + 5
            hidden(x) = x + 6
            include("matrixops.jl")
            end
            """,
        "src/commons/matrixops.jl" => matrixops,
        "src/alpha/Alpha.jl" => """
            module Alpha
            import ..Commons
            uses(m) = Commons.shared(1) + Commons.lonely(1) + Commons.undocumented(1) +
                Commons.untested(1) + Commons.unreserved(1) + Commons.ideal_transposition!(m)[1] +
                Commons.hidden_wrapper([1])[1]
            const EPS0 = 8.8541878128e-12
            transposed!(m) = m
            function permeability(T)
                μ0 = 4π * 1e-7
                return T(μ0)
            end
            _tiny(x) = x + 1
            twice(x) = 2 * _tiny(x)
            Commons.extended(x::Int) = x
            end
            """,
        "src/beta/Beta.jl" => """
            module Beta
            import ..Commons
            uses(m) = Commons.shared(1) + Commons.undocumented(1) + Commons.untested(1) +
                Commons.unreserved(1) + Commons.ideal_transposition!(m)[1] +
                Commons.hidden_wrapper([1])[1]
            c0() = 299792458
            mu(T) = 4 * (one(T) * π) * (one(T) * 10)^(-7)
            function average_cyclic!(values::AbstractMatrix)
                size = checksquare(values)
                sums = similar(diag(values))
                @inbounds for shift in 0:(size - 1)
                    accumulator = zero(eltype(values))
                    for i in 1:size
                        accumulator += values[i, 1 + mod(i - 1 + shift, size)]
                    end
                    sums[shift + 1] = accumulator / size
                end
                @inbounds for i in 1:size, j in 1:size
                    values[i, j] = sums[mod1(j - i + 1, size)]
                end
                return values
            end
            _forward(x; scale = 1) = Commons.shared(x; scale = scale)
            first_use(x) = _forward(x) + 1
            second_use(x) = _forward(x; scale = 2) + 1
            end
            """,
        "test/unit/commons/commons.jl" =>
            "shared, lonely, undocumented, unreserved, ideal_transposition!, extended, inner_only\n",
        # Guard sources name definitions without testing them.
        "test/quality/guards.jl" => "_tiny, _forward\n")
    planted_expected = Dict{String, Any}(
        "commons" => Dict("lonely | owners" => 1, "undocumented | docstring" => 1,
            "untested | tests" => 1, "unreserved | vocabulary" => 1,
            "hidden | public" => 1, "hidden | docstring" => 1, "hidden | owners" => 1,
            "hidden | tests" => 1, "hidden | vocabulary" => 1,
            "welcome | tests" => 1, "extended | owners" => 1, "inner_only | owners" => 1,
            "hidden_wrapper | public" => 1, "hidden_wrapper | docstring" => 1,
            "hidden_wrapper | tests" => 1, "hidden_wrapper | vocabulary" => 1),
        "vocabulary" => Dict("src/alpha/Alpha.jl | EPS0" => 1,
            "src/alpha/Alpha.jl | transposed!" => 1, "src/alpha/Alpha.jl | μ0" => 1),
        "fingerprints" => Dict("src/alpha/Alpha.jl" => 2, "src/beta/Beta.jl" => 2),
        "clones" => Dict("src/beta/Beta.jl | average_cyclic! ~ ideal_transposition!" => 1),
        "helpers" => Dict("src/alpha/Alpha.jl" => 1, "src/beta/Beta.jl" => 1),
        "root" => Dict("CommonsProbePlanted" => 2))

    # Every Commons definition qualifies, and none of the near misses counts.
    clean_files = Dict(
        "src/CommonsProbeClean.jl" => """
            module CommonsProbeClean
            include("commons/Commons.jl")
            using .Commons: greet, declared
            export greet
            public declared
            include("alpha/Alpha.jl")
            include("beta/Beta.jl")
            end
            """,
        "src/commons/Commons.jl" => """
            module Commons
            export shared, vacuum_permittivity, ideal_transposition!, greet, accumulate_into!
            export front, kernel, scale_of
            public declared
            "Applied inside the closure of `front`."
            kernel(x) = x + 1
            "Applied inside the keyword body of `front`."
            scale_of(x) = 2x
            "Used by two owners. Its closure and keyword body use the other two."
            function front(values; factor = 1)
                weight = scale_of(factor)
                return map(value -> kernel(value) * weight, values)
            end
            "Used by two owners."
            shared(x) = x + 1
            "User API that no owner uses."
            greet(x) = x
            "User API that the root declares public."
            declared(x) = x
            "A hook whose default method a hook method of `Alpha` repeats."
            function accumulate_into!(values::AbstractVector)
                total = zero(eltype(values))
                for index in eachindex(values)
                    total += values[index]
                    values[index] = total
                end
                return values
            end
            include("consts.jl")
            include("matrixops.jl")
            end
            """,
        "src/commons/consts.jl" => """
            "Vacuum permittivity in F/m."
            vacuum_permittivity(::Type{T}) where {T} = T(88541878128) / T(10)^22
            """,
        "src/commons/matrixops.jl" => matrixops,
        "src/alpha/Alpha.jl" => """
            module Alpha
            import ..Commons
            function capacitance(m, T)
                ε0 = Commons.vacuum_permittivity(T)
                tolerance = 1e-7
                return Commons.shared(ε0) * Commons.ideal_transposition!(m)[1] * tolerance
            end
            phase(x) = 2π * x
            front_use(values) = Commons.front(values; factor = 2)
            "The vacuum permittivity is 8.8541878128e-12 F/m."
            documented(x) = x
            _tested(x) = x + 1
            calls_tested(x) = _tested(x) * 2
            function _long(x)
                a = x + 1
                b = a * 2
                c = b - 3
                return c
            end
            calls_long(x) = _long(x) + 1
            _reorder(x, y) = Commons.shared(y, x)
            first_reorder(x) = _reorder(x, 1)
            second_reorder(x) = _reorder(1, x)
            function Commons.accumulate_into!(values::Vector{Int})
                total = zero(eltype(values))
                for index in eachindex(values)
                    total += values[index]
                    values[index] = total
                end
                return values
            end
            end
            """,
        "src/beta/Beta.jl" => """
            module Beta
            import ..Commons
            admittance(m, T) = Commons.shared(1) * Commons.vacuum_permittivity(T) *
                Commons.ideal_transposition!(m)[1]
            Commons.accumulate_into!(values::Vector{Float64}) = values
            front_use(values) = Commons.front(values; factor = 3)
            end
            """,
        "test/unit/commons/commons.jl" =>
            "shared, vacuum_permittivity, ideal_transposition!, greet, declared, accumulate_into!, " *
            "front, kernel, scale_of\n",
        "test/unit/alpha.jl" => "_tested\n")

    planted = A.probe_inventory(planted_files, "CommonsProbePlanted"; vocabulary,
        tables = A.C_TABLES)
    clean = A.probe_inventory(clean_files, "CommonsProbeClean"; vocabulary,
        tables = A.C_TABLES)
    for table in A.C_TABLES
        @testset "$table" begin
            @test planted[table] == planted_expected[table]
            @test isempty(clean[table])
        end
    end
    @test A.TOML.parse(A.render(planted)) == planted

    # Throw guards are left out of the shingles. Other `||` statements stay in.
    body = A.JuliaSyntax.parseall(A.SyntaxNode, """
        function guarded(x)
            x > 0 || throw(ArgumentError("x must be positive"))
            isempty(x) && error("x is empty")
            x === nothing || isfinite(x) || Base.throw(DomainError(x))
            y = f(x) || g(x)
            return y
        end
        """)[1][2]
    @test A.body_tokens(body) == ["NAME", "=", "CALL", "(", "NAME", ")", "||", "CALL", "(",
        "NAME", ")", "return", "NAME"]
end

@testitem "Quality / architecture / baseline ratchet follows git renames" tags=[:quality] begin
    ratchet = Module(:BaselineRatchet)
    Base.include(ratchet, joinpath(pkgdir(LineCableModels), "test", "tools", "baseline_ratchet.jl"))
    mktempdir() do repository
        git(arguments...) = run(pipeline(Base.invokelatest(ratchet.git, repository,
            "-c", "user.name=Probe", "-c", "user.email=probe@example.invalid",
            arguments...); stdout = devnull, stderr = devnull))
        definitions = join(("f$i(x) = x + $i" for i in 1:20), "\n")
        function write_file(path, text)
            mkpath(dirname(joinpath(repository, path)))
            write(joinpath(repository, path), text)
        end
        baseline(rows...; table = "ownership") = write_file(ratchet.BASELINE,
            "[$table]\n" * join(("\"$key\" = $n" for (key, n) in rows), "\n") * "\n")
        grown() = Base.invokelatest(ratchet.grown, repository, "HEAD").lines

        git("init", "-q")
        write_file("src/early/Early.jl", "module Early\n$definitions\nend\n")
        write_file("src/early/part.jl", definitions * "\n")
        baseline("Early | Early.f1 | src/early/part.jl" => 2,
            "Late | Early.f2 | src/late/Late.jl" => 1)
        git("add", "-A")
        git("commit", "-q", "--no-verify", "--no-gpg-sign", "-m", "start")
        # The folder and the entry file move, and the module is renamed with them.
        git("mv", "src/early", "src/first")
        git("mv", "src/first/Early.jl", "src/first/First.jl")
        write_file("src/first/First.jl", "module First\n$definitions\nend\n")

        baseline("First | First.f1 | src/first/part.jl" => 2,
            "Late | First.f2 | src/late/Late.jl" => 1)
        @test grown() == String[]
        baseline("First | First.f1 | src/first/part.jl" => 3)
        @test grown() == ["ownership | First | First.f1 | src/first/part.jl: 3 (2 at HEAD)"]
        # Keys renamed without a matching git rename are added keys.
        baseline("First | First.f1 | src/other/part.jl" => 2,
            "Other | First.f2 | src/late/Late.jl" => 1)
        @test grown() == [
            "ownership | First | First.f1 | src/other/part.jl: 2 (absent at HEAD)",
            "ownership | Other | First.f2 | src/late/Late.jl: 1 (absent at HEAD)"]
        # A renamed file that is not a module entry file renames no module.
        git("mv", "src/first/part.jl", "src/first/Piece.jl")
        baseline("Piece | First.f1 | src/first/Piece.jl" => 2)
        @test grown() == ["ownership | Piece | First.f1 | src/first/Piece.jl: 2 (absent at HEAD)"]
        # A table absent at the reference belongs to a new guard.
        baseline("src/first/Piece.jl" => 1; table = "helpers")
        result = Base.invokelatest(ratchet.grown, repository, "HEAD")
        @test result.lines == String[]
        @test result.tables == ["helpers"]
    end
end
