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
    # descendants. The root module and package extensions are exempt.
    const ORDER = (:Units, :Grammar, :TextDisplay, :InputValidation, :PlotBuilder,
        :Materials, :Earth, :DataModel, :Engine, :ModalAnalysis, :ParametricBuilder,
        :UQ, :ReportBuilder, :ImportExport, :PSCAD)
    const EXTENSIONS = (:LineCableModelsMeasurementsExt, :LineCableModelsDistributionsExt,
        :LineCableModelsGmshExt, :LineCableModelsXLSXExt, :LineCableModelsMakieExt,
        :LineCableModelsCairoMakieExt)
    const DEPENDENCIES = (Base, DataFrames, LinearAlgebra, Statistics, Random,
        Base.require(LineCableModels, :Dates), Logging)
    # Per-family formula selections are the only name shared across modules.
    const SHARED_NAMES = (:Formula,)
    const RESERVED_VERB = r"^_*(validate|check|require|assert|verify|ensure)_"
    const TABLES = ("ownership", "placement", "direction", "names", "shadowing",
        "validate", "reserved_verbs", "switches")
    const SWITCHES = ("applicable", "eval", "kind")
    const BASELINE = joinpath(@__DIR__, "architecture_baseline.toml")

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
    relative(path, directory) = join(splitpath(relpath(path, directory)), "/")

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
            root = joinpath(directory, base)
            isdir(root) || continue
            for (path, _, names) in walkdir(root), name in names
                endswith(name, ".jl") && push!(files, joinpath(path, name))
            end
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

    function mentioned!(found::Set{DataType}, T)
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
            foreach(p -> mentioned!(found, p), T.parameters)
        end
        return found
    end

    # A1. A core method that extends a function owned by another package
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
            types = Set{DataType}()
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

    # A2. No core method is defined under `ext/`. The nearest home directory
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
            nearest = nothing
            for directory in keys(homes)
                inside(file, directory) || continue
                (nearest === nothing || length(directory) > length(nearest)) &&
                    (nearest = directory)
            end
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

    # A3. No method body of a core submodule references a value owned by a
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

    # A4a. Distinct functions or types, owned by different package modules,
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

    # A4b. A public package-owned name that a dependency exports for a
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

    function validate_returns_subject(node)
        name = first_argument(signature_call(node))
        name === nothing && return false
        return returns_subject(node[2], name) && early_returns_subject(node[2], name)
    end

    is_required_block(node) = kind(node) == K"macrocall" && numchildren(node) > 0 &&
        syntax_name(node[1]) in ((Symbol("@required"),),
            (:RequiredInterfaces, Symbol("@required")))

    # A5. Every `validate` method names its first positional argument and
    # returns it on every path, or throws.
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

    # A6. No function name spells an input check with a reserved verb prefix.
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

    # A7. Per file, calls to `applicable`, `@eval` calls and comparisons of a
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

    function inventory(tree::PackageTree; order = ORDER, dependencies = DEPENDENCIES)
        methods = package_methods(tree)
        files = source_files(tree.directory)
        return Dict{String, Any}(
            "ownership" => ownership(methods, tree),
            "placement" => placement(methods, tree, home_directories(tree)),
            "direction" => direction(methods, tree, order),
            "names" => shared_names(tree),
            "shadowing" => shadowing(tree, dependencies),
            "validate" => validate_returns(files, tree.directory),
            "reserved_verbs" => reserved_verbs(files, tree.directory),
            "switches" => switches(files, tree.directory))
    end

    function live_tree()
        extensions = map(EXTENSIONS) do name
            extension = Base.get_extension(LineCableModels, name)
            extension === nothing && error("Extension $name is not loaded")
            extension
        end
        return PackageTree(LineCableModels, extensions, pkgdir(LineCableModels))
    end

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

@testitem "Quality / architecture / A1 ownership" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("ownership")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / A2 placement" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("placement")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / A3 direction" tags=[:quality] setup=[ArchitectureGuards] begin
    A = ArchitectureGuards
    # Each top-level submodule has a position in the declared order.
    top = Set(nameof(m) for m in A.live_tree().core
        if m !== LineCableModels && parentmodule(m) === LineCableModels)
    @test top == Set(A.ORDER)
    result = A.check("direction")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / A4 names" tags=[:quality] setup=[ArchitectureGuards] begin
    for table in ("names", "shadowing")
        result = ArchitectureGuards.check(table)
        @test result.added == String[]
        @test result.stale == String[]
    end
end

@testitem "Quality / architecture / A5 validate returns its subject" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("validate")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / A6 reserved verbs" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("reserved_verbs")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / A7 symbol switches and probes" tags=[:quality] setup=[ArchitectureGuards] begin
    result = ArchitectureGuards.check("switches")
    @test result.added == String[]
    @test result.stale == String[]
end

@testitem "Quality / architecture / A8 negative controls" tags=[:quality] setup=[ArchitectureGuards] begin
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
        "validate" => Dict("src/sources.jl" => 5),
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
                    x
                catch
                    rethrow()
                end
            end
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

    function probe_inventory(files, name, extension, order)
        mktempdir() do directory
            for (path, text) in files
                mkpath(dirname(joinpath(directory, path)))
                write(joinpath(directory, path), text)
            end
            root = Base.include(Main, joinpath(directory, "src", name * ".jl"))
            extensions = extension === nothing ? Module[] :
                [Base.include(Main, joinpath(directory, "ext", extension * ".jl"))]
            tree = Base.invokelatest(A.PackageTree, root, extensions, directory)
            Base.invokelatest(A.inventory, tree; order, dependencies = (Base,))
        end
    end
    planted = probe_inventory(planted_files, "ArchitectureProbePlanted", nothing,
        (:Early, :Late, :Consumer))
    clean = probe_inventory(clean_files, "ArchitectureProbeClean",
        "ArchitectureProbeCleanExt", (:Early, :Late))
    for table in A.TABLES
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
