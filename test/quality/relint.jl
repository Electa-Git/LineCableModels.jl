@testitem "Quality / ReLint / advisory source patterns" tags=[:quality] default_imports=false begin
    using Test, ReLint, Argus, JuliaSyntax
    using JuliaSyntax: children, numchildren, is_leaf, is_prefix_call

    # Resolve spelling only, never the ownership of a source binding.
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

    function execution_context(node)
        ancestors = typeof(node)[]
        parent = node.parent
        while parent !== nothing
            push!(ancestors, parent)
            parent = parent.parent
        end
        quoted, in_function, opaque, in_type = 0, false, false, false
        for i in length(ancestors):-1:1
            ancestor = ancestors[i]
            child = i == 1 ? node : ancestors[i - 1]
            k = kind(ancestor)
            if k == K"quote"
                quoted += 1
            elseif k == K"$" && quoted > 0
                quoted -= 1
            elseif quoted == 0
                if k in (K"function", K"->", K"do") && numchildren(ancestor) >= 2
                    # Signatures and default expressions are outside this rule's scope.
                    in_function |= child === children(ancestor)[end]
                elseif k == K"macrocall"
                    # Docstrings wrap declarations without generating their bodies.
                    name = syntax_name(ancestor[1])
                    opaque |= name ∉ ((Symbol("@doc"),), (:Core, Symbol("@doc")))
                elseif k == K"struct"
                    in_type = true
                end
            end
        end
        return (; quoted=quoted > 0, in_function, opaque, in_type)
    end

    function runtime_eval(node)
        k = kind(node)
        k in (K"call", K"macrocall") && numchildren(node) > 0 || return false
        k == K"call" && !is_prefix_call(node) && return false
        syntax_name(node[1]) in ((:eval,), (:Base, :eval), (:Core, :eval),
            (Symbol("@eval"),), (:Base, Symbol("@eval")), (:Core, Symbol("@eval"))) || return false
        context = execution_context(node)
        return context.in_function && !context.quoted && !context.opaque
    end

    function constant_catch_result(node)
        kind(node) == K"catch" || return false
        context = execution_context(node)
        (context.quoted || context.opaque) && return false
        body = children(node)[end]
        while kind(body) == K"block"
            numchildren(body) == 0 && return true
            numchildren(body) == 1 || return false
            body = body[1]
        end
        if kind(body) == K"return"
            numchildren(body) == 0 && return true
            numchildren(body) == 1 || return false
            body = body[1]
        end
        syntax_name(body) == (:nothing,) && return true
        is_leaf(body) && body.val isa Union{Bool,Number} && return true
        # Negative numeric literals have a unary operator node in Julia syntax.
        return kind(body) == K"call" && numchildren(body) == 2 &&
            syntax_name(body[1]) in ((:+,), (:-,)) &&
            is_leaf(body[2]) && body[2].val isa Number
    end

    # Only declarations visible in the enclosing source scope identify local
    # constructors. No package code is loaded to resolve names or aliases.
    function declares_type(node, name)
        if kind(node) == K"macrocall" &&
                syntax_name(node[1]) in ((Symbol("@doc"),), (:Core, Symbol("@doc")))
            return declares_type(children(node)[end], name)
        end
        if kind(node) in (K"struct", K"abstract", K"primitive")
            declaration = node[1]
            while kind(declaration) in (K"curly", K"<:")
                declaration = declaration[1]
            end
            return syntax_name(declaration) == (name,)
        end
        kind(node) in (K"toplevel", K"block", K"if", K"elseif", K"else") || return false
        return any(child -> declares_type(child, name), children(node))
    end

    function private_forwarder(node)
        kind(node) == K"function" && numchildren(node) == 2 || return false
        context = execution_context(node)
        (context.quoted || context.opaque || context.in_type) && return false
        signature = node[1]
        while kind(signature) == K"where"
            signature = signature[1]
        end
        kind(signature) == K"call" && is_prefix_call(signature) || return false
        name = syntax_name(signature[1])
        length(name) == 1 && startswith(String(only(name)), "_") || return false
        bindings = Symbol[]
        for argument in children(signature)[2:end]
            if kind(argument) == K"::" && numchildren(argument) == 2
                argument = argument[1]
            end
            kind(argument) == K"Identifier" || return false
            all(==('_'), String(argument.val)) && return false
            push!(bindings, argument.val)
        end
        allunique(bindings) || return false
        body = node[2]
        while kind(body) == K"block" && numchildren(body) == 1
            body = body[1]
        end
        kind(body) == K"return" && numchildren(body) == 1 && (body = body[1])
        kind(body) == K"call" && is_prefix_call(body) || return false
        target = syntax_name(body[1])
        isempty(target) && return false
        last(target) == only(name) && return false
        # This pattern excludes known conversion spellings and built-in type
        # constructors without resolving arbitrary source aliases.
        last(target) in (:convert, :oftype, :reinterpret, :cconvert, :unsafe_convert) && return false
        if length(target) == 1 || first(target) in (:Base, :Core)
            any(owner -> isdefined(owner, last(target)) && getfield(owner, last(target)) isa Type,
                (Base, Core)) && return false
        end
        scope = node.parent
        while scope.parent !== nothing && kind(scope.parent) != K"module"
            scope = scope.parent
        end
        (declares_type(scope, only(name)) || declares_type(scope, last(target))) && return false
        forwarded = children(body)[2:end]
        return length(forwarded) == length(bindings) &&
            all(pair -> kind(pair[1]) == K"Identifier" && pair[1].val == pair[2],
                zip(forwarded, bindings))
    end

    # Argus compiles conditions in its own evaluation context. Explicit references
    # keep the predicates in this test item's module, without registering globals
    # or calling dependency internals. ReLint supplies parsing and rule traversal.
    rules = [
        Argus.Rule("runtime-eval",
            "Runtime evaluation candidate; spelling does not establish binding ownership.",
            Argus.Pattern(Argus.SyntaxPatternNode(:(~and({node},
                ~when([:node], $(GlobalRef(@__MODULE__, :runtime_eval))(node.src))))))),
        Argus.Rule("constant-catch-result",
            "Swallowed-exception candidate: empty or constant-result catch.",
            Argus.Pattern(Argus.SyntaxPatternNode(:(~and({node},
                ~when([:node], $(GlobalRef(@__MODULE__, :constant_catch_result))(node.src))))))),
        Argus.Rule("private-forwarder",
            "Unchanged forwarding candidate; dispatch purpose has not been determined.",
            Argus.Pattern(Argus.SyntaxPatternNode(:(~and({node},
                ~when([:node], $(GlobalRef(@__MODULE__, :private_forwarder))(node.src))))))),
    ]
    context = ReLint.LintContext(rules)

    @testset "ReLint rule controls" begin
        examples = [
            ("runtime-eval", true, "function f(x); eval(x); end"),
            ("runtime-eval", true, "f(x) = Base.eval(@__MODULE__, x)"),
            ("runtime-eval", true, "f(x) = Core.eval(@__MODULE__, x)"),
            ("runtime-eval", true, "f(x) = @eval generated() = 1"),
            ("runtime-eval", true, "f = x -> eval(x)"),
            ("runtime-eval", true, "map(xs) do x; eval(x); end"),
            ("runtime-eval", true, raw"f(x) = :($(eval(x)))"),
            ("runtime-eval", false, "@eval generated(x) = x"),
            ("runtime-eval", false, "f(x) = :(eval(x))"),
            ("runtime-eval", false, "quote; f(x) = eval(x); end"),
            ("runtime-eval", false, "f(x) = Other.eval(x)"),
            ("runtime-eval", false, "f(x) = evaluator(x)"),
            ("constant-catch-result", true, "f() = try work() catch; end"),
            ("constant-catch-result", true, "function f(); try work() catch err; return nothing; end; end"),
            ("constant-catch-result", true, "f() = try work() catch; false; end"),
            ("constant-catch-result", true, "f() = try work() catch; return -1; end"),
            ("constant-catch-result", false, "f() = try work() catch; rethrow(); end"),
            ("constant-catch-result", false, "f() = try work() catch; cleanup(); nothing; end"),
            ("constant-catch-result", false, "f() = try work() catch; fallback(); end"),
            ("constant-catch-result", false, "quote; try work() catch; nothing; end; end"),
            ("constant-catch-result", false, "stage(::Nothing) = nothing"),
            ("private-forwarder", true, "_f(x, y) = g(x, y)"),
            ("private-forwarder", true, "function _f(x::T, y) where T; return g(x, y); end"),
            ("private-forwarder", true, "_f(α::Real) = Owner.g(α)"),
            ("private-forwarder", false, "Base.show(io::IO, x::T) = render(io, x)"),
            ("private-forwarder", false, "_f(x) = g(convert(Float64, x))"),
            ("private-forwarder", false, "_f(x) = Float64(x)"),
            ("private-forwarder", false, "_f(x)::Float64 = g(x)"),
            ("private-forwarder", false, "_f(kind::Symbol, x) = _f(Val(kind), x)"),
            ("private-forwarder", false, "_f(::Val{:known}, x) = g(x)"),
            ("private-forwarder", false, "_f(x, y) = g(y, x)"),
            ("private-forwarder", false, "_f(x=1) = g(x)"),
            ("private-forwarder", false, "_f(x; y=1) = g(x; y=y)"),
            ("private-forwarder", false, "_f(xs...) = g(xs...)"),
            ("private-forwarder", false, "_f((x, y)) = g(x, y)"),
            ("private-forwarder", false, "struct _T; x; _T(x) = new(x); end"),
            ("private-forwarder", false, "struct _T; x; end; _T(x) = g(x)"),
            ("private-forwarder", false, "\"A type.\"\nstruct _T; x; end; _T(x) = g(x)"),
            ("private-forwarder", false, "struct _T{T}; x::T; end; _T{T}(x) where T = g(x)"),
            ("private-forwarder", false, "struct LocalType; x; end; _f(x) = LocalType(x)"),
            ("private-forwarder", false, "@inline _f(x) = g(x)"),
            ("private-forwarder", false, "quote; _f(x) = g(x); end"),
            ("private-forwarder", true, "\"\"\"Value \$(1).\"\"\"\n_f(α) = g(α)"),
        ]
        for (rule, expected, source) in examples
            @testset "$rule: $source" begin
                findings = ReLint.lint_text(source; context)
                @test count(f -> f.rule.name == rule, findings) == Int(expected)
            end
        end
        unicode = ReLint.lint_text("_α(β) = γ(β)"; context, filename="unicode.jl")
        @test only(unicode).file == "unicode.jl"
        @test (only(unicode).line, only(unicode).column) == (1, 1)
        @test_throws JuliaSyntax.ParseError ReLint.lint_text("function broken("; context)
        @test_throws SystemError ReLint.lint_file(joinpath(tempdir(), "absent-relint-$(getpid()).jl"), context)
    end

    @testset "ReLint production scan" begin
        root = normpath(joinpath(@__DIR__, "..", ".."))
        files = String[]
        for name in ("src", "ext")
            directory = joinpath(root, name)
            isdir(directory) || error("Missing production source root: $directory")
            for (path, _, names) in walkdir(directory), name in names
                endswith(name, ".jl") && push!(files, joinpath(path, name))
            end
        end
        sort!(files)
        isempty(files) && error("Empty production source inventory")
        println("ReLint ", pkgversion(ReLint), "; Argus ", pkgversion(Argus),
            "; JuliaSyntax ", pkgversion(JuliaSyntax), "; requested source files: ", length(files))
        counts = Dict(rule.name => 0 for rule in rules)
        failures = 0
        for file in files
            try
                for finding in ReLint.lint_file(file, context)
                    println("ReLint advisory [", finding.rule.name, "] ", relpath(file, root),
                        ":", finding.line, ":", finding.column, ": ", finding.msg)
                    counts[finding.rule.name] += 1
                end
            catch exception
                failures += 1
                println("ReLint scan failure ", relpath(file, root), ": ", sprint(showerror, exception))
            end
        end
        for rule in sort(collect(keys(counts)))
            println("ReLint advisory count [", rule, "]: ", counts[rule])
        end
        println("ReLint production scan failures: ", failures)
        @test failures == 0
    end
end
