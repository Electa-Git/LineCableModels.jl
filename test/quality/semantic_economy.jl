@testmodule SemanticEconomyChecks begin
    using LineCableModels

    module ForwardingProbe
        _kernel(x)=2x
        _duplicate(x)=_kernel(x)
        _dispatch(x::Number)=_kernel(x)
        _dispatch(x::Tuple)=map(_kernel, x)
    end

    function module_tree(root::Module)
        owners=Module[root]
        for owner in owners, name in names(owner; all=true)
            isdefined(owner, name) || continue
            value=getfield(owner, name)
            value isa Module && value !== owner && parentmodule(value) === owner &&
                value ∉ owners && push!(owners, value)
        end
        return owners
    end

    function referenced_globals(value)
        found=GlobalRef[]
        function visit(node)
            if node isa GlobalRef
                push!(found, node)
            elseif node isa Expr
                foreach(visit, node.args)
            elseif node isa Core.CodeInfo
                foreach(visit, node.code)
            end
        end
        visit(value)
        return found
    end

    # Inspect actual bindings, including aliases. An underscore on a local
    # numerical kernel is fine; referring to another owner's private binding is not.
    function private_edges(owner, code, owners)
        violations=Tuple{Module, Symbol}[]
        for reference in referenced_globals(code)
            isdefined(reference.mod, reference.name) || continue
            value=getfield(reference.mod, reference.name)
            value isa Union{Function, DataType, UnionAll} || continue
            value isa UnionAll && !(Base.unwrap_unionall(value) isa DataType) && continue
            target=parentmodule(value)
            name=nameof(value)
            target in owners && target !== owner && startswith(String(name), "_") &&
                push!(violations, (target, name))
        end
        return unique(violations)
    end

    function forwarded_argument(value)
        value isa Symbol && return value
        value isa Expr || return nothing
        value.head in (:(::), :kw) && return forwarded_argument(first(value.args))
        value.head === :(...) && return Expr(:(...), forwarded_argument(only(value.args)))
        if value.head === :parameters
            return Expr(:parameters, map(value.args) do argument
                name=forwarded_argument(argument)
                name isa Symbol ? Expr(:kw, name, name) : name
            end...)
        end
        return nothing
    end

    function exact_forwarder(expression)
        expression isa Expr && expression.head in (:function, :(=)) || return nothing
        length(expression.args) == 2 || return nothing
        signature, body=expression.args
        while signature isa Expr && signature.head in (:where, :(::))
            signature=first(signature.args)
        end
        signature isa Expr && signature.head === :call || return nothing
        name=first(signature.args)
        name isa Symbol && startswith(String(name), "_") || return nothing
        if body isa Expr && body.head === :block
            statements=filter(value -> !(value isa LineNumberNode), body.args)
            length(statements) == 1 || return nothing
            body=only(statements)
        end
        body isa Expr && body.head === :return && (body=only(body.args))
        body isa Expr && body.head === :call || return nothing
        expected=forwarded_argument.(signature.args[2:end])
        any(isnothing, expected) && return nothing
        actual=map(body.args[2:end]) do argument
            if argument isa Expr && argument.head === :parameters
                Expr(:parameters, map(argument.args) do keyword
                    keyword isa Symbol ? Expr(:kw, keyword, keyword) : keyword
                end...)
            else
                argument
            end
        end
        expected == actual || return nothing
        return (name, first(body.args))
    end

    function binding(owner, expression)
        expression isa Symbol && return isdefined(owner, expression) ? getfield(owner, expression) : nothing
        expression isa Expr && expression.head === :. && length(expression.args) == 2 || return nothing
        parent=binding(owner, first(expression.args))
        name=last(expression.args)
        name=name isa QuoteNode ? name.value : name
        return parent isa Module && name isa Symbol && isdefined(parent, name) ?
               getfield(parent, name) : nothing
    end

    function inspect(roots)
        owners=unique(vcat(module_tree.(roots)...))
        definitions=Set{Method}()
        private=String[]
        files=Dict{String, Set{Module}}()
        # Base methods authored by an owner (show, getindex, iterate, ...) count too.
        for owner in [Base; Core; owners], name in names(owner; all=true, imported=true)
            isdefined(owner, name) || continue
            value=getfield(owner, name)
            value isa Union{Function, DataType, UnionAll} || continue
            value isa UnionAll && !(Base.unwrap_unionall(value) isa DataType) && continue
            for method in methods(value, owners)
                method in definitions && continue
                push!(definitions, method)
                path=abspath(String(method.file))
                isfile(path) || continue
                push!(get!(Set{Module}, files, path), method.module)
                if startswith(String(nameof(value)), "_") && parentmodule(value) !== method.module
                    push!(private, "$path:$(method.line): $(method.module) extends private $value")
                end
                for (target, symbol) in private_edges(method.module, Base.uncompressed_ast(method), owners)
                    push!(private, "$path:$(method.line): $(method.module) uses $target.$symbol")
                end
            end
        end
        declarations=Dict{Function, Vector{Tuple{String, Any}}}()
        for (path, namespaces) in files
            # A shared macro-definition file can emit methods in several owners;
            # those generated methods were checked above, not guessed lexically.
            length(namespaces) == 1 || continue
            owner=only(namespaces)
            function visit(expression)
                expression isa Expr || return
                if expression.head in (:function, :(=)) && length(expression.args) == 2
                    signature=first(expression.args)
                    while signature isa Expr && signature.head in (:where, :(::))
                        signature=first(signature.args)
                    end
                    if signature isa Expr && signature.head === :call
                        name=first(signature.args)
                        value=binding(owner, name)
                        if name isa Symbol && startswith(String(name), "_") && value isa Function
                            found=exact_forwarder(expression)
                            target=found === nothing ? nothing : binding(owner, last(found))
                            push!(get!(Vector{Tuple{String, Any}}, declarations, value), (path, target))
                        end
                    end
                end
                # Local functions are not module-level semantic authorities.
                expression.head === :function && return
                expression.head === :(=) && first(expression.args) isa Expr &&
                    first(expression.args).head in (:call, :where) && return
                foreach(visit, expression.args)
            end
            visit(Meta.parseall(read(path, String)))
        end
        forwarders=String[]
        for (wrapper, bodies) in declarations
            target=last(first(bodies))
            # A dispatch branch selecting a kernel is not a duplicate authority.
            # Reject a wrapper only when all its declarations pass every argument
            # unchanged to the very same package-owned function.
            target isa Function && parentmodule(target) in owners || continue
            all(body -> last(body) === target, bodies) || continue
            for (path, _) in bodies
                push!(forwarders, "$path: $wrapper forwards unchanged to $target")
            end
        end
        return (; private=sort!(unique(private)), forwarders=sort!(unique(forwarders)))
    end
end

@testitem "Quality / semantic economy / owned calls and forwarding" tags=[:quality] setup=[SemanticEconomyChecks] begin
    using Gmsh
    using Measurements
    using Calculus
    using Distributions
    using XLSX

    roots=Module[LineCableModels]
    for name in (:LineCableModelsGmshExt, :LineCableModelsMeasurementsExt,
            :LineCableModelsDistributionsExt, :LineCableModelsXLSXExt)
        extension=Base.get_extension(LineCableModels, name)
        @test extension !== nothing
        extension === nothing || push!(roots, extension)
    end
    found=SemanticEconomyChecks.inspect(roots)
    @test isempty(found.private)
    @test isempty(found.forwarders)

    # Prove the detector catches fresh names, keywords and splats, without
    # banning local kernels, argument adaptation, or arithmetic.
    for source in ("_alias(x) = target(x)",
            "_alias(x; y=1) = target(x; y=y)",
            "_alias(args...; kwargs...) = target(args...; kwargs...)",
            "function _alias(x::T) where T; return target(x); end")
        @test SemanticEconomyChecks.exact_forwarder(Meta.parse(source)) !== nothing
    end
    for source in ("_kernel(x) = target(2x)", "_adapt(x) = target(x, 1)",
            "_kernel(x) = x * x", "public_action(x) = target(x)")
        @test SemanticEconomyChecks.exact_forwarder(Meta.parse(source)) === nothing
    end
    owner=LineCableModels.Engine
    other=LineCableModels.ReportBuilder
    probe=Expr(:call, GlobalRef(other, :_is_diagonal), :table, :columns)
    @test SemanticEconomyChecks.private_edges(owner, probe, [owner, other]) ==
          [(other, :_is_diagonal)]
    @test isempty(SemanticEconomyChecks.private_edges(other, probe, [owner, other]))

    probe_results=SemanticEconomyChecks.inspect([SemanticEconomyChecks.ForwardingProbe])
    @test isempty(probe_results.private)
    @test length(probe_results.forwarders) == 1
    @test occursin("_duplicate", only(probe_results.forwarders))
end

@testitem "Makie addons / semantic economy / owned calls and forwarding" tags=[:visual] setup=[SemanticEconomyChecks] begin
    using LineCableModels
    using CairoMakie

    roots=Module[LineCableModels,
        Base.get_extension(LineCableModels, :LineCableModelsMakieExt),
        Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt)]
    found=SemanticEconomyChecks.inspect(roots)
    @test isempty(found.private)
    @test isempty(found.forwarders)
end
