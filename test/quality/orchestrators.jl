@testitem "Quality / declarative orchestrators / enforced contracts" tags = [:quality] begin
    using RequiredInterfaces
    using XLSX

    const Grammar = LineCableModels.Grammar
    const PlotBuilder = LineCableModels.PlotBuilder
    const ReportBuilder = LineCableModels.ReportBuilder

    @test parentmodule(Grammar.orchestrator_root) === Grammar
    @test parentmodule(Grammar.orchestrator_method) === Grammar
    @test Base.ispublic(Grammar, :orchestrator_root)
    @test Base.ispublic(Grammar, :orchestrator_method)
    @test !isdefined(LineCableModels, :orchestrator_root)
    @test !isdefined(LineCableModels, :orchestrator_method)

    function signature_mentions(value, owner)
        value === owner && return true
        value isa TypeVar && return signature_mentions(value.lb, owner) ||
               signature_mentions(value.ub, owner)
        value isa UnionAll && return signature_mentions(value.var, owner) ||
               signature_mentions(value.body, owner)
        value isa DataType && return any(
            parameter -> signature_mentions(parameter, owner),
            value.parameters
        )
        value isa Union && return signature_mentions(value.a, owner) ||
               signature_mentions(value.b, owner)
        return false
    end

    function referenced_functions(method::Method)
        found = Set{Any}()
        visited = Set{Method}()
        function visit_method(candidate)
            candidate in visited && return nothing
            push!(visited, candidate)
            foreach(visit, Base.uncompressed_ast(candidate).code)
            return nothing
        end
        function visit(value)
            if value isa GlobalRef
                if isdefined(value.mod, value.name)
                    resolved = getproperty(value.mod, value.name)
                    if resolved isa Function
                        push!(found, resolved)
                        # Keyword methods use one compiler-generated function.
                        # Follow that local wrapper to inspect the choreography
                        # recorded in the declared public method.
                        if value.mod === method.module &&
                           startswith(String(value.name), "#")
                            foreach(visit_method, methods(resolved))
                        end
                    end
                end
            elseif value isa Expr
                foreach(visit, value.args)
            elseif value isa QuoteNode
                visit(value.value)
            end
            return nothing
        end
        visit_method(method)
        return found
    end

    function test_orchestrator(
            action,
            root,
            owned_implementors;
            required_hooks = nothing,
            optional_hooks = ()
    )
        findings = String[]
        hooks = if required_hooks === nothing
            RequiredInterfaces.isInterface(root) ?
            Tuple(RequiredInterfaces.functions(
                RequiredInterfaces.getInterface(root),
            )) : ()
        else
            required_hooks
        end
        declared_root = try
            Grammar.orchestrator_root(action, root)
        catch error
            push!(findings, "missing root metadata: $(sprint(showerror, error))")
            nothing
        end
        declared_root === root || push!(
            findings,
            "declared root is $declared_root instead of $root"
        )
        root isa Type && isabstracttype(root) || push!(
            findings,
            "declared owner $root is not an abstract type"
        )
        parentmodule(root) === parentmodule(action) || push!(
            findings,
            "declared owner $root does not belong to $(parentmodule(action))"
        )

        declared_method = try
            Grammar.orchestrator_method(action, root)
        catch error
            push!(findings, "missing method metadata: $(sprint(showerror, error))")
            nothing
        end
        if declared_method isa Method
            declared_method.module === parentmodule(action) || push!(
                findings,
                "declared method belongs to $(declared_method.module)"
            )
            declared_method.line > 0 || push!(
                findings,
                "declared method has no source line"
            )
            signature_mentions(declared_method.sig, root) || push!(
                findings,
                "declared method signature $(declared_method.sig) does not bind $root"
            )
            specializations = Method[method
                                     for method in methods(action)
                                     if method !== declared_method &&
                Base.morespecific(method.sig, declared_method.sig)]
            isempty(specializations) || push!(
                findings,
                "more-specific public methods: $(join(specializations, "; "))"
            )

            referenced = referenced_functions(declared_method)
            missing_hooks = Tuple(
                hook for hook in (hooks..., optional_hooks...)
            if hook ∉ referenced
            )
            isempty(missing_hooks) || push!(
                findings,
                "declared method does not call hooks: $(join(missing_hooks, ", "))"
            )
        end

        invalid_implementors = Tuple(
            implementor for implementor in owned_implementors
        if !(implementor <: root)
        )
        isempty(invalid_implementors) || push!(
            findings,
            "owned implementors outside $root: $(join(invalid_implementors, ", "))"
        )
        missing_implementations = Tuple(
            (implementor, hook) for implementor in owned_implementors
        for hook in hooks
        if !any(method -> signature_mentions(method.sig, implementor), methods(hook))
        )
        isempty(missing_implementations) || push!(
            findings,
            "missing required hook implementations: " * join(
                ("$(implementor).$(nameof(hook))"
                for
                (implementor, hook) in missing_implementations),
                ", "
            )
        )

        @testset "$(nameof(action)) with $root" begin
            @test isempty(findings)
        end
        return findings
    end

    plot_implementors = Tuple(RequiredInterfaces.nonabstract_subtypes(
        PlotBuilder.AbstractPlotDefinition,
    ))
    plot_findings = test_orchestrator(
        PlotBuilder.make_render,
        PlotBuilder.AbstractPlotDefinition,
        plot_implementors;
        required_hooks = (PlotBuilder.entitle, PlotBuilder.resolve, PlotBuilder.fetch),
        optional_hooks = (PlotBuilder.parse, PlotBuilder.finish)
    )

    report_implementors = Tuple(RequiredInterfaces.nonabstract_subtypes(
        ReportBuilder.AbstractReportDefinition,
    ))
    report_findings = test_orchestrator(
        ReportBuilder.report,
        ReportBuilder.AbstractReportDefinition,
        report_implementors;
        optional_hooks = (ReportBuilder.finish,)
    )

    # Observation publication is one closed Grammar-owned method, but its
    # sources deliberately have no common abstract owner. `observe` is the
    # open result hook. Treating either generic as an AbstractProblemResult
    # orchestrator would exclude standalone matrix wrappers and external
    # result types, so CI guards the sole publication method directly.
    publication_methods = Method[method
                                 for method in methods(LineCableModels.observables)
                                 if method.nargs == 3]
    @test length(publication_methods) == 1
    publication_method = only(publication_methods)
    @test publication_method.module === Grammar
    @test publication_method.sig ==
          Tuple{typeof(LineCableModels.observables), Any, Tuple}
    @test_throws MethodError Grammar.orchestrator_root(
        LineCableModels.observables,
        Grammar.AbstractProblemResult
    )
    @test_throws MethodError Grammar.orchestrator_root(
        LineCableModels.observe,
        Grammar.AbstractProblemResult
    )

    @test plot_findings isa Vector{String}
    @test report_findings isa Vector{String}
end
