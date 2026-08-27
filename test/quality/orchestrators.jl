@testitem "Quality / declarative orchestrators / advisory contracts" tags = [:quality] begin
    using RequiredInterfaces

    const Grammar = LineCableModels.Grammar
    const PlotBuilder = LineCableModels.PlotBuilder
    const ReportBuilder = LineCableModels.ReportBuilder
    const ParametricBuilder = LineCableModels.ParametricBuilder

    @test parentmodule(Grammar.orchestrator_root) === Grammar
    @test parentmodule(Grammar.orchestrator_method) === Grammar
    @test Base.ispublic(Grammar, :orchestrator_root)
    @test Base.ispublic(Grammar, :orchestrator_method)
    @test !isdefined(LineCableModels, :orchestrator_root)
    @test !isdefined(LineCableModels, :orchestrator_method)

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
            specializations = Method[method
                                     for method in methods(action)
                                     if method !== declared_method &&
                Base.morespecific(method.sig, declared_method.sig)]
            isempty(specializations) || push!(
                findings,
                "more-specific public methods: $(join(specializations, "; "))"
            )

            hooks = if required_hooks === nothing
                RequiredInterfaces.isInterface(root) ?
                Tuple(RequiredInterfaces.functions(
                    RequiredInterfaces.getInterface(root),
                )) : ()
            else
                required_hooks
            end
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

        @testset "$(nameof(action)) with $root" begin
            if isempty(findings)
                @test isempty(findings)
            else
                for finding in findings
                    @warn "Advisory orchestrator finding" finding
                end
                @test_skip isempty(findings)
            end
        end
        return findings
    end

    plot_implementors = (
        LineCableModels.Engine.LineParameterPlotDefinition,
        LineCableModels.Engine.LineParametersBenchmarkPlotDefinition,
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        LineCableModels.DataModel.CableCollectionPreviewPlotDefinition,
        LineCableModels.DataModel.SystemPreviewPlotDefinition,
        LineCableModels.DataModel.MaterialScalePlotDefinition
    )
    plot_findings = test_orchestrator(
        PlotBuilder.make_render,
        PlotBuilder.AbstractPlotDefinition,
        plot_implementors;
        required_hooks = (PlotBuilder.entitle, PlotBuilder.resolve, PlotBuilder.fetch),
        optional_hooks = (PlotBuilder.parse, PlotBuilder.finish)
    )

    report_implementors = (
        ReportBuilder.TableReportDefinition,
        ReportBuilder.CableConstantsTableDefinition,
        ReportBuilder.LineParametersTableDefinition,
        ReportBuilder.BenchmarkTableDefinition,
        ReportBuilder.MonteCarloTableDefinition,
        ReportBuilder.XLSXReportDefinition
    )
    report_findings = test_orchestrator(
        ReportBuilder.report,
        ReportBuilder.AbstractReportDefinition,
        report_implementors;
        optional_hooks = (ReportBuilder.finish,)
    )

    projection_findings = test_orchestrator(
        ParametricBuilder.project,
        ParametricBuilder.AbstractProjectionDefinition,
        ();
        optional_hooks = (ParametricBuilder.finish,)
    )

    @test plot_findings isa Vector{String}
    @test report_findings isa Vector{String}
    @test projection_findings isa Vector{String}
end
