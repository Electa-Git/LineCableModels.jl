@testitem "Quality / report orchestration / direct contracts" tags = [:quality] begin
    using RequiredInterfaces

    const Grammar = LineCableModels.Grammar
    const ReportBuilder = LineCableModels.ReportBuilder

    @test !isdefined(Grammar, Symbol("@orchestrator"))
    @test !isdefined(Grammar, :orchestrator_root)
    @test !isdefined(Grammar, :orchestrator_method)
    @test !isdefined(ReportBuilder, :entitle)
    @test !isdefined(ReportBuilder, :finish)

    root = ReportBuilder.AbstractReportDefinition
    @test RequiredInterfaces.isInterface(root)
    hooks = Set(RequiredInterfaces.functions(RequiredInterfaces.getInterface(root)))
    @test hooks == Set((ReportBuilder.select, ReportBuilder.tabulate))

    report_methods = collect(methods(ReportBuilder.report))
    @test length(report_methods) == 1
    report_method = only(report_methods)
    @test report_method.module === ReportBuilder
    @test report_method.sig == Tuple{
        typeof(ReportBuilder.report),
        ReportBuilder.AbstractReportDefinition,
        Any
    }

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

    implementors = Tuple(RequiredInterfaces.nonabstract_subtypes(root))
    missing = Tuple(
        (implementor, hook) for implementor in implementors
    for hook in hooks
    if !any(method -> signature_mentions(method.sig, implementor), methods(hook))
    )
    @test isempty(missing)

    publication_methods = Method[method
                                 for method in methods(LineCableModels.observables)
                                 if method.nargs == 3]
    @test length(publication_methods) == 1
    publication_method = only(publication_methods)
    @test publication_method.module === Grammar
    @test publication_method.sig ==
          Tuple{typeof(LineCableModels.observables), Any, Tuple}
end
