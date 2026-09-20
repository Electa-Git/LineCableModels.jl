@testitem "Quality / TextDisplay / ownership and side-effect boundaries" tags=[:quality] begin
    using DataFrames

    # Current show/summary methods are owned by the domain. Observable
    # publication and display side effects are exercised by the behavioral
    # ReportBuilder/TextDisplay tests, not inferred from source tokens.
    owner_modules=(
        LineCableModels.Units,
        LineCableModels.InputValidation,
        LineCableModels.Materials,
        LineCableModels.Earth,
        LineCableModels.DataModel,
        LineCableModels.Engine,
        LineCableModels.ParametricBuilder,
        LineCableModels.UQ,
        LineCableModels.ReportBuilder
    )
    display_modules=(owner_modules..., LineCableModels.TextDisplay)
    for binding in names(LineCableModels; all = false, imported = true)
        isdefined(LineCableModels, binding) || continue
        owned_type=getfield(LineCableModels, binding)
        owned_type isa Union{DataType, UnionAll} || continue
        Base.isabstracttype(owned_type) && continue
        parentmodule(Base.unwrap_unionall(owned_type)) in owner_modules || continue
        @test which(summary, (IO, owned_type)).module in display_modules
        @test which(show, (IO, owned_type)).module in display_modules
        @test which(show, (IO, MIME"text/plain", owned_type)).module in display_modules
    end
    for owner in owner_modules
        for binding in names(owner; all = false, imported = false)
            isdefined(owner, binding) || continue
            owned_type=getfield(owner, binding)
            owned_type isa Union{DataType, UnionAll} || continue
            Base.isabstracttype(owned_type) && continue
            parentmodule(Base.unwrap_unionall(owned_type)) in owner_modules || continue
            @test which(summary, (IO, owned_type)).module in display_modules
            @test which(show, (IO, owned_type)).module in display_modules
            @test which(show, (IO, MIME"text/plain", owned_type)).module in display_modules
        end
    end

    report_methods=filter(method -> method.module === LineCableModels.ReportBuilder,
        collect(methods(DataFrame)))
    @test length(report_methods) == 1
    signature=Base.unwrap_unionall(only(report_methods).sig)
    @test signature.parameters[2] <: LineCableModels.ObservedResult

end
