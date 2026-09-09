@testitem "Gauntlet / one selected comparison and rejection before calculation" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels, JLD2
    using LineCableModels.Engine: compare
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using .GauntletSupport.Gauntlet
    using LineCableModels.Grammar: AbstractCoreResult
    const executions=Bool[]
    const comparisons=Tuple[]
    struct AlternativeResult <: AbstractCoreResult
        matrices::Tuple{Array{ComplexF64,3},Array{ComplexF64,3}}
        samples::Vector{Float64}
    end
    LineCableModels.observe(value::AlternativeResult, ::typeof(Z))=value.matrices[1]
    LineCableModels.observe(value::AlternativeResult, ::typeof(Y))=value.matrices[2]
    LineCableModels.observe(value::AlternativeResult, ::typeof(frequencies))=value.samples
    LineCableModels.frequencies(value::AlternativeResult)=value.samples
    LineCableModels.basis(::AlternativeResult)=:pul
    LineCableModels.domain(::AlternativeResult)=PhaseDomain
    LineCableModels.details(::AlternativeResult)=(;)
    function LineCableModels.Engine.compare(a::AlternativeResult, b::AlternativeResult, quantity::typeof(Z); kwargs...)
        push!(comparisons, (quantity, kwargs[:band], kwargs[:normalization]))
        return invoke(compare, Tuple{AbstractCoreResult,AbstractCoreResult,typeof(quantity)}, a, b, quantity; kwargs...)
    end
    struct AuditFormulation <: LineCableModels.Grammar.AbstractFormulation
        candidate::Bool
    end
    function LineCableModels.compute(problem::LineParametersProblem, formula::AuditFormulation; options=(;))
        push!(executions, formula.candidate)
        z=ones(ComplexF64,2,2,length(problem.frequencies))
        y=fill(1im,2,2,length(problem.frequencies))
        if formula.candidate
            z[:,:,1].=100
            z[:,:,2:end].=2
        end
        return AlternativeResult((z,ComplexF64.(y)),copy(problem.frequencies))
    end
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,10.,100.]))
    reference=BenchmarkCalculation(:reference,model.problem,AuditFormulation(false))
    candidate=BenchmarkCalculation(:candidate,model.problem,AuditFormulation(true))
    settings=(quantities=(:Z,),bands=((10.,100.),))
    definition=benchmark_definition(:selected,model.id,:fixture,@__FILE__,model,reference,candidate,settings,(;))
    mktempdir() do parent
        directory=joinpath(parent,"selected")
        result=run_benchmark(definition;directory)
        @test executions == [false,true]
        @test comparisons == [(Z,(10.,100.),:reference_rms)]
        @test length(result.comparison)==1
        @test only(result.comparison).quantity===:Z
        @test all(==(1.),only(result.comparison).error.relative)
        @test !hasproperty(result,:configured_comparisons)
        saved=read_benchmark(result.artifact;load_results=true)
        record=only(saved.analyses)
        @test record["comparison_settings"] == definition.comparison_settings
        @test isequal(only(record["reference_comparison"]).relative,only(result.comparison).error.relative)
        tables=report(BenchmarkTableDefinition(false),saved).table
        @test all(==(100.),tables.terms.relative_rms_percent)
        @test all(==(:Z),tables.terms.quantity)
        # A different result representation uses the same report sequence and public observations.
        alternative=merge(saved,(reference=merge(saved.reference,(result=result.reference,)),
            candidate=merge(saved.candidate,(result=result.candidate,))))
        @test isequal(report(BenchmarkTableDefinition(false),alternative).table.terms,tables.terms)
        @test length(comparisons)==1
        @test length(executions)==2
        observe(result.candidate,Z)[1,1,2]+=1
        @test_throws r"modified after loading" report(BenchmarkTableDefinition(false),alternative)
    end
    for invalid in ((bands=(:invalid_band,),),(bands=((100.,10.),),),
            (normalizations=(:invalid,),),(fundamental=-1.,),(harmonics=0,),
            (atol=-1.,),(atol=(G=-1.,),),(unsupported=(G="",),),
            (statistics=(:unknown,),))
        empty!(executions);empty!(comparisons)
        declaration=benchmark_definition(:invalid,model.id,:fixture,@__FILE__,model,
            reference,candidate,invalid,(;))
        mktempdir() do parent
            directory=joinpath(parent,"invalid")
            @test_throws ArgumentError run_benchmark(declaration;directory)
            @test !ispath(directory)
            @test isempty(executions)
            @test isempty(comparisons)
            @test_throws ArgumentError run_campaign(directory,[declaration])
            @test !ispath(directory)
        end
    end
    empty!(executions)
    limits=(Z=(absolute=1.,relative=1.),Y=(absolute=1.,relative=1.))
    unsupported_limits=benchmark_definition(:unsupported_limits,model.id,:fixture,
        @__FILE__,model,reference,candidate,(;),(reference=limits,))
    @test_throws r"only to declared moment comparisons" run_benchmark(unsupported_limits)
    @test isempty(executions)
    @test !isdefined(GauntletSupport.Gauntlet,:LineParametersPolicy)
    @test !isdefined(GauntletSupport.Gauntlet,:UQMomentPolicy)
    @test !isdefined(GauntletSupport.Gauntlet,:benchmark_calculation)
end
