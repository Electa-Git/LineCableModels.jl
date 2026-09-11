@testitem "Gauntlet / captured declarations survive source edits and collections archive" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels, JLD2, SHA, TOML
model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
reference=BenchmarkCalculation(:first,model.problem,Formulation(;options))
candidate=BenchmarkCalculation(:second,model.problem,Formulation(earth_impedance=:Pollaczek1926;options))
definition=benchmark_definition(:portable,model.id,:fixture,@__FILE__,model,
    reference,candidate,(; quantities=(:Z,:Y,:G)),(;))
mktempdir() do parent
    source=joinpath(parent,"execution.jl");write(source,"global declared_calculation_factor=1\n")
    original_source=read(source)
    campaign=joinpath(parent,"campaign")
    outcome=only(run_campaign(campaign,[definition];on_error=:fail,
        execution_sources=[(path=source,module_name=:Main)]))
    @test outcome.state===:complete
    @test read(source)==original_source
    state=TOML.parsefile(joinpath(campaign,"portable","state.toml"))
    attempt=joinpath(campaign,"portable",state["current"])
    @test isfile(joinpath(attempt,"reference","calculation.jld2"))
    @test isfile(joinpath(attempt,"candidate","calculation.jld2"))
    @test only(resume_campaign(campaign)).result.timings.execution.reference.reused
    @test read(source)==original_source
    bundle=lock_campaign(campaign,joinpath(parent,"bundle"))
    moved=joinpath(parent,"moved");mv(bundle.path,moved)
    @test only(read_campaign(moved)).reference.result.Z==outcome.result.reference_result.Z
    before=readdir(moved)
    @test only(campaign_status(moved)).state===:complete
    @test readdir(moved)==before
    write(source,"global declared_calculation_factor=2\n")
    @test only(resume_campaign(campaign)).result.timings.execution.reference.reused
    @test Base.invokelatest(getproperty,Main,:declared_calculation_factor) == 1
    @test read(source,String)=="global declared_calculation_factor=2\n"
    rm(source)
    @test only(resume_campaign(campaign)).result.timings.execution.reference.reused
    package=package_collection(:fixture,v"1.0.0";reason="Offline transport verification",
        bundles=[moved],output=joinpath(parent,"release"))
    @test bytes2hex(open(sha256,package.archive))==package.archive_sha256
    @test isfile(package.archive)
end
end
