@testitem "Gauntlet / mutable drafts, explicit acceptance and immutable release packages" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels, SHA, TOML
    using .GauntletSupport.Gauntlet
    using LineCableModels.ReportBuilder
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1.,37.]))
    options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    formulation=Formulation(;options)
    definition=benchmark_definition(model;id=:accepted,source_file=@__FILE__,reference=formulation,
        formulations=Formulation(earth_impedance=Grid((:default,:Pollaczek1926));options),collection=:fixture)
    sibling=benchmark_definition(model;id=:sibling,source_file=@__FILE__,reference=formulation,
        formulations=formulation,collection=:fixture)
    many_reference=benchmark_definition(model;id=:many_reference,source_file=@__FILE__,
        reference=Formulation(earth_impedance=Grid((:default,:Pollaczek1926));options),formulations=formulation)
    @test_throws r"before execution" run_benchmark(many_reference)
    mktempdir() do root
        campaign=joinpath(root,"staging")
        initial=run_campaign(campaign,[definition,sibling];on_error=:fail)
        sibling_state=read(joinpath(campaign,"sibling","state.toml"))
        rerun=only(run_campaign(campaign,[definition];on_error=:fail))
        @test read(joinpath(campaign,"sibling","state.toml"))==sibling_state
        @test !rerun.result.timings.execution.reference.reused
        retained_draft=read_benchmark(joinpath(campaign,"accepted"))
        restored_labels=report(BenchmarkTableDefinition(),retained_draft).table.formulations.label
        @test restored_labels==rerun.result.report.table.formulations.label
        @test retained_draft.reference.metadata.repository.commit == repository_revision().commit
        @test retained_draft.reference.metadata.active_project == Base.active_project()
        @test any(entry -> endswith(entry.path,"Project.toml"),retained_draft.reference.metadata.implementation)
        @test_throws r"inspected" lock_campaign(campaign,joinpath(root,"wrong");benchmarks=:accepted,expected="wrong")
        @test !ispath(joinpath(root,"wrong"))
        manifest=TOML.parsefile(joinpath(campaign,"campaign.toml"))
        push!(manifest["benchmarks"],"pending")
        open(io -> TOML.print(io,manifest),joinpath(campaign,"campaign.toml"),"w")
        @test_throws r"no completed attempt" lock_campaign(campaign,joinpath(root,"whole"))
        @test !ispath(joinpath(root,"whole"))
        bundle=lock_campaign(campaign,joinpath(root,"vault");benchmarks=:accepted,expected=rerun.identity)
        @test_throws r"already a locked bundle" lock_campaign(bundle.path,joinpath(root,"relock"))
        figure=joinpath(root,"selected.svg")
        Base.write(figure,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"10\" height=\"10\"/>")
        @test_throws r"selection" lock_campaign(campaign,joinpath(root,"invalid_figure");benchmarks=:accepted,
            illustrations=[(path=figure,benchmark=:accepted,caption="Fixture")])
        illustrated=lock_campaign(campaign,joinpath(root,"illustrated");benchmarks=:accepted,
            illustrations=[(path=figure,benchmark=:accepted,caption="Fixture",
                selection=(problem=1,quantities=["R"],formulations=[1,2]))])
        document=TOML.parsefile(joinpath(illustrated.path,"bundle.toml"))
        @test only(document["illustrations"])["selection"]["formulations"]==[1,2]
        @test read(joinpath(illustrated.path,only(document["illustrations"])["path"]))==read(figure)
        @test length(read_campaign(bundle.path))==1
        @test only(read_campaign(bundle.path)).id==:accepted
        @test_throws r"immutable" cleanup_work(;work_root=bundle.path)
        @test_throws r"immutable" cleanup_work(;work_root=root)
        @test length(read_campaign(bundle.path))==1
        guarded=joinpath(root,"old_staging")
        mkpath(guarded)
        cp(bundle.path,joinpath(guarded,"staging"))
        @test_throws r"immutable" prepare_staging(;artifact_root=guarded,force=true)
        @test length(read_campaign(joinpath(guarded,"staging")))==1

        @test_throws r"immutable" run_benchmark(definition;directory=joinpath(bundle.path,"another"))
        @test_throws r"immutable" lock_campaign(campaign,bundle.path)
        @test_throws r"locked bundles" package_collection(:fixture,v"1.0.0";
            bundles=[campaign],reason="draft rejected",output=joinpath(root,"bad_release"))
        package=package_collection(:fixture,v"1.0.0";bundles=[bundle.path],reason="Accepted comparisons",output=joinpath(root,"release"))
        archive=read(package.archive)
        repeated=package_collection(:fixture,v"1.0.0";bundles=[bundle.path],reason="Accepted comparisons",output=package.path)
        @test repeated.tree_hash==package.tree_hash
        @test read(package.archive)==archive
        @test_throws r"new version" package_collection(:fixture,v"1.0.0";
            bundles=[bundle.path],reason="different definition",output=package.path)
        bindings=joinpath(root,"Artifacts.toml")
        publication=bind_published_artifact(package.path,"file://"*package.archive;artifacts_toml=bindings)
        @test publication.tree_hash==package.tree_hash
        @test Set(keys(TOML.parsefile(bindings)))==Set(("gauntlet_fixture_v1_0_0",))
        before=read(bindings)
        @test_throws r"immutable" bind_published_artifact(package.path,"file://"*package.archive;
            artifacts_toml=joinpath(bundle.path,"Artifacts.toml"))
        wrong=joinpath(root,"wrong.tar.gz");Base.write(wrong,"wrong")
        @test_throws r"checksum" bind_published_artifact(package.path,"file://"*wrong;artifacts_toml=bindings)
        @test read(bindings)==before
        rm(campaign;recursive=true)
        restored=only(read_campaign(bundle.path))
        @test length(report(BenchmarkTableDefinition(),restored).published.comparisons)==60
        extra=joinpath(bundle.path,"unexpected");Base.write(extra,"unexpected")
        @test_throws r"unexpected" read_campaign(bundle.path)
        rm(extra)
        @test length(read_campaign(bundle.path))==1
        document=TOML.parsefile(joinpath(bundle.path,"bundle.toml"));document["note"]="changed"
        open(io->TOML.print(io,document),joinpath(bundle.path,"bundle.toml"),"w")
        @test_throws r"inventory changed" read_campaign(bundle.path)
    end
end
