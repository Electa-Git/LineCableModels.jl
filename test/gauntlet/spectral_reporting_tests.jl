@testitem "Gauntlet / spectral construction choices survive grids and retained reports" tags=[:gauntlet_toolkit] setup=[GauntletSupport,TestFixtures] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1e3]))
    physical=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    controls=[(method,options=(;samples)) for method in (:trapz,:cim) for samples in (nothing,100_000)]
    selections=[formula(:default;options=(integration=choice,)) for choice in controls]
    formulations=Formulation(earth_impedance=Grid(selections),earth_admittance=Grid(selections);
        combine=:zip,options=physical)
    # Exhaustive labels and controls use current transport records. Their
    # numerical payload is distinguishable test data, not a spectral calculation.
    records=Gauntlet.formulation_record.(collect(formulations))
    source=TestFixtures.two_conductor_results(frequencies=[1e3])
    points=[LineParameters(Z(source).*index,Y(source),frequencies(source);
        details=ComputationDetails(merge(details(source).data,(formulations=record,)))) for (index,record) in enumerate(records)]
    transported=ParametricResult(nothing,points,(problems=[model.problem],formulations=records), ComputationDetails((;)))
    table=report(BenchmarkTableDefinition(),(reference=source,candidate=transported)).table
    @test allunique(table.formulations.label[table.formulations.role.===:candidate])
    @test Set(table.comparisons.candidate_point)==Set(eachindex(controls))
    for (index,choice) in enumerate(controls), slot in (:earth_impedance,:earth_admittance)
        @test getproperty(records[index].requested,slot).options.integration==choice
    end

    # Two adaptive constructions keep the actual compute/persist/report path.
    # Fixed 100,000-sample work does not belong to a metadata assertion.
    controls=filter(choice->choice.options.samples===nothing,controls)
    selections=[formula(:default;options=(integration=choice,)) for choice in controls]
    formulations=Formulation(earth_impedance=Grid(selections),earth_admittance=Grid(selections);
        combine=:zip,options=physical)
    definition=benchmark_definition(model;id=:spectral_choices,source_file=@__FILE__,
        reference=Formulation(;options=physical),formulations)
    mktempdir() do directory
        value=run_benchmark(definition;directory)
        @test length(value.candidate_result)==length(controls)
        saved=read_benchmark(directory)
        @test isconcretetype(eltype(saved.candidate.result))
        artifact=report(BenchmarkTableDefinition(),saved)
        rows=filter(row->row.role===:candidate,artifact.table.formulations)
        @test allunique(rows.label)
        @test Set(artifact.table.comparisons.candidate_point)==Set(eachindex(controls))
        for (index,choice) in enumerate(controls)
            current=value.candidate_result[index]
            retained=saved.candidate.result[index]
            @test Z(retained)==Z(current)
            @test Y(retained)==Y(current)
            for slot in (:earth_impedance,:earth_admittance)
                @test getproperty(details(retained).data.formulations.requested,slot).options.integration==choice
                for interaction in getproperty(details(current).data.formulations.numerical,slot)
                    @test interaction.options.integration.method===Val(choice.method)
                    @test interaction.options.integration.options.samples===choice.options.samples
                end
                for interaction in getproperty(details(retained).data.formulations.numerical,slot)
                    @test interaction.options.integration.method==(
                        type="Base.Val{$(repr(choice.method))}",fields=(;))
                    @test interaction.options.integration.options.samples===choice.options.samples
                end
            end
        end
        @test run_benchmark(definition;directory).timings.execution.candidate.reused
    end
end
