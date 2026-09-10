@testitem "Gauntlet / spectral construction choices survive grids and retained reports" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    model=load_case(:two_insulated_wires;variation=ExactOverrides(frequencies=[1e3]))
    physical=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    controls=[(method,options=(samples,)) for method in (:trapz,:cim) for samples in (nothing,100_000)]
    selections=[formula(:default;options=(integration=choice,)) for choice in controls]
    formulations=Formulation(earth_impedance=Grid(selections),earth_admittance=Grid(selections);
        combine=:zip,options=physical)
    definition=benchmark_definition(model;id=:spectral_choices,source_file=@__FILE__,
        reference=Formulation(;options=physical),formulations)
    mktempdir() do directory
        value=run_benchmark(definition;directory)
        @test length(value.candidate_result)==length(controls)
        saved=read_benchmark(directory)
        artifact=report(BenchmarkTableDefinition(),saved)
        rows=filter(row->row.role===:candidate,artifact.table.formulations)
        @test allunique(rows.label)
        @test Set(artifact.table.comparisons.formulation_index)==Set(eachindex(controls))
        for (index,choice) in enumerate(controls)
            current=value.candidate_result[index]
            retained=saved.candidate.result[index]
            @test Z(retained)==Z(current)
            @test Y(retained)==Y(current)
            for slot in (:earth_impedance,:earth_admittance)
                @test getproperty(rows.record[index].requested,slot).options.integration==choice
                @test getproperty(details(retained).formulations.requested,slot).options.integration==choice
                for interaction in getproperty(details(retained).formulations.numerical,slot)
                    @test interaction.options.integration.method===Val(choice.method)
                    @test interaction.options.integration.options.samples===choice.options.samples
                end
            end
            for quantity in (Z,Y), component in (real,imag)
                @test component.(quantity(current))≈component.(quantity(value.reference_result)) rtol=1e-5 atol=1e-10
            end
        end
        @test run_benchmark(definition;directory).timings.execution.candidate.reused
    end
    for method in (:trapz,:cim)
        limited=formula(:default;options=(integration=(method,options=(samples=16,)),))
        @test_throws r"construction sample budget exhausted" compute(model.problem,
            Formulation(earth_impedance=limited,earth_admittance=limited;options=physical))
    end
end
