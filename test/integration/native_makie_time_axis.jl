@testitem "Makie addons / retained transient time axes and frequency scales" tags=[:visual, :makie, :slow] setup=[
    UseNativePlotSupport, TestFixtures
] begin
    get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false")=="true" ||
        error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual test")
    using CairoMakie

    times=[0.0,1e-6,2e-6]
    voltages=reshape(collect(1.0:12.0),2,2,3)
    response=TransientResponse(times,voltages,-voltages,
        ["cable:1:core","cable:2:core"],ComputationDetails())
    observed=ObservedResult(response,(voltage,))
    for (options,multiplier,label) in (((;),1.0,"Time [s]"),
            ((time_unit=:micro,),1e6,"Time [μs]"))
        page=LineCableModels.plot(observed;options...,
            backend=:cairo,display_plot=false,open_export=false)
        @test !isempty(Makie.colorbuffer(page.figure))
        @test length(page.axes)==4
        @test all(axis -> axis.xscale[]===identity,page.axes)
        @test all(axis -> first(split(String(axis.xlabel[]),"  × "))==label,page.axes)
        curves=[series for binding in page.plot_state.axis_bindings for series in binding.series]
        @test length(curves)==4
        for series in curves
            @test series.xdata≈times.*multiplier
        end
    end
    @test observed.quantities[1].coordinates.times==times

    parameters=TestFixtures.two_conductor_results(;frequencies=[50.0,100.0,500.0])
    frequency_page=LineCableModels.plot(ObservedResult(parameters,(R,X)),(R,);
        backend=:cairo,display_plot=false,open_export=false)
    @test !isempty(Makie.colorbuffer(frequency_page.figure))
    @test all(axis -> axis.xscale[]===log10,frequency_page.axes)
end
