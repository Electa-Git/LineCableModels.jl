@testitem "TransientAnalysis / response arrays" tags=[:unit, :engine] begin
    import LineCableModels.Units: units

    t = [0.0, 0.1, 0.2]
    voltages = reshape(collect(1.0:12.0), 2, 2, 3)
    currents = -voltages
    coordinates = ["cable:1:core", "cable:2:core"]
    retained = ComputationDetails(; formula=:nlt, window=(tmax=0.2, step=0.1))
    response = TransientResponse(t, voltages, currents, coordinates, retained)
    @test response.t === t
    @test response.voltage === voltages
    @test response.current === currents
    @test response.coordinates === coordinates
    @test details(response) === retained
    @test computation_details(TransientFormulation, response) === retained
    @test observe(response, voltage, 1, 2, :) == voltages[1, 2, :]
    @test observe(response, current, "receiving", "cable:1:core", :) == currents[2, 1, :]
    time = LineCableModels.Units.Quantity{:time}()
    @test label(native_unit(time)) == "s"
    @test label(display_unit(time)) == "s"
    for (identifier, name, unit, notation) in
            ((:terminal_voltage, "Terminal voltage", :volt, "V"),
             (:terminal_current, "Terminal current", :ampere, "I"))
        quantity = LineCableModels.Units.Quantity{identifier}()
        @test native_unit(quantity) == units(:base, unit)
        @test display_unit(quantity) == units(:base, unit)
        @test label(quantity) == name
        @test LineCableModels.Units.symbol(quantity) == notation
        @test label(units(:milli, unit)) == "m$(label(units(:base, unit)))"
        @test label(units(:kilo, unit)) == "k$(label(units(:base, unit)))"
    end
end
