@testitem "UnitHandler / units / labels and scaling" tags = [:unit] begin
    const UH = LineCableModels.UnitHandler

    resistance_per_km = UH.units(:base, :ohm; per = (:kilo, :meter))
    resistance_per_m = UH.units(:base, :ohm; per = (:base, :meter))
    @test UH.get_label(resistance_per_km) == "Ω/km"
    @test UH.get_exp(resistance_per_km) == -3
    @test UH.scale_factor(resistance_per_m, resistance_per_km) == 1_000.0
    @test UH.scale_factor(resistance_per_km, resistance_per_m) == 0.001

    composite = UH.Units(
        base = [UH.Unit(name = :ohm), UH.Unit(name = :meter)],
        per = [UH.Unit(name = :second), UH.Unit(name = :meter)]
    )
    @test UH.get_label(composite) == "Ω.m/(second.m)"
    @test UH.get_label(UH.Units()) == ""
    @test_throws ArgumentError UH.get_exp(UH.units(:unsupported, :ohm))
end

@testitem "UnitHandler / quantities / semantic selector contracts" tags = [:unit] begin
    const UH = LineCableModels.UnitHandler

    @test UH.get_symbol(UH.quantity(R)) == "R"
    @test UH.get_symbol(UH.quantity(L)) == "L"
    @test UH.get_label(UH.quantity(C)) == "Shunt capacitance"
    @test UH.get_symbol(UH.quantity(Z, :angle)) == "Z"
    @test UH.get_label(UH.line_component_unit(:resistance, :per_length).units) == "Ω/km"
    @test UH.line_component_unit(:inductance, :per_length).scale == 1e6
    @test UH.line_component_unit(:capacitance, :total).scale == 1e6

    custom = UH.line_component_unit(
        :capacitance,
        :per_length;
        length_unit = :base,
        quantity_units = (capacitance = :nano,)
    )
    @test UH.get_label(custom.units) == "nF/m"
    @test custom.scale == 1e9

    @test_throws ArgumentError UH.quantity(R, :abs)
    @test_throws ArgumentError UH.quantity(identity)
    @test_throws ArgumentError UH.line_component_quantity(:unknown)
    @test_throws ArgumentError UH.line_component_unit(:resistance, :invalid)
end

@testitem "UnitHandler / quantities / complete display policy" tags = [:unit] begin
    const UH = LineCableModels.UnitHandler

    @test UH.get_label(UH.Unit(name = :custom, prefix = :milli)) == "mcustom"
    @test_throws ArgumentError UH.get_label(UH.Unit(name = :ohm, prefix = :unsupported))
    @test UH.get_label(UH.Units(
        base = UH.Unit[],
        per = [UH.Unit(name = :meter), UH.Unit(name = :hertz)]
    )) == "1/(m.Hz)"

    tag = UH.QuantityTag(Val(:custom))
    @test tag isa UH.QuantityTag{:custom}
    @test UH.QuantityTag(typeof(tag)) isa UH.QuantityTag{:custom}
    @test UH.get_label(tag) == "custom"
    @test UH.get_symbol(tag) == "custom"
    @test UH.get_label(UH.default_unit(tag)) == ""
    @test UH.get_label(UH.display_unit(tag)) == UH.get_label(UH.default_unit(tag))

    component_expectations = Dict(
        :frequency => ("Hz", "f"),
        :resistance => ("Ω/km", "R"), :reactance => ("Ω/km", "X"),
        :inductance => ("mH/km", "L"), :conductance => ("S/km", "G"),
        :susceptance => ("S/km", "B"), :capacitance => ("μF/km", "C"),
        :series_impedance => ("Ω/km", "Z"),
        :shunt_admittance => ("S/km", "Y"),
        :angle => ("°", "∠"),
        :series_impedance_absolute_error => ("Ω/km", "ΔZ"),
        :series_impedance_relative_error => ("", "εZ"),
        :shunt_admittance_absolute_error => ("S/km", "ΔY"),
        :shunt_admittance_relative_error => ("", "εY")
    )
    for (component, (unit_label, symbol)) in component_expectations
        specification = UH.line_component_quantity(component)
        resolved = UH.line_component_unit(component, :per_length)
        @test UH.get_label(resolved.units) == unit_label
        @test UH.get_symbol(specification.tag) == symbol
        @test resolved.scale > 0
    end

    @test UH.quantity(Z) isa UH.QuantityTag{:series_impedance}
    @test UH.quantity(Y) isa UH.QuantityTag{:shunt_admittance}
    @test UH.get_label(UH.quantity(Z, :abs)) == "Series impedance magnitude"
    @test UH.get_label(UH.quantity(Y, :abs)) == "Shunt admittance magnitude"
    @test UH.get_label(UH.quantity(Z, :angle)) == "Series impedance angle"
    @test UH.get_label(UH.quantity(Y, :angle)) == "Shunt admittance angle"
    @test_throws ArgumentError UH.quantity(Z, :unsupported)

    @test UH.get_label(UH.default_unit(Val(:frequency))) == "Hz"
    @test UH.get_label(UH.display_unit(:frequency)) == "Hz"
    @test UH.get_label(Val(:frequency)) == "Frequency"
    @test UH.get_label(:angle) == "Angle"
    @test UH.scale_factor(UH.display_unit(:resistance)) == 1_000.0
    @test UH.scale_factor(
        UH.QuantityTag{:resistance}(),
        :total,
        UH.units(:milli, :ohm)
    ) == 1_000.0
    @test UH.get_label(UH.display_unit(UH.QuantityTag{:resistance}(), :total)) == "Ω"
    @test UH.get_label(UH.display_unit(
        UH.QuantityTag{:resistance}(),
        :per_length;
        length_prefix = :mega
    )) == "Ω/Mm"

    by_component = UH.line_component_unit(
        :resistance, :per_length; quantity_units = Dict(:resistance =>
            :milli))
    by_semantic = UH.line_component_unit(
        :resistance,
        :per_length;
        quantity_units = Dict(:resistance => :micro)
    )
    @test UH.get_label(by_component.units) == "mΩ/km"
    @test UH.get_label(by_semantic.units) == "μΩ/km"
    @test UH.get_label(UH.line_component_unit(
        :resistance,
        :per_length;
        quantity_units = Dict(:unrelated => :nano)
    ).units) == "Ω/km"
    @test_throws ArgumentError UH.line_component_unit(:resistance, :per_length; quantity_units = 1)
    @test_throws ArgumentError UH.line_component_unit(
        :resistance,
        :per_length;
        quantity_units = Dict(:resistance => "milli")
    )

    quantity_expectations = Dict(
        :frequency => ("Frequency", "f", "Hz"),
        :inductance => ("Series inductance", "L", "mH/km"),
        :capacitance => ("Shunt capacitance", "C", "μF/km"),
        :conductance => ("Shunt conductance", "G", "S/km"),
        :series_impedance => ("Series impedance", "Z", "Ω/km"),
        :shunt_admittance => ("Shunt admittance", "Y", "S/km"),
        :reactance => ("Inductive reactance", "X", "Ω/km"),
        :susceptance => ("Capacitive susceptance", "B", "S/km"),
        :angle => ("Angle", "∠", "°")
    )
    for (name, (label, symbol, unit_label)) in quantity_expectations
        quantity = UH.QuantityTag{name}()
        @test UH.get_label(quantity) == label
        @test UH.get_symbol(quantity) == symbol
        @test UH.get_label(UH.display_unit(quantity)) == unit_label
    end
    @test UH.default_unit(:resistance) isa UH.Units
    @test UH.scale_factor(
        UH.QuantityTag{:resistance}(),
        UH.display_unit(UH.QuantityTag{:resistance}())
    ) == 1000.0
    @test UH.get_label(UH.default_unit(UH.QuantityTag{(:series_impedance, :re)}())) == "Ω/m"
    @test UH.get_label(UH.default_unit(UH.QuantityTag{(:series_impedance, :im)}())) == "Ω/m"
    @test UH.get_label(UH.default_unit(UH.QuantityTag{(:shunt_admittance, :re)}())) == "S/m"
    @test UH.get_label(UH.default_unit(UH.QuantityTag{(:shunt_admittance, :im)}())) == "S/m"
end
