@testitem "QuantityUnits / units / labels and scaling" tags = [:unit] begin
    const QU = LineCableModels.QuantityUnits

    resistance_per_km = QU.units(:base, :ohm; per = (:kilo, :meter))
    resistance_per_m = QU.units(:base, :ohm; per = (:base, :meter))
    @test QU.get_label(resistance_per_km) == "Ω/km"
    @test QU.get_exp(resistance_per_km) == -3
    @test QU.scale_factor(resistance_per_m, resistance_per_km) == 1_000.0
    @test QU.scale_factor(resistance_per_km, resistance_per_m) == 0.001

    composite = QU.Units(
        base = [QU.Unit(name = :ohm), QU.Unit(name = :meter)],
        per = [QU.Unit(name = :second), QU.Unit(name = :meter)]
    )
    @test QU.get_label(composite) == "Ω.m/(second.m)"
    @test QU.get_label(QU.Units()) == ""
    @test_throws ArgumentError QU.get_exp(QU.units(:unsupported, :ohm))
end

@testitem "QuantityUnits / quantities / semantic selector contracts" tags = [:unit] begin
    const QU = LineCableModels.QuantityUnits

    @test QU.get_symbol(QU.quantity(R)) == "R"
    @test QU.get_symbol(QU.quantity(L)) == "L"
    @test QU.get_label(QU.quantity(C)) == "Shunt capacitance"
    @test QU.get_symbol(QU.quantity(Z, :angle)) == "Z"
    @test QU.get_label(QU.line_component_unit(:resistance, :per_length).units) == "Ω/km"
    @test QU.line_component_unit(:inductance, :per_length).scale == 1e6
    @test QU.line_component_unit(:capacitance, :total).scale == 1e6

    custom = QU.line_component_unit(
        :capacitance,
        :per_length;
        length_unit = :base,
        quantity_units = (capacitance = :nano,)
    )
    @test QU.get_label(custom.units) == "nF/m"
    @test custom.scale == 1e9

    @test_throws ArgumentError QU.quantity(R, :abs)
    @test_throws ArgumentError QU.quantity(identity)
    @test_throws ArgumentError QU.line_component_quantity(:unknown)
    @test_throws ArgumentError QU.line_component_unit(:resistance, :invalid)
end

@testitem "QuantityUnits / quantities / complete display policy" tags = [:unit] begin
    const QU = LineCableModels.QuantityUnits

    @test QU.get_label(QU.Unit(name = :custom, prefix = :milli)) == "mcustom"
    @test_throws ArgumentError QU.get_label(QU.Unit(name = :ohm, prefix = :unsupported))
    @test QU.get_label(QU.Units(
        base = QU.Unit[],
        per = [QU.Unit(name = :meter), QU.Unit(name = :hertz)]
    )) == "1/(m.Hz)"

    tag = QU.QuantityTag(Val(:custom))
    @test tag isa QU.QuantityTag{:custom}
    @test QU.QuantityTag(typeof(tag)) isa QU.QuantityTag{:custom}
    @test QU.get_label(tag) == "custom"
    @test QU.get_symbol(tag) == "custom"
    @test QU.get_label(QU.default_unit(tag)) == ""
    @test QU.get_label(QU.display_unit(tag)) == QU.get_label(QU.default_unit(tag))

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
        specification = QU.line_component_quantity(component)
        resolved = QU.line_component_unit(component, :per_length)
        @test QU.get_label(resolved.units) == unit_label
        @test QU.get_symbol(specification.tag) == symbol
        @test resolved.scale > 0
    end

    @test QU.quantity(Z) isa QU.QuantityTag{:series_impedance}
    @test QU.quantity(Y) isa QU.QuantityTag{:shunt_admittance}
    @test QU.quantity(:R) === QU.quantity(R)
    @test QU.quantity(:Z_abs) === QU.quantity(Z)
    @test QU.quantity(:Y_angle) === QU.quantity(:angle)
    @test QU.get_label(QU.quantity(Z, :abs)) == "Series impedance magnitude"
    @test QU.get_label(QU.quantity(Y, :abs)) == "Shunt admittance magnitude"
    @test QU.get_label(QU.quantity(Z, :angle)) == "Series impedance angle"
    @test QU.get_label(QU.quantity(Y, :angle)) == "Shunt admittance angle"
    @test_throws ArgumentError QU.quantity(Z, :unsupported)

    @test QU.get_label(QU.default_unit(Val(:frequency))) == "Hz"
    @test QU.get_label(QU.display_unit(:frequency)) == "Hz"
    @test QU.get_label(Val(:frequency)) == "Frequency"
    @test QU.get_label(:angle) == "Angle"
    @test QU.scale_factor(QU.display_unit(:resistance)) == 1_000.0
    @test QU.scale_factor(
        QU.QuantityTag{:resistance}(),
        :total,
        QU.units(:milli, :ohm)
    ) == 1_000.0
    @test QU.get_label(QU.display_unit(QU.QuantityTag{:resistance}(), :total)) == "Ω"
    @test QU.get_label(QU.display_unit(
        QU.QuantityTag{:resistance}(),
        :per_length;
        length_prefix = :mega
    )) == "Ω/Mm"

    by_component = QU.line_component_unit(
        :resistance, :per_length; quantity_units = Dict(:resistance =>
            :milli))
    by_semantic = QU.line_component_unit(
        :resistance,
        :per_length;
        quantity_units = Dict(:resistance => :micro)
    )
    @test QU.get_label(by_component.units) == "mΩ/km"
    @test QU.get_label(by_semantic.units) == "μΩ/km"
    @test QU.get_label(QU.line_component_unit(
        :resistance,
        :per_length;
        quantity_units = Dict(:unrelated => :nano)
    ).units) == "Ω/km"
    @test_throws ArgumentError QU.line_component_unit(:resistance, :per_length; quantity_units = 1)
    @test_throws ArgumentError QU.line_component_unit(
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
        quantity = QU.QuantityTag{name}()
        @test QU.get_label(quantity) == label
        @test QU.get_symbol(quantity) == symbol
        @test QU.get_label(QU.display_unit(quantity)) == unit_label
    end
    @test QU.default_unit(:resistance) isa QU.Units
    @test QU.scale_factor(
        QU.QuantityTag{:resistance}(),
        QU.display_unit(QU.QuantityTag{:resistance}())
    ) == 1000.0
    @test QU.get_label(QU.default_unit(QU.QuantityTag{(:series_impedance, :re)}())) == "Ω/m"
    @test QU.get_label(QU.default_unit(QU.QuantityTag{(:series_impedance, :im)}())) == "Ω/m"
    @test QU.get_label(QU.default_unit(QU.QuantityTag{(:shunt_admittance, :re)}())) == "S/m"
    @test QU.get_label(QU.default_unit(QU.QuantityTag{(:shunt_admittance, :im)}())) == "S/m"
end
