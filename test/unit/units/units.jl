@testitem "Units / immutable expressions and compatible conversion" tags = [:unit] begin
    const U = LineCableModels.Units

    per_meter = U.units(:base, :ohm; per = (:base, :meter))
    per_kilometer = U.units(:base, :ohm; per = (:kilo, :meter))
    @test per_meter.numerator isa Tuple
    @test per_meter.denominator isa Tuple
    @test U.label(per_kilometer) == "Ω/km"
    @test U.scale_factor(per_meter, per_kilometer) == 1_000.0
    @test U.scale_factor(per_kilometer, per_meter) == 0.001
    @test hash(per_meter) == hash(U.units(:base, :ohm; per = (:base, :meter)))

    composite = U.UnitExpr(
        (U.Unit(:ohm), U.Unit(:meter)),
        (U.Unit(:second), U.Unit(:meter))
    )
    @test U.label(composite) == "Ω.m/(s.m)"
    @test U.label(U.UnitExpr()) == ""
    @test U.label(inv(U.units(:base, :hertz))) == "1/Hz"

    @test_throws ArgumentError U.Unit(:custom)
    @test_throws ArgumentError U.Unit(:ohm, :unsupported)
    @test_throws ArgumentError U.Unit(:degree, :milli)
    @test_throws ArgumentError U.scale_factor(
        U.units(:base, :ohm),
        U.units(:base, :farad)
    )
end

@testitem "Units / quantity metadata and basis" tags = [:unit] begin
    const U = LineCableModels.Units

    expected = (
        (R, U.QuantityTag{:series_resistance}, "Series resistance", "R", "Ω/km"),
        (X, U.QuantityTag{:series_reactance}, "Series reactance", "X", "Ω/km"),
        (L, U.QuantityTag{:series_inductance}, "Series inductance", "L", "mH/km"),
        (G, U.QuantityTag{:shunt_conductance}, "Shunt conductance", "G", "S/km"),
        (B, U.QuantityTag{:shunt_susceptance}, "Shunt susceptance", "B", "S/km"),
        (C, U.QuantityTag{:shunt_capacitance}, "Shunt capacitance", "C", "μF/km"),
        (Z, U.QuantityTag{:series_impedance}, "Series impedance", "Z", "Ω/km"),
        (Y, U.QuantityTag{:shunt_admittance}, "Shunt admittance", "Y", "S/km"),
        (frequencies, U.QuantityTag{:frequency}, "Frequency", "f", "Hz")
    )

    for (selector, tag_type, name, notation, unit_text) in expected
        tag = U.quantity(selector)
        @test tag isa tag_type
        @test U.label(tag) == name
        @test U.symbol(tag) == notation
        @test U.label(U.display_unit(tag, :per_length)) == unit_text
        @test U.label(tag, U.display_unit(tag, :per_length)) == "$name [$unit_text]"
    end

    @test U.label(U.native_unit(U.quantity(R), :per_length)) == "Ω/m"
    @test U.label(U.display_unit(U.quantity(R), :total)) == "Ω"
    @test U.label(U.display_unit(U.quantity(C), :total)) == "μF"
    @test_throws ArgumentError U.display_unit(U.quantity(R), :invalid)
end

@testitem "Units / transforms and angular conversion" tags = [:unit] begin
    const U = LineCableModels.Units

    z_magnitude = U.quantity(Z, abs)
    z_angle = U.quantity(Z, angle)
    y_magnitude = U.quantity(Y, abs)
    y_angle = U.quantity(Y, angle)

    @test z_magnitude isa U.QuantityTag{(:series_impedance, :magnitude)}
    @test z_angle isa U.QuantityTag{(:series_impedance, :phase_angle)}
    @test y_magnitude isa U.QuantityTag{(:shunt_admittance, :magnitude)}
    @test y_angle isa U.QuantityTag{(:shunt_admittance, :phase_angle)}
    @test U.label(z_magnitude) == "Series impedance magnitude"
    @test U.symbol(y_angle) == "∠Y"
    @test U.label(U.native_unit(z_angle)) == "rad"
    @test U.label(U.display_unit(z_angle)) == "°"
    @test U.scale_factor(U.native_unit(z_angle), U.display_unit(z_angle)) ≈ 180 / π
    @test U.quantity(Z, :abs) === z_magnitude
    @test_throws ArgumentError U.quantity(Z, :unsupported)
    @test_throws MethodError U.quantity(identity)
end

@testitem "Units / locked public vocabulary" tags = [:unit] begin
    const U = LineCableModels.Units

    for name in (
            :Unit, :UnitExpr, :QuantityTag, :units, :quantity, :native_unit,
            :display_unit, :scale_factor, :label, :symbol
    )
        @test name in names(U)
    end
    @test !isdefined(LineCableModels, :QuantityUnits)
    @test !isdefined(U, :default_unit)
    @test !isdefined(U, :get_label)
    @test !isdefined(U, :get_symbol)
end
