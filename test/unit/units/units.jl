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
        (R, U.Quantity{:series_resistance}, "Series resistance", "R", "Ω/km"),
        (X, U.Quantity{:series_reactance}, "Series reactance", "X", "Ω/km"),
        (L, U.Quantity{:series_inductance}, "Series inductance", "L", "mH/km"),
        (G, U.Quantity{:shunt_conductance}, "Shunt conductance", "G", "S/km"),
        (B, U.Quantity{:shunt_susceptance}, "Shunt susceptance", "B", "S/km"),
        (C, U.Quantity{:shunt_capacitance}, "Shunt capacitance", "C", "μF/km"),
        (Z, U.Quantity{:series_impedance}, "Series impedance", "Z", "Ω/km"),
        (Y, U.Quantity{:shunt_admittance}, "Shunt admittance", "Y", "S/km"),
        (frequencies, U.Quantity{:frequency}, "Frequency", "f", "Hz")
    )

    for (selector, tag_type, name, notation, unit_text) in expected
        tag = U.quantity(selector)
        @test tag isa tag_type
        @test U.label(tag) == name
        @test U.symbol(tag) == notation
        @test U.label(U.display_unit(tag, :pul)) == unit_text
        @test U.label(tag, U.display_unit(tag, :pul)) == "$name [$unit_text]"
    end

    @test U.label(U.native_unit(U.quantity(R), :pul)) == "Ω/m"
    @test U.label(U.display_unit(U.quantity(R), :total)) == "Ω"
    @test U.label(U.display_unit(U.quantity(C), :total)) == "μF"
    @test U.label(U.display_unit(U.quantity(R), :pul; prefix = :milli)) ==
          "mΩ/km"
    @test U.label(U.display_unit(U.quantity(C), :total; prefix = :nano)) == "nF"
    @test U.display_unit(U.quantity(R), :pul, nothing) ==
          U.display_unit(U.quantity(R), :pul)
    @test U.label(U.display_unit(U.quantity(R), :pul, :micro)) == "μΩ/km"
    explicit=U.units(:nano, :ohm; per = (:base, :meter))
    @test U.display_unit(U.quantity(R), :pul, explicit) === explicit
    @test_throws ArgumentError U.display_unit(U.quantity(R), :pul, 1)
    @test_throws ArgumentError U.display_unit(U.quantity(R), :invalid)

    @test quantity(R) isa U.Quantity{:series_resistance}
    @test label(R) == "Series resistance"
    @test symbol(R) == "R"
    @test label(native_unit(R, :pul)) == "Ω/m"
    @test label(display_unit(R, :pul)) == "Ω/km"
    @test label(display_unit(R, :total)) == "Ω"
    @test scale_factor(R, display_unit(R)) == 1_000.0
    @test scale_factor(R, :pul, display_unit(R, :pul)) == 1_000.0
end

@testitem "Grammar / observable unit targets / aligned normalization" tags=[:unit] begin
    const Grammar=LineCableModels.Grammar
    const U=LineCableModels.Units

    requests=(resistance = R, phase = (Z, angle), capacitance = C)
    named=Grammar.unit_targets(
        requests,
        :pul;
        length_prefix = :base,
        overrides = (
            resistance = :milli,
            phase = U.units(:base, :radian)
        )
    )
    tupled=Grammar.unit_targets(
        values(requests),
        :pul;
        length_prefix = :base,
        overrides = Dict(
            R => :milli,
            (Z, angle) => U.units(:base, :radian)
        )
    )
    @test values(named) == tupled
    @test keys(named) == keys(requests)
    @test all(target -> target isa U.UnitExpr, values(named))
    @test named.resistance == U.units(:milli, :ohm; per = (:base, :meter))
    @test named.phase == U.units(:base, :radian)
    @test named.capacitance == U.units(:micro, :farad; per = (:base, :meter))
    @test Base.ispublic(Grammar, :unit_targets)
    @test !isdefined(LineCableModels, :unit_targets)
    @test_throws ArgumentError Grammar.unit_targets((R,), :pul; overrides = 1)
end

@testitem "Units / transforms and angular conversion" tags = [:unit] begin
    const U = LineCableModels.Units

    z_magnitude = U.quantity(Z, abs)
    z_angle = U.quantity(Z, angle)
    y_magnitude = U.quantity(Y, abs)
    y_angle = U.quantity(Y, angle)

    @test z_magnitude isa U.Quantity{(:series_impedance, :magnitude)}
    @test z_angle isa U.Quantity{(:series_impedance, :phase_angle)}
    @test y_magnitude isa U.Quantity{(:shunt_admittance, :magnitude)}
    @test y_angle isa U.Quantity{(:shunt_admittance, :phase_angle)}
    @test U.label(z_magnitude) == "Series impedance magnitude"
    @test U.symbol(y_angle) == "∠Y"
    @test U.label(U.native_unit(z_angle)) == "rad"
    @test U.label(U.display_unit(z_angle)) == "°"
    @test U.scale_factor(U.native_unit(z_angle), U.display_unit(z_angle)) ≈ 180 / π
    @test quantity(Z, abs) isa U.Quantity{(:series_impedance, :magnitude)}
    @test label(Z, abs) == "Series impedance magnitude"
    @test symbol(Z, abs) == "|Z|"
    @test symbol(Z, angle) == "∠Z"
    @test label(native_unit(Z, angle)) == "rad"
    @test label(display_unit(Z, angle)) == "°"
    @test label(native_unit(Z, angle, :pul)) == "rad"
    @test label(display_unit(Z, angle, :pul)) == "°"
    @test scale_factor(Z, angle, display_unit(Z, angle)) ≈ 180 / π
    @test scale_factor(Z, angle, :pul, display_unit(Z, angle, :pul)) ≈ 180 / π
    @test_throws MethodError U.quantity(Z, :abs)
    @test_throws MethodError U.quantity(identity)
end

@testitem "Units / selector metadata / extension locality" tags = [:unit] begin
    const U = LineCableModels.Units

    profile_response() = nothing
    U.quantity(::typeof(profile_response)) = U.Quantity{:profile_response}()
    U.native_unit(::U.Quantity{:profile_response}) = U.units(:base, :ohm)
    U.display_unit(::U.Quantity{:profile_response}) = U.units(:milli, :ohm)
    U.label(::U.Quantity{:profile_response}) = "Profile response"
    U.symbol(::U.Quantity{:profile_response}) = "u"

    @test quantity(profile_response) isa U.Quantity{:profile_response}
    @test label(profile_response) == "Profile response"
    @test symbol(profile_response) == "u"
    @test native_unit(profile_response) == U.units(:base, :ohm)
    @test native_unit(profile_response, :pul) == U.units(:base, :ohm)
    @test display_unit(profile_response) == U.units(:milli, :ohm)
    @test display_unit(profile_response, :pul) == U.units(:milli, :ohm)
    @test scale_factor(profile_response, U.units(:milli, :ohm)) == 1_000.0
    @test scale_factor(profile_response, :pul, U.units(:milli, :ohm)) == 1_000.0

    unregistered_selector() = nothing
    @test_throws MethodError quantity(unregistered_selector)
    @test_throws MethodError label(unregistered_selector)
    @test_throws MethodError display_unit(unregistered_selector, :pul)
end

@testitem "Units / locked public vocabulary" tags = [:unit] begin
    const U = LineCableModels.Units

    for name in (
        :Unit, :UnitExpr, :Quantity, :units, :quantity, :native_unit,
        :display_unit, :scale_factor, :label, :symbol
    )
        @test name in names(U)
    end
    @test isbitstype(U.Quantity{:series_resistance})
    @test fieldcount(U.Quantity{:series_resistance}) == 0
    @test !isdefined(U, :QuantityTag)
    for name in (:quantity, :native_unit, :display_unit, :scale_factor, :label, :symbol)
        @test name in names(LineCableModels)
        @test getfield(LineCableModels, name) === getfield(U, name)
    end
    for name in (:Quantity, :Unit, :UnitExpr)
        @test name ∉ names(LineCableModels)
        @test !isdefined(LineCableModels, name)
    end
    @test !isdefined(LineCableModels, :QuantityTag)
    @test !isdefined(LineCableModels, :QuantityUnits)
    @test !isdefined(U, :default_unit)
    @test !isdefined(U, :get_label)
    @test !isdefined(U, :get_symbol)
    @test !isdefined(U, :line_component_quantity)
    @test !isdefined(U, :line_component_unit)
end
