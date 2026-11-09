# An earth-return identifier that the test declares, whose formula gives neither quantity.
function LineCableModels.defaults(::Type{LineCableModels.Engine.EarthReturn.Formula{:declared}})
    (;)
end

@testset "Engine formulations" begin
    Engine = LineCableModels.Engine
    backend = LineCableModels.backend

    @testset "every formula slot accepts its three forms" begin
        for (constructor, slots) in (
            (CoaxialFormulation,
                (
                    (:internal_impedance, Engine.InternalImpedance,
                        (:schelkunoff1934, :wedepohl1973)),
                    (:insulation_impedance, Engine.InsulationImpedance, (:ametani1980,)),
                    (:insulation_admittance,
                        Engine.InsulationAdmittance, (:lossless, :lossy)),
                    (:semicon_admittance, Engine.SemiconAdmittance, (:lossless, :lossy)),
                    (:earth_impedance, Engine.EarthReturn,
                        (:martins2026, :gary1976, :wedepohl1973, :saad1996, :lucca1994,
                            :ametani2009, :carson1926, :pollaczek1926, :papadopoulos2010,
                            :wise1934, :xue2018)),
                    (:earth_admittance, Engine.EarthReturn,
                        (:martins2026, :ideal, :papadopoulos2010, :wise1948, :xue2018)),
                    (:shunt, Engine.Shunt, (:series,)),
                    (:resistivity, Engine.Resistivity, (:linear,)),
                    (:earth_properties, Engine.EarthProperties,
                        (:constant, :alipio2014, :cigre2019, :datsios2019, :longmire1975,
                            :messier1985, :portela1999, :scott1967, :visacro1987, :visacro2012)),
                    (:equivalent_earth, Engine.EarthProperties,
                        (:bottommost, :martins2020, :xue2021)))),
            (FEM,
                (
                    (:insulation_admittance,
                        Engine.InsulationAdmittance, (:lossless, :lossy)),
                    (:semicon_admittance, Engine.SemiconAdmittance, (:lossless, :lossy)),
                    (:resistivity, Engine.Resistivity, (:linear,)),
                    (:earth_properties, Engine.EarthProperties,
                        (:constant, :alipio2014, :cigre2019, :datsios2019, :longmire1975,
                            :messier1985, :portela1999, :scott1967, :visacro1987, :visacro2012)),
                    (:equivalent_earth, Engine.EarthProperties,
                        (:bottommost, :martins2020, :xue2021)))))
            for (slot, Family, identifiers) in slots, id in identifiers

                expected = constructor(; slot => id)
                @test constructor(; slot => (id => (;))) == expected
                @test constructor(; slot => Family.Formula(Val(id))) == expected
                if slot in (:earth_impedance, :earth_admittance)
                    record = constructor(; slot => (earth = id,))
                    @test getproperty(record, slot).earth ==
                          getproperty(expected, slot).earth
                end
            end
        end
        section = (integral = (rtol = 1e-10,),)
        law = Engine.EarthReturn.Formula(Val(:martins2026); section...)
        formulation = CoaxialFormulation(; earth_impedance = :martins2026 => section,
            earth_admittance = :martins2026 => section)
        @test formulation == CoaxialFormulation(; earth_impedance = law,
            earth_admittance = law)
        @test formulation.earth_impedance.earth.options.integral.rtol == 1e-10
        @test formulation.earth_admittance.earth.options.integral.rtol == 1e-10
    end

    @testset "unknown formula identifiers name the available formulas" begin
        for (constructor, slot, Family) in (
            (CoaxialFormulation, :internal_impedance, Engine.InternalImpedance),
            (CoaxialFormulation, :insulation_impedance, Engine.InsulationImpedance),
            (CoaxialFormulation, :insulation_admittance, Engine.InsulationAdmittance),
            (CoaxialFormulation, :semicon_admittance, Engine.SemiconAdmittance),
            (CoaxialFormulation, :shunt, Engine.Shunt),
            (CoaxialFormulation, :earth_properties, Engine.EarthProperties),
            (FEM, :insulation_admittance, Engine.InsulationAdmittance))
            @test_throws ArgumentError LineCableModels.select(constructor, Val(slot), :unknown)
            @test_throws "formula :unknown; the formulas are $(Family.FORMULAS)" LineCableModels.select(
                constructor, Val(slot), :unknown)
        end
    end

    @testset "a backend runs each formulation type" begin
        @test backend(CoaxialFormulation) === LineCableModels.LCMBackend()
        @test backend(typeof(FEM())) === LineCableModels.GmshBackend()
    end

    @testset "the coaxial formulation and its earth slots" begin
        formulation = CoaxialFormulation()
        @test (formulation.internal_impedance, formulation.insulation_impedance,
            formulation.insulation_admittance, formulation.semicon_admittance,
            formulation.earth_impedance.air, formulation.shunt, formulation.resistivity,
            formulation.earth_properties) ===
              (Engine.InternalImpedance.Formula(Val(:schelkunoff1934)),
            Engine.InsulationImpedance.Formula(Val(:ametani1980)),
            Engine.InsulationAdmittance.Formula(Val(:lossless)),
            Engine.SemiconAdmittance.Formula(Val(:lossless)),
            Engine.EarthReturn.Formula(Val(:martins2026)), Engine.Shunt.Formula(Val(:series)),
            Engine.Resistivity.Formula(Val(:linear)),
            Engine.EarthProperties.Formula(Val(:constant)))
        @test (formulation.options.Γ, formulation.options.bundle, formulation.options.kron,
            formulation.options.transposition) === (complex(0.0), true, true, false)
        # A formula of one identifier in both earth slots is one formula.
        section = (integral = (rtol = 1e-10,),)
        shared = CoaxialFormulation(; earth_impedance = :martins2026 => section,
            earth_admittance = :martins2026 => section)
        @test length(CoaxialFormulation().formulas) == 9 &&
              count(==(formulation.earth_impedance.air), CoaxialFormulation().formulas) == 1
        other = Engine.EarthReturn.Formula(Val(:martins2026); section...)
        @test_throws ArgumentError(
            "the earth slots hold two different :martins2026 formulas; a formula that gives Z and P " *
            "together is one formula in both slots") CoaxialFormulation(
            formulation.internal_impedance, formulation.insulation_impedance,
            formulation.insulation_admittance, formulation.semicon_admittance,
            formulation.earth_impedance,
            (air = other, earth = other, mixed = other), formulation.shunt,
            formulation.resistivity, formulation.earth_properties, formulation.equivalent_earth,
            formulation.options)
        @test_throws "two different :martins2026 formulas" CoaxialFormulation(;
            earth_impedance = :martins2026 => section)
        @test_throws DomainError CoaxialFormulation(; earth_impedance = :declared)
        # The default of each slot, by omission and by name, is the formula the slot names.
        explicit = CoaxialFormulation(; internal_impedance = :schelkunoff1934,
            insulation_impedance = :ametani1980, insulation_admittance = :lossless,
            semicon_admittance = :lossless, earth_impedance = :martins2026,
            earth_admittance = :martins2026, shunt = :series, resistivity = :linear,
            earth_properties = :constant, equivalent_earth = :bottommost)
        routed = CoaxialFormulation(; internal_impedance = :default,
            insulation_impedance = :default, insulation_admittance = :default,
            semicon_admittance = :default, earth_impedance = :default, earth_admittance = :default,
            shunt = :default, resistivity = :default, earth_properties = :default,
            equivalent_earth = :default)
        @test CoaxialFormulation() === explicit
        @test CoaxialFormulation() === routed
        beside = CoaxialFormulation(; earth_impedance = :default, earth_admittance = :ideal)
        martins = Engine.EarthReturn.Formula(Val(:martins2026))
        ideal = Engine.EarthReturn.Formula(Val(:ideal))
        @test beside.earth_impedance === (air = martins, earth = martins, mixed = martins)
        @test beside.earth_admittance === (air = ideal, earth = ideal, mixed = ideal)
        @test !occursin("default", string(typeof(beside)))
        # A selection per domain leaves the other domains to the default. The default takes
        # options too.
        @test CoaxialFormulation(; earth_impedance = (air = :martins2026,)) ===
              CoaxialFormulation(; earth_impedance = :martins2026)
        @test CoaxialFormulation(; earth_impedance = :default => section,
            earth_admittance = :default => section) === shared
        @test CoaxialFormulation(; earth_admittance = ideal) === beside
        # A law describes a layer, and a model reduces a stack.
        @test_throws "describes no layer" CoaxialFormulation(;
            earth_properties = :bottommost)
        @test_throws "reduces no earth model" CoaxialFormulation(;
            equivalent_earth = :constant)
        @test CoaxialFormulation(; Γ = [0.1im]).options.Γ == [0.1im]
        # A formula that does not read Γ holds at Γ = 0, and a computation at another Γ warns.
        # A Γ per frequency presents the range of its moduli.
        @test CoaxialFormulation(; Γ = 1.0e-3im).options.Γ == 1.0e-3im
        scan = LineParametersProblem(Dict(REFERENCES)["overhead_bare_wires"].system,
            Dict(REFERENCES)["overhead_bare_wires"].earth; frequencies = [50.0, 60.0])
        ideal = CoaxialFormulation(; earth_admittance = :ideal, Γ = 1.0e-3im)
        @test_logs (:warn,
            "the formula :ideal holds for Γ in 0.0 .. 0.0; " *
            "the computation has Γ in 0.001 .. 0.001") LineCableModels.validate(
            LineCableModels.domain(scan, ideal.options), ideal)
        ideal = CoaxialFormulation(;
            earth_admittance = (air = :ideal,), Γ = [0.0im, 1.0e-3im])
        @test_logs (:warn,
            "the formula :ideal holds for Γ in 0.0 .. 0.0; " *
            "the computation has Γ in 0.0 .. 0.001") LineCableModels.validate(
            LineCableModels.domain(scan, ideal.options), ideal)
        @test CoaxialFormulation(;
            earth_admittance = :ideal, Γ = [0.0im]).options.Γ == [0.0im]
        @test_throws DomainError CoaxialFormulation(; Γ = complex(Inf))
    end

    @testset "the finite-element formulation and its options" begin
        fem = FEM()
        @test fem.formulas === (fem.insulation_admittance, fem.semicon_admittance,
            fem.resistivity, fem.earth_properties, fem.equivalent_earth)
        for selection in ((internal_impedance = :schelkunoff1934,),
            (insulation_impedance = :ametani1980,),
            (earth_impedance = :martins2026,), (earth_admittance = :martins2026,), (shunt = :series,))
            @test CoaxialFormulation(; selection...) isa CoaxialFormulation
            @test_throws "unknown options" FEM(; selection...)
        end
        @test (fem.insulation_admittance, fem.semicon_admittance, fem.resistivity,
            fem.earth_properties) === (Engine.InsulationAdmittance.Formula(Val(:lossless)),
            Engine.SemiconAdmittance.Formula(Val(:lossless)),
            Engine.Resistivity.Formula(Val(:linear)),
            Engine.EarthProperties.Formula(Val(:constant)))
        @test FEM(; insulation_admittance = :default, semicon_admittance = :default,
            resistivity = :default, earth_properties = :default,
            equivalent_earth = :default) === fem
        @test_throws "describes no layer" FEM(; earth_properties = :bottommost)
        @test FEM(; equivalent_earth = :bottommost) === fem
        @test_throws "reduces no earth model" FEM(; equivalent_earth = :constant)
        for material in
            (Material(:insulator, 1.0e14, 2.3), Material(:semicon, 1000.0, 1000.0),
            Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393))
            @test Engine.admittivity(FEM(), material, 50.0, 70.0) ==
                  Engine.admittivity(CoaxialFormulation(), material, 50.0, 70.0)
        end
        @test FEM(; Γ = [complex(0.0, 1.0)]).options.Γ == [complex(0.0, 1.0)]
        @test_throws DomainError FEM(; Γ = complex(NaN))
        @test_throws DomainError FEM(; Γ = ComplexF64[])
    end

    @testset "the admittivity of a medium in one home" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        semicon = Material(:semicon, 1000.0, 1000.0)
        s = complex(0.0, 2π * 50.0)
        # The air, an earth layer, an equivalent earth of negative permittivity and a metal.
        for (rho, eps_r) in
            ((Inf, 1.0), (100.0, 10.0), (62.05, -48961.76), (1.7241e-8, 1.0))
            @test Engine.admittivity(rho, eps_r, s) === inv(rho) + s * (Engine.ε₀ * eps_r)
        end
        layer = EarthLayer(100.0, 10.0, 1.0, 5.0)
        @test Engine.admittivity(layer, s) === Engine.admittivity(layer.rho, layer.eps_r, s)
        for formulation in (CoaxialFormulation(), FEM())
            ρ = formulation.resistivity(copper, 70.0)
            for f in (0.0, 50.0, 1.0e6)
                κ = Engine.admittivity(formulation, copper, f, 70.0)
                @test κ === Engine.admittivity(ρ, copper.eps_r, complex(0.0, 2π * f))
                @test real(κ) === inv(ρ)
                for θ in (20.0, 90.0)
                    @test Engine.admittivity(formulation, polymer, f, θ) ===
                          formulation.insulation_admittance(polymer, f)
                    @test Engine.admittivity(formulation, semicon, f, θ) ===
                          formulation.semicon_admittance(semicon, f)
                end
            end
            @test Engine.admittivity(formulation, copper, 0.0, 70.0) === complex(inv(ρ))
            @test real(Engine.admittivity(formulation, copper, 50.0, 20.0)) >
                  real(Engine.admittivity(formulation, copper, 50.0, 70.0))
        end
        @test @inferred(Engine.admittivity(100.0, 10.0, s)) isa ComplexF64
        @test @inferred(Engine.admittivity(layer, s)) isa ComplexF64
        @test @inferred(Engine.admittivity(FEM(), copper, 50.0, 70.0)) isa ComplexF64
        @test only(map((g, l, z) -> @allocated(g(l, z)), (Engine.admittivity,), (layer,), (s,))) ==
              0
        @test only(map((g, m, c, f, θ) -> @allocated(g(m, c, f, θ)),
            (Engine.admittivity,), (FEM(),), (copper,), (50.0,), (70.0,))) == 0
    end

    @testset "a block lowers to the constructor" begin
        domains = (air = :gary1976, earth = :saad1996, mixed = :lucca1994)
        @test (@formulation begin
            earth_impedance = @air :gary1976 @earth :saad1996 @mixed :lucca1994
            earth_admittance = :ideal
            transposition = true
        end) === CoaxialFormulation(;
            earth_impedance = domains, earth_admittance = :ideal, transposition = true)
        @test (@formulation begin
            earth_impedance = @mixed :lucca1994 @air :gary1976 @earth :saad1996
        end) === CoaxialFormulation(; earth_impedance = domains)
        @test (@formulation begin
            earth_admittance = @air :ideal
        end) === CoaxialFormulation(; earth_admittance = (air = :ideal,))
        section = (integral = (rtol = 1e-6,),)
        @test (@formulation begin
            earth_impedance = @air :martins2026 => section @earth :martins2026 => section @mixed :martins2026 => section
            earth_admittance = :martins2026 => section
        end) === CoaxialFormulation(; earth_impedance = :martins2026 => section,
            earth_admittance = :martins2026 => section)
        saad = Engine.EarthReturn.Formula(Val(:saad1996))
        @test (@formulation begin
            earth_impedance = @earth saad
        end) === CoaxialFormulation(; earth_impedance = (earth = saad,))
        @test (@formulation FEM begin
            earth_properties = :alipio2014
            Γ = 1e-3im
        end) === FEM(; earth_properties = :alipio2014, Γ = 1e-3im)
        @test (@formulation begin end) === CoaxialFormulation()

        @test_throws "unknown options (:pipe,)" @formulation begin
            pipe = :none
        end
        @test_throws "unknown options (:earth_impedance,)" @formulation FEM begin
            earth_impedance = :saad1996
        end
        @test (@formulation begin
            earth_impedance = @air :gary1976 @earth :saad1996 @mixed :lucca1994
            earth_admittance = :ideal
            Γ = 1e-3im
        end) === CoaxialFormulation(;
            earth_impedance = (air = :gary1976, earth = :saad1996, mixed = :lucca1994),
            earth_admittance = :ideal, Γ = 1e-3im)
        @test_throws "the earth slots hold two different :martins2026 formulas" @formulation begin
            earth_impedance = @air :martins2026 => section
        end
        @test_throws "invalid keyword argument syntax" eval(:(@formulation begin
            earth_admittance:ideal
        end))
    end

    @testset "the finite-element computation needs the Gmsh extension" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        wire = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        system = build(
            CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        @test_throws "using Gmsh, getdp_jll" compute(problem, FEM())
    end
end
