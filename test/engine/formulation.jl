# An earth-return identifier that the test declares, whose formula gives neither quantity.
function LineCableModels.defaults(::Type{LineCableModels.Engine.EarthReturn.Formula{:declared}})
    (;)
end

@testset "Engine formulations" begin
    Engine = LineCableModels.Engine
    backend = LineCableModels.backend

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
        @test_throws DomainError CoaxialFormulation(; Γ = complex(Inf))
    end

    @testset "the finite-element formulation and its options" begin
        fem = FEM()
        @test (fem.insulation_admittance, fem.semicon_admittance, fem.resistivity,
            fem.earth_properties) === (Engine.InsulationAdmittance.Formula(Val(:lossless)),
            Engine.SemiconAdmittance.Formula(Val(:lossless)),
            Engine.Resistivity.Formula(Val(:linear)),
            Engine.EarthProperties.Formula(Val(:constant)))
        @test FEM(; insulation_admittance = :default, semicon_admittance = :default,
            resistivity = :default, earth_properties = :default) === fem
        @test_throws "describes no layer" FEM(; earth_properties = :bottommost)
        for material in
            (Material(:insulator, 1.0e14, 2.3), Material(:semicon, 1000.0, 1000.0))
            @test Engine.admittivity(FEM(), material, 50.0) ==
                  Engine.admittivity(CoaxialFormulation(), material, 50.0)
        end
        @test FEM(; Γ = [complex(0.0, 1.0)]).options.Γ == [complex(0.0, 1.0)]
        @test_throws DomainError FEM(; Γ = complex(NaN))
        @test_throws DomainError FEM(; Γ = ComplexF64[])
    end

    @testset "the finite-element computation needs the Gmsh extension" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        wire = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)))
        system = build(
            CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        @test_throws "using Gmsh" compute(problem, FEM())
    end
end
