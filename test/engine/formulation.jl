# An earth-return identifier that the test declares, whose formula gives neither quantity.
function LineCableModels.defaults(::Type{LineCableModels.Engine.EarthReturn.Formula{:declared}})
    (;)
end

@testset "Engine formulations" begin
    Engine = LineCableModels.Engine
    backend = LineCableModels.backend

    @testset "every formula slot accepts its three forms" begin
        for constructor in (CoaxialFormulation, FEM)
            formulation = constructor()
            for slot in fieldnames(typeof(formulation))[1:(end - 1)]
                law = getproperty(formulation, slot)
                family = parentmodule(typeof(law))
                id = only(family.FORMULAS)
                expected = constructor(; (slot => id,)...)
                @test constructor(; (slot => (id => (;)),)...) == expected
                @test constructor(; (slot => law,)...) == expected
            end
        end
        section = (integral = (rtol = 1e-10,),)
        law = Engine.EarthReturn.Formula(Val(:unified); section...)
        formulation = CoaxialFormulation(; earth_impedance = :unified => section,
            earth_admittance = :unified => section)
        @test formulation == CoaxialFormulation(; earth_impedance = law,
            earth_admittance = law)
        @test formulation.earth_impedance.options.integral.rtol == 1e-10
        @test formulation.earth_admittance.options.integral.rtol == 1e-10
    end

    @testset "unknown formula identifiers name the available formulas" begin
        for family in (Engine.InternalImpedance, Engine.InsulationImpedance,
            Engine.InsulationAdmittance, Engine.SemiconAdmittance, Engine.Shunt,
            Engine.EarthProperties)
            @test_throws ArgumentError LineCableModels.select(family.Formula, :unknown)
            @test_throws "formula :unknown; the formulas are $(family.FORMULAS)" LineCableModels.select(
                family.Formula, :unknown)
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
            formulation.earth_impedance, formulation.shunt, formulation.resistivity,
            formulation.earth_properties) ===
              (Engine.InternalImpedance.Formula(Val(:schelkunoff1934)),
            Engine.InsulationImpedance.Formula(Val(:ametani1980)),
            Engine.InsulationAdmittance.Formula(Val(:lossless)),
            Engine.SemiconAdmittance.Formula(Val(:lossless)),
            Engine.EarthReturn.Formula(Val(:unified)), Engine.Shunt.Formula(Val(:series)),
            Engine.Resistivity.Formula(Val(:linear)),
            Engine.EarthProperties.Formula(Val(:constant)))
        @test (formulation.options.Γ, formulation.options.bundle, formulation.options.kron,
            formulation.options.transposition) === (complex(0.0), true, true, false)
        # A formula of one identifier in both earth slots is one formula.
        section = (integral = (rtol = 1e-10,),)
        shared = CoaxialFormulation(; earth_impedance = :unified => section,
            earth_admittance = :unified => section)
        @test length(Engine.formulas(CoaxialFormulation())) == 8 &&
              count(==(formulation.earth_impedance), Engine.formulas(CoaxialFormulation())) ==
              1
        @test_throws ArgumentError(
            "the earth slots hold two different :unified formulas; a formula that gives Z and P " *
            "together is one formula in both slots") CoaxialFormulation(
            formulation.internal_impedance, formulation.insulation_impedance,
            formulation.insulation_admittance, formulation.semicon_admittance,
            formulation.earth_impedance,
            Engine.EarthReturn.Formula(Val(:unified); section...), formulation.shunt,
            formulation.resistivity, formulation.earth_properties, formulation.options)
        @test_throws "two different :unified formulas" CoaxialFormulation(;
            earth_impedance = :unified => section)
        @test_throws DomainError CoaxialFormulation(; earth_impedance = :declared)
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
        @test_throws "using Gmsh, getdp_jll" compute(problem, FEM())
    end
end
