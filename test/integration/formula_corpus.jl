@testitem "Engine / literature formula corpus / registered routes" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const EN=LineCableModels.Engine
    const EI=EN.EarthImpedance
    const EA=EN.EarthAdmittance
    const II=EN.InternalImpedance
    const IZ=EN.InsulationImpedance
    const IA=EN.InsulationAdmittance
    const SA=EN.SemiconAdmittance

    expected_internal=(
        :Ametani2004,
        :AmetaniFuse1992,
        :Schelkunoff1934,
        :WedepohlWilcox1973,
        :Zhao2020
    )
    expected_earth_impedance=(
        :AlvaradoBetancourt1983,
        :Ametani1974,
        :Ametani2009,
        :Bridges1995,
        :Carson1926,
        :Gary1976,
        :Lucca1994,
        :Magalhaes2018,
        :MartinsBritto2024,
        :Nakagawa1973,
        :Noda2006,
        :Papadopoulos2009,
        :Papadopoulos2010,
        :Papadopoulos2011,
        :Petrache2005,
        :Pettersson1994,
        :Pollaczek1926,
        :Saad1996,
        :Sunde1968,
        :Theethayi2007,
        :Theodoulidis2015,
        :Tsiamitros2008,
        :Vance1978,
        :WedepohlWilcox1973,
        :Wise1934,
        :Xue2018
    )
    expected_earth_admittance=(
        :Ametani2021,
        :IdealGround,
        :Magalhaes2018,
        :MartinsBritto2024,
        :Papadopoulos2009,
        :Papadopoulos2010,
        :Papadopoulos2011,
        :Pettersson1994,
        :Pollaczek1926,
        :Theethayi2007,
        :Wise1948,
        :Xue2018,
        :Xue2021
    )

    @test II.formulas() == expected_internal
    @test EI.formulas() == expected_earth_impedance
    @test EA.formulas() == expected_earth_admittance

    allowed_templates=(
        EI => Set((:earth_impedance, :propagation_constant)),
        EA => Set((
            :earth_potential_coefficient,
            :earth_impedance,
            :propagation_constant
        ))
    )
    for (catalogue, templates) in allowed_templates, id in catalogue.formulas()
        formula=catalogue.Formula(id)
        for route in values(catalogue.routes(formula))
            @test route isa EN.FormulaMethod
            @test typeof(route).parameters[1] === id
            @test nameof(route.method) in templates
            @test all(argument -> argument isa Val, route.arguments)
        end
    end

    @test all(id -> !isempty(description(II.Formula(id))), II.formulas())
    @test all(id -> !isempty(description(EI.Formula(id))), EI.formulas())
    @test all(id -> !isempty(description(EA.Formula(id))), EA.formulas())
    @test all(id -> !isempty(description(IZ.Formula(id))), IZ.formulas())
    @test all(id -> !isempty(description(IA.Formula(id))), IA.formulas())
    @test all(id -> !isempty(description(SA.Formula(id))), SA.formulas())

    epsilon0=8.8541878128e-12
    mu0=4π*1.0e-7
    rho=Float64[Inf, 0.1]
    epsilon=Float64[epsilon0, epsilon0]
    mu=Float64[mu0, mu0]
    s=ComplexF64(0.0, 2π*50.0)
    pair=EarthPair(1, 2, (-1.0, -1.2), 1.0, (2, 2))
    formula=EI.Formula(:Papadopoulos2010)
    functor=@inferred formula(rho, epsilon, mu, s, nothing)
    @test isfinite(@inferred functor(Val(:mutual), pair))
    admittance=EA.Formula(:Papadopoulos2010)
    admittance_functor=@inferred admittance(rho, epsilon, mu, s, nothing)
    @test isfinite(@inferred admittance_functor(Val(:mutual), pair))
end

@testitem "Engine / literature formula corpus / recovered closed forms" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const EI=LineCableModels.Engine.EarthImpedance
    const EA=LineCableModels.Engine.EarthAdmittance

    epsilon0=8.8541878128e-12
    mu0=4π*1.0e-7
    rho=Float64[Inf, 1000.0]
    epsilon=Float64[epsilon0, 10epsilon0]
    mu=Float64[mu0, mu0]
    s=ComplexF64(0.0, 2π*1.0e5)
    overhead=EarthPair(1, 2, (10.0, 12.0), 5.0, (1, 1))
    high_angle=EarthPair(1, 2, (10.0, 10.0), 30.0, (1, 1))

    carson=EI.Formula(:Carson1926)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)
    wise_z=EI.Formula(:Wise1934)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)
    wise_p=EA.Formula(:Wise1948)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)
    alvarado=EI.Formula(:AlvaradoBetancourt1983)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)
    noda=EI.Formula(:Noda2006)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)
    pettersson_z=EI.Formula(:Pettersson1994)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)
    pettersson_p=EA.Formula(:Pettersson1994)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), overhead)

    @test abs(alvarado / carson - 1) < 0.003
    @test abs(noda / carson - 1) < 0.001
    @test abs(pettersson_z / wise_z - 1) < 0.02
    @test abs(pettersson_p / wise_p - 1) < 0.03
    @test isfinite(EI.Formula(:Noda2006)(
        rho, epsilon, mu, s, nothing
    )(Val(:mutual), high_angle))

    mixed=EarthPair(1, 2, (10.0, -1.0), 100.0, (1, 2))
    s_mixed=ComplexF64(0.0, 2π*50.0)
    exact=EI.Formula(:Pollaczek1926)(
        rho, epsilon, mu, s_mixed, nothing
    )(Val(:mutual), mixed)
    ametani=EI.Formula(:Ametani2009)(
        rho, epsilon, mu, s_mixed, nothing
    )(Val(:mutual), mixed)
    lucca=EI.Formula(:Lucca1994)(
        rho, epsilon, mu, s_mixed, nothing
    )(Val(:mutual), mixed)
    @test abs(ametani / exact - 1) < 0.04
    @test abs(lucca / exact - 1) < 0.01

    function evaluate(formula, rho, epsilon, mu, s, pair)
        functor=@inferred formula(rho, epsilon, mu, s, nothing)
        return @inferred functor(Val(:mutual), pair)
    end

    for T in (Float32, Float64, BigFloat)
        epsilon0_t=T(8.8541878128e-12)
        mu0_t=T(4)*T(π)*T(1.0e-7)
        rho_t=T[Inf, 100]
        epsilon_t=T[epsilon0_t, 10epsilon0_t]
        mu_t=T[mu0_t, mu0_t]
        s_t=Complex{T}(zero(T), T(2)*T(π)*T(50))
        overhead_t=EarthPair(
            1, 2, (T(10), T(12)), T(5), (1, 1)
        )
        mixed_t=EarthPair(
            1, 2, (T(10), T(-1)), T(5), (1, 2)
        )
        underground_t=EarthPair(
            1, 2, (T(-1), T(-1.2)), T(1), (2, 2)
        )
        for formula in (
            EI.Formula(:Carson1926),
            EI.Formula(:AlvaradoBetancourt1983),
            EI.Formula(:Noda2006),
            EI.Formula(:Pettersson1994),
            EA.Formula(:Pettersson1994)
        )
            @test evaluate(
                formula, rho_t, epsilon_t, mu_t, s_t, overhead_t
            ) isa Complex{T}
        end
        for formula in (
            EI.Formula(:Ametani2009),
            EI.Formula(:Lucca1994)
        )
            @test evaluate(
                formula, rho_t, epsilon_t, mu_t, s_t, mixed_t
            ) isa Complex{T}
        end
        @test evaluate(
            EI.Formula(:Saad1996),
            rho_t,
            epsilon_t,
            mu_t,
            s_t,
            underground_t
        ) isa Complex{T}
    end
end

@testitem "Engine / literature formula corpus / two-wire numerical smoke" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const EN=LineCableModels.Engine

    function two_bare_wires(;
            y = -1.0,
            frequencies = Float64[1.0, 50.0, 1.0e6],
            earth = Earth(rho = 0.1, eps_r = 1.0, mu_r = 1.0)
    )
        conductor=Material(
            kind = :conductor,
            rho = eps(Float64),
            eps_r = 1.0,
            mu_r = 1.0,
            T0 = 20.0,
            alpha = 0.0
        )
        insulation=Material(
            kind = :insulator,
            rho = 1.97e14,
            eps_r = 2.3,
            mu_r = 1.0,
            T0 = 20.0,
            alpha = 0.0
        )
        design=build(
            CableDesign,
            "two_bare_wires",
            Stack(
                Group(:core, Region(:core_metal, Disk(0.0425), conductor)),
                Region(:core_insulation, Shell(1.0e-3), insulation)
            )
        )
        system=build(
            LineCableSystem,
            [design, design],
            [Pose2(0.0, y), Pose2(1.0, y)];
            connections = [Dict(:core=>1), Dict(:core=>2)],
            system_id = "two_bare_wires",
            line_length = 1.0
        )
        return LineParametersProblem(
            system;
            temperature = 20.0,
            earth_props = earth,
            frequencies
        )
    end

    function finite_reciprocal(result)
        return all(isfinite, result.Z.values)&&
               all(isfinite, result.Y.values)&&
               all(
            index->isapprox(
                result.Z.values[1, 2, index],
                result.Z.values[2, 1, index];
                rtol = 1.0e-11
            ),
            eachindex(result.f))&&
               all(
            index->isapprox(
                result.Y.values[1, 2, index],
                result.Y.values[2, 1, index];
                rtol = 1.0e-11
            ),
            eachindex(result.f))
    end

    underground=two_bare_wires()
    overhead=two_bare_wires(y = 1.0)
    underground_impedance=(
        :Ametani2009,
        :Lucca1994,
        :Magalhaes2018,
        :Papadopoulos2010,
        :Petrache2005,
        :Pollaczek1926,
        :Saad1996,
        :Theethayi2007,
        :Tsiamitros2008,
        :WedepohlWilcox1973,
        :Xue2018
    )
    overhead_impedance=(
        :AlvaradoBetancourt1983,
        :Ametani2009,
        :Carson1926,
        :Gary1976,
        :Lucca1994,
        :Noda2006,
        :Pettersson1994,
        :Pollaczek1926,
        :Sunde1968,
        :Theodoulidis2015,
        :Tsiamitros2008,
        :Wise1934
    )
    underground_admittance=(
        :Magalhaes2018,
        :Papadopoulos2010,
        :Pollaczek1926,
        :Theethayi2007,
        :Xue2018,
        :Xue2021
    )
    overhead_admittance=(
        :Ametani2021,
        :Pettersson1994,
        :Pollaczek1926,
        :Wise1948
    )

    for identifier in underground_impedance
        result=compute(underground, Formulation(
            earth_impedance = identifier,
            earth_admittance = :IdealGround
        ))
        @test finite_reciprocal(result)
    end
    for identifier in overhead_impedance
        result=compute(overhead, Formulation(
            earth_impedance = identifier,
            earth_admittance = :IdealGround
        ))
        @test finite_reciprocal(result)
    end
    for identifier in underground_admittance
        result=compute(underground, Formulation(
            earth_impedance = :Petrache2005,
            earth_admittance = identifier
        ))
        @test finite_reciprocal(result)
    end
    for identifier in overhead_admittance
        result=compute(overhead,
            Formulation(
                earth_impedance = :Pollaczek1926,
                earth_admittance = identifier
            ))
        @test finite_reciprocal(result)
    end

    for identifier in (:Schelkunoff1934, :AmetaniFuse1992)
        result=compute(underground,
            Formulation(
                internal_impedance = identifier,
                earth_impedance = :Petrache2005,
                earth_admittance = :IdealGround
            ))
        @test finite_reciprocal(result)
    end
    @test finite_reciprocal(compute(underground,
        Formulation(
            insulation_impedance = :Ametani1980,
            earth_impedance = :Petrache2005,
            earth_admittance = :IdealGround
        )))
    for identifier in (:Gustavsen2013, :Marti2001)
        result=compute(underground,
            Formulation(
                insulation_admittance = identifier,
                earth_impedance = :Petrache2005,
                earth_admittance = :IdealGround
            ))
        @test finite_reciprocal(result)
    end
end

@testitem "Engine / literature formula corpus / stratified and mixed limits" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    const EI=LineCableModels.Engine.EarthImpedance
    const EA=LineCableModels.Engine.EarthAdmittance

    epsilon0=8.8541878128e-12
    mu0=4π*1.0e-7
    s=ComplexF64(0.0, 2π*50.0)
    rho2=Float64[Inf, 0.1]
    epsilon2=Float64[epsilon0, epsilon0]
    mu2=Float64[mu0, mu0]
    rho3=Float64[Inf, 0.1, 0.1]
    epsilon3=Float64[epsilon0, epsilon0, epsilon0]
    mu3=Float64[mu0, mu0, mu0]
    thickness=Float64[Inf, 4.0, Inf]
    homogeneous_thickness=Float64[Inf, Inf]
    overhead=EarthPair(1, 2, (1.0, 1.2), 1.0, (1, 1))
    underground=EarthPair(1, 2, (-1.0, -1.2), 1.0, (2, 2))
    vertical_underground=EarthPair(1, 2, (-1.0, -1.2), 0.0, (2, 2))
    mixed=EarthPair(1, 2, (1.0, -1.0), 1.0, (1, 2))

    sunde=EI.Formula(:Sunde1968)(
        rho2, epsilon2, mu2, s, nothing
    )(Val(:mutual), overhead)
    wise=EI.Formula(:Wise1934)(
        rho2, epsilon2, mu2, s, nothing
    )(Val(:mutual), overhead)
    sunde_two=EI.Formula(:Sunde1968)(
        rho3, epsilon3, mu3, s, nothing, nothing, thickness
    )(Val(:mutual), overhead)
    @test sunde_two ≈ sunde rtol=1.0e-10
    for identifier in (
        :Ametani1974,
        :Nakagawa1973,
        :Papadopoulos2009
    )
        value=EI.Formula(identifier)(
            rho3, epsilon3, mu3, s, nothing, nothing, thickness
        )(Val(:mutual), overhead)
        @test value ≈ wise rtol=1.0e-10
    end

    homogeneous_z=EI.Formula(:Papadopoulos2010)(
        rho2, epsilon2, mu2, s, nothing
    )(Val(:mutual), underground)
    stratified_z=EI.Formula(:Papadopoulos2011)(
        rho3, epsilon3, mu3, s, nothing, nothing, thickness
    )(Val(:mutual), underground)
    @test stratified_z ≈ homogeneous_z rtol=1.0e-9
    homogeneous_y=EA.Formula(:Papadopoulos2010)(
        rho2, epsilon2, mu2, s, nothing
    )(Val(:mutual), underground)
    stratified_y=EA.Formula(:Papadopoulos2011)(
        rho3, epsilon3, mu3, s, nothing, nothing, thickness
    )(Val(:mutual), underground)
    @test stratified_y ≈ homogeneous_y rtol=1.0e-8
    overhead_y=EA.Formula(:Wise1948)(
        rho2, epsilon2, mu2, s, nothing
    )(Val(:mutual), overhead)
    stratified_overhead_y=EA.Formula(:Papadopoulos2009)(
        rho3, epsilon3, mu3, s, nothing, nothing, thickness
    )(Val(:mutual), overhead)
    @test stratified_overhead_y ≈ overhead_y rtol=1.0e-9

    rho_mixed=Float64[Inf, 0.1]
    epsilon_mixed=Float64[epsilon0, epsilon0]
    mu_mixed=Float64[mu0, mu0]
    radial_self=EarthPair(1, 1, (-1.0, -1.0), 0.0425, (2, 2))
    for identifier in (:Bridges1995, :Vance1978)
        radial=EI.Formula(identifier)(
            rho_mixed, epsilon_mixed, mu_mixed, s, nothing
        )
        @test isfinite(radial(Val(:self), radial_self))
        @test_throws ArgumentError radial(Val(:mutual), underground)
    end
    for identifier in (
            :Petrache2005,
            :Saad1996,
            :Theethayi2007,
            :WedepohlWilcox1973
        )
        radial=EI.Formula(identifier)(
            rho_mixed, epsilon_mixed, mu_mixed, s, nothing
        )
        @test_throws DomainError radial(Val(:mutual), vertical_underground)
    end
    for identifier in (:Theethayi2007, :Xue2021)
        radial=EA.Formula(identifier)(
            rho_mixed, epsilon_mixed, mu_mixed, s, nothing
        )
        @test_throws DomainError radial(Val(:mutual), vertical_underground)
    end
    for identifier in (
        :Ametani2009,
        :Lucca1994,
        :MartinsBritto2024,
        :Pollaczek1926
    )
        value=EI.Formula(identifier)(
            rho_mixed, epsilon_mixed, mu_mixed, s, nothing
        )(Val(:mutual), mixed)
        @test isfinite(value)
    end
    carson=EI.Formula(:Carson1926)(
        rho_mixed, epsilon_mixed, mu_mixed, s, nothing
    )(Val(:mutual), overhead)
    pollaczek_overhead=EI.Formula(:Pollaczek1926)(
        rho_mixed, epsilon_mixed, mu_mixed, s, nothing
    )(Val(:mutual), overhead)
    pollaczek_underground=EI.Formula(:Pollaczek1926)(
        rho_mixed, epsilon_mixed, mu_mixed, s, nothing
    )(Val(:mutual), underground)
    @test carson == pollaczek_overhead
    for identifier in (:Ametani2009, :Lucca1994)
        recipe=EI.Formula(identifier)(
            rho_mixed, epsilon_mixed, mu_mixed, s, nothing
        )
        @test recipe(Val(:mutual), overhead) == carson
        @test recipe(Val(:mutual), underground) == pollaczek_underground
    end
    for identifier in (:Pollaczek1926, :MartinsBritto2024)
        value=EA.Formula(identifier)(
            rho_mixed, epsilon_mixed, mu_mixed, s, nothing
        )(Val(:mutual), mixed)
        @test isfinite(value)
    end
    for owner in (EI, EA)
        formula=owner.Formula(:MartinsBritto2024)
        functor=formula(rho_mixed, epsilon_mixed, mu_mixed, s, nothing)
        @test all(isfinite,
            (
                functor(Val(:mutual), overhead),
                functor(Val(:mutual), underground),
                functor(Val(:mutual), mixed)
            ))
        @test all(
            name->hasproperty(owner.routes(formula), name),
            (:overhead, :underground, :mixed)
        )
    end
    xue=EA.Formula(:Xue2018)
    xue_functor=xue(rho_mixed, epsilon_mixed, mu_mixed, s, nothing)
    for name in (:infinite, :surface, :penetration)
        @test isfinite(getproperty(EA.routes(xue), name)(xue_functor, underground))
    end
    surface=EA.Formula(
        :Xue2018;
        self = EA.routes(xue).surface,
        mutual = EA.routes(xue).surface
    )
    @test surface(rho_mixed, epsilon_mixed, mu_mixed, s, nothing)(
        Val(:mutual), underground
    ) == EA.routes(xue).surface(xue_functor, underground)
    magalhaes=EA.Formula(:Magalhaes2018)(
        rho_mixed, epsilon_mixed, mu_mixed, s, nothing
    )
    for (kind, pair) in (
            (:mutual, underground),
            (:self, EarthPair(1, 1, (-1.0, -1.0), 0.0425, (2, 2)))
        )
        @test magalhaes(Val(kind), pair) ≈
              EA.routes(xue).surface(xue_functor, pair) rtol=1.0e-9
    end
    tsiamitros=@inferred EI.Formula(:Tsiamitros2008)(
        rho3, epsilon3, mu3, s, nothing, nothing, thickness
    )(Val(:mutual), mixed)
    @test isfinite(tsiamitros)

    homogeneous_tsiamitros=EI.Formula(:Tsiamitros2008)(
        rho2, epsilon2, mu2, s, nothing, nothing, homogeneous_thickness
    )
    stratified_tsiamitros=EI.Formula(:Tsiamitros2008)(
        rho3, epsilon3, mu3, s, nothing, nothing, thickness
    )
    for pair in (overhead, underground, mixed)
        homogeneous_value=homogeneous_tsiamitros(Val(:mutual), pair)
        stratified_value=stratified_tsiamitros(Val(:mutual), pair)
        @test stratified_value ≈ homogeneous_value rtol=1.0e-12
    end
    self_pair=EarthPair(1, 1, (-1.0, -1.0), 0.0425, (2, 2))
    @test homogeneous_tsiamitros(Val(:self), self_pair) ==
          homogeneous_tsiamitros(Val(:mutual), self_pair)

    for identifier in (
        :Schelkunoff1934,
        :WedepohlWilcox1973,
        :Zhao2020
    )
        functor=LineCableModels.Engine.InternalImpedance.Formula(identifier)(
            0.008, 0.010, 1.7241e-8, 1.0, s
        )
        @test all(isfinite, (
            functor(Val(:inner)),
            functor(Val(:outer)),
            functor(Val(:mutual))
        ))
    end
    aggregate=LineCableModels.Engine.InternalImpedance.Formula(
        :Ametani2004
    )(
        0.005, 0.010, 0.010, 0.0105,
        1.7241e-8, 1.0e-4, 1.0, 1.0, s
    )
    @test all(isfinite, (
        aggregate(Val(:inner)),
        aggregate(Val(:outer)),
        aggregate(Val(:mutual))
    ))
end
