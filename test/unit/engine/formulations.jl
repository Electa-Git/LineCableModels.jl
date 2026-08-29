@testitem "Engine / formula selector / owner-local lowering" tags=[:unit] setup=[
    UseEngineSupport,
] begin
    const EP=LineCableModels.EarthProps
    const EH=EP.EHEM
    const EN=LineCableModels.Engine

    inner=state->oftype(state.jω, 7)
    insulation_impedance_route=(r_in, r_ex, mu_r, s, values)->zero(s)
    insulation_admittance_route=(r_in, r_ex, rho, eps_r, s, values)->zero(s)
    formulation=@inferred Formulation(
        internal_impedance=formula(:Schelkunoff; inner),
        insulation_impedance=formula(
            :Lossless;
            route=insulation_impedance_route
        ),
        earth_impedance=formula(:Papadopoulos2010),
        insulation_admittance=formula(
            :Lossless;
            route=insulation_admittance_route
        ),
        earth_admittance=formula(:Papadopoulos2010),
        earth_properties=formula(:CIGRE2019; epsilon_infinity=10.0),
        equivalent_earth=formula(:Xue2021; order=:before, layer=2)
    )

    @test formula_id(formulation.methods.internal_impedance) === :Schelkunoff
    @test EN.InternalImpedance.routes(
        formulation.methods.internal_impedance
    ).inner === inner
    @test formula_id(formulation.methods.insulation_impedance) === :Lossless
    @test formulation.methods.insulation_impedance.route ===
          insulation_impedance_route
    @test formula_id(formulation.methods.earth_impedance) === :Papadopoulos2010
    @test formula_id(formulation.methods.insulation_admittance) === :Lossless
    @test formulation.methods.insulation_admittance.route ===
          insulation_admittance_route
    @test formula_id(formulation.methods.earth_admittance) === :Papadopoulos2010
    @test formula_id(formulation.methods.earth_properties) === :CIGRE2019
    @test EP.FD.assumptions(
        formulation.methods.earth_properties
    ).epsilon_infinity == 10.0
    @test formulation.methods.equivalent_earth isa EH.BeforeFD
    @test formula_id(EH.rule(formulation.methods.equivalent_earth)) === :Xue2021
    @test EH.assumptions(EH.rule(
        formulation.methods.equivalent_earth
    )).layer == 2

    defaults=@inferred Formulation()
    bare=Formulation(
        insulation_impedance=:Lossless,
        insulation_admittance=:ParallelRC
    )
    @test formula_id(defaults.methods.insulation_impedance) === :Lossless
    @test formula_id(defaults.methods.insulation_admittance) === :Lossless
    @test formula_id(bare.methods.insulation_impedance) === :Lossless
    @test formula_id(bare.methods.insulation_admittance) === :ParallelRC

    after_default=@inferred Formulation(equivalent_earth=formula(:Xue2021))
    after_explicit=@inferred Formulation(
        equivalent_earth=formula(:Xue2021; order=:after)
    )
    layer=@inferred Formulation(
        equivalent_earth=formula(:Layer; order=:before, layer=2)
    )
    @test after_default.methods.equivalent_earth isa EH.AfterFD
    @test after_explicit.methods.equivalent_earth isa EH.AfterFD
    @test layer.methods.equivalent_earth isa EH.BeforeFD
    @test EH.rule(layer.methods.equivalent_earth) == EH.Layer(2)

    selection=formula(:Xue2021; order=:before)
    @test formula_id(selection) === :Xue2021
    @test isconcretetype(typeof(selection))
    @test_throws ArgumentError formula(:Xue2021; order=:sideways)
    @test_throws ArgumentError Formulation(
        earth_properties=formula(:CIGRE2019; order=:before)
    )
    @test_throws ArgumentError Formulation(
        insulation_impedance=formula(:Lossless; order=:before)
    )
    @test_throws ArgumentError Formulation(
        insulation_admittance=formula(:Lossless; order=:after)
    )
    @test_throws ArgumentError Formulation(
        equivalent_earth=formula(:Layer; unknown=true)
    )
end

@testitem "Engine / external earth-impedance vocabulary" tags=[:unit] setup=[
    UseEngineSupport,
    TestFixtures
] begin
    earth_impedance=LineCableModels.Engine.EarthImpedance
    formulations=(
        earth_impedance.Formula(:Deri)=>"Deri-Semlyen",
        earth_impedance.Formula(:Wedepohl)=>"Wedepohl",
        earth_impedance.Formula(:Saad)=>"Saad",
        earth_impedance.Formula(:Ametani)=>"Ametani",
        earth_impedance.Formula(:Lucca)=>"Lucca"
    )
    for (formulation, label) in formulations
        @test supertype(typeof(formulation)) ===
              LineCableModels.Engine.EarthImpedanceFormulation
        @test description(formulation) == label
    end
    @test !isdefined(earth_impedance, :ReferenceEarthImpedance)
    @test !isdefined(earth_impedance, :DirectNumericalIntegration)

    problem=TestFixtures.line_parameters_problem()
    formulation=Formulation(earth_impedance = :Saad)
    @test_throws MethodError compute(problem, formulation)
end

@testitem "Engine / insulation formulations / analytical limits across precision" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics
] begin
    using TOML

    impedance_formulation=InsulationImpedance.Formula(:Lossless)
    admittance_formulation=InsulationAdmittance.Formula(:Lossless)
    @test description(impedance_formulation) == "Lossless insulation (ideal dielectric)"
    @test description(admittance_formulation) ==
          "Lossless insulation (ideal dielectric)"
    @test formula_id(impedance_formulation) === :Lossless
    @test formula_id(admittance_formulation) === :Lossless
    @test InsulationImpedance.formulas() == (:Lossless,)
    @test InsulationAdmittance.formulas() == (:Lossless, :ParallelRC)
    @test !isdefined(InsulationImpedance, :Lossless)
    @test !isdefined(InsulationAdmittance, :Lossless)
    @test !isdefined(InsulationAdmittance, :ParallelRC)
    reference=TOML.parsefile(joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "coaxial_capacitance.toml"
    ))["lossless_formulation"]

    for T in (Float32, Float64, BigFloat)
        setprecision(BigFloat, 128) do
            r_in=parse(T, "0.01")
            r_ex=parse(T, "0.02")
            relative_permeability=parse(T, "1.2")
            relative_permittivity=parse(T, "2.3")
            angular_frequency=T(2)*T(π)*T(50)
            s=Complex{T}(zero(T), angular_frequency)

            impedance=impedance_formulation(
                r_in,
                r_ex,
                relative_permeability,
                s
            )
            potential=admittance_formulation(
                r_in,
                r_ex,
                T(Inf),
                relative_permittivity,
                s
            )
            @test impedance isa Complex{T}
            @test potential isa Complex{T}
            @test TestNumerics.isapprox_scaled(
                impedance,
                Complex{T}(
                    zero(T),
                    T(parse(BigFloat, reference["series_impedance_imaginary"]))
                )
            )
            @test TestNumerics.isapprox_scaled(
                potential,
                Complex{T}(T(parse(BigFloat, reference["potential_coefficient"])))
            )
            @test iszero(impedance_formulation(zero(T), r_ex, one(T), s))
            @test iszero(impedance_formulation(r_ex, r_ex, one(T), s))
            @test iszero(admittance_formulation(
                zero(T), r_ex, T(Inf), one(T), s
            ))
            @test iszero(admittance_formulation(
                r_ex, r_ex, T(Inf), one(T), s
            ))
        end
    end

    impedance_route=(r_in, r_ex, mu_r, s, values)->values.scale*s
    admittance_route=(r_in, r_ex, rho, eps_r, s, values)->
        Complex(values.scale)
    experimental_impedance=InsulationImpedance.Formula(
        :Experiment,
        impedance_route,
        (scale = 2.0,)
    )
    experimental_admittance=InsulationAdmittance.Formula(
        :Experiment,
        admittance_route,
        (scale = 3.0,)
    )
    @test @inferred(experimental_impedance(0.01, 0.02, 1.0, 2.0im)) == 4.0im
    @test @inferred(experimental_admittance(
        0.01, 0.02, 1.0e12, 2.3, 2.0im
    )) == 3.0+0.0im
    @test_throws ArgumentError InsulationImpedance.Formula(:Lossless; bad=true)
    @test_throws ArgumentError InsulationAdmittance.Formula(:ParallelRC; bad=true)
end

@testitem "Engine / internal impedance / passivity and solid-conductor limits" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    formulation=InternalImpedance.Formula(:Schelkunoff)
    @test description(formulation) == "Schelkunoff"
    @test InternalImpedance.formula_id(formulation) === :Schelkunoff
    @test InternalImpedance.formulas() == (:Schelkunoff,)

    r_in=0.005
    r_ex=0.01
    rho=1.7241e-8
    relative_permeability=1.0
    s=ComplexF64(0.0, 2π*50.0)
    interaction=@inferred formulation(
        r_in, r_ex, rho, relative_permeability, s)
    inner=@inferred interaction(Val(:inner))
    outer=@inferred interaction(Val(:outer))
    mutual=@inferred interaction(Val(:mutual))
    @test all(isfinite, (inner, outer, mutual))
    @test real(inner) > 0
    @test real(outer) > 0
    @test real(mutual) > 0
    @test imag(inner) >= 0
    @test imag(outer) >= 0

    solid=@inferred formulation(
        0.0, r_ex, rho, relative_permeability, s)
    @test iszero(solid(Val(:inner)))
    @test iszero(solid(Val(:mutual)))
    solid_outer=solid(Val(:outer))
    @test isfinite(solid_outer)
    @test real(solid_outer) > 0
    @test_throws ArgumentError interaction(:unsupported)

    custom_inner=state->oftype(state.jω, 7)
    experiment=InternalImpedance.Formula(
        :Schelkunoff; inner = custom_inner)
    experimental=@inferred experiment(
        r_in, r_ex, rho, relative_permeability, s)
    @test experimental(Val(:inner)) == 7
    @test experimental(Val(:outer)) == outer
    @test experimental(Val(:mutual)) == mutual
end

@testitem "Engine / homogeneous earth / reciprocity and boundary validation" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    using QuadGK: quadgk

    impedance_module=LineCableModels.Engine.EarthImpedance
    admittance_module=LineCableModels.Engine.EarthAdmittance
    epsilon0=8.8541878128e-12
    mu0=4π*1.0e-7
    rho=[Inf, 100.0]
    epsilon=[epsilon0, 10*epsilon0]
    mu=[mu0, mu0]
    s=ComplexF64(0.0, 2π*50.0)

    impedance_cases=(
        (impedance_module.Formula(:Papadopoulos2010), (-1.0, -1.2), (2, 2)),
        (impedance_module.Formula(:Pollaczek), (-1.0, -1.2), (2, 2)),
        (impedance_module.Formula(:Carson), (10.0, 12.0), (1, 1))
    )
    @test description(first(impedance_cases)[1]) == "Papadopoulos"
    @test description(impedance_cases[2][1]) == "Pollaczek"
    @test description(impedance_cases[3][1]) == "Carson"
    @test impedance_module.formula_id(first(impedance_cases)[1]) ===
          :Papadopoulos2010
    @test LineCableModels.Engine.formula_id(first(impedance_cases)[1]) ===
          :Papadopoulos2010
    @test keys(impedance_module.assumptions(first(impedance_cases)[1])) == (
        :air,
        :earth,
        :permeability
    )
    @test Formulation(earth_impedance = :Papadopoulos2010).methods.earth_impedance ==
          first(impedance_cases)[1]
    @test isconcretetype(typeof(
        Formulation(earth_impedance = :Papadopoulos2010).methods.earth_impedance
    ))
    @test impedance_module.formulas() == (
        :Ametani, :Carson, :Deri, :Lucca, :Papadopoulos2010,
        :Pollaczek, :Saad, :Wedepohl)
    for (formulation, heights, layers) in impedance_cases
        functor=@inferred formulation(rho, epsilon, mu, s, nothing)
        pair=EarthPair(1, 2, heights, 0.25, layers)
        reciprocal_pair=EarthPair(2, 1, reverse(heights), 0.25, reverse(layers))
        mutual=@inferred functor(Val(:mutual), pair)
        reciprocal=@inferred functor(Val(:mutual), reciprocal_pair)
        self=@inferred functor(Val(:self), pair)
        @test isfinite(mutual)
        @test mutual ≈ reciprocal
        @test self == mutual
    end

    let regression_formula=impedance_module.Formula(:Papadopoulos2010),
        regression_rho=Float64[Inf, 100.0],
        regression_epsilon=Float64[epsilon0, 10 * epsilon0],
        regression_mu=Float64[mu0, mu0], regression_frequency=ComplexF64(0.0, 2π*50.0)

        regression_functor=@inferred regression_formula(
            regression_rho,
            regression_epsilon,
            regression_mu,
            regression_frequency,
            nothing
        )
        regression_integrand=impedance_module.Integrand(
            regression_functor, 2.2, 0.25
        )
        finite_part, _=quadgk(
            regression_integrand, 0.0, 1.0; rtol = 1.0e-10
        )
        infinite_tail, _=quadgk(
            regression_integrand, 1.0, Inf; rtol = 1.0e-10
        )
        @test impedance_module.formula_id(regression_formula) ===
              :Papadopoulos2010
        @test impedance_module.propagation(regression_formula) === Val(:explicit)
        @test impedance_module._integral(regression_functor, 2.2, 0.25) ≈
              2 * (finite_part + infinite_tail) rtol=1.0e-10
        @test abs(infinite_tail) > 0.005abs(finite_part + infinite_tail)
        pair=EarthPair(1, 2, (-1.0, -1.2), 0.25, (2, 2))
        @test regression_functor(Val(:mutual), pair) == ComplexF64(
            4.947519285382972e-5,
            5.010204379584122e-4
        )

        explicit=impedance_module.Γ(regression_functor)
        explicit_functor=@inferred regression_formula(
            regression_rho,
            regression_epsilon,
            regression_mu,
            regression_frequency,
            explicit
        )
        @test impedance_module.Γ(explicit_functor) === explicit
        @test explicit_functor(Val(:mutual), pair) ≈
              regression_functor(Val(:mutual), pair) rtol=1.0e-12
    end

    admittance_cases=(
        (admittance_module.Formula(:Papadopoulos2010), (-1.0, -1.2), (2, 2)),
        (admittance_module.Formula(:Pollaczek), (-1.0, -1.2), (2, 2)),
        (admittance_module.Formula(:ElectrostaticImages), (10.0, 12.0), (1, 1))
    )
    @test description(first(admittance_cases)[1]) == "Papadopoulos"
    @test description(admittance_cases[2][1]) == "Pollaczek"
    @test description(admittance_cases[3][1]) == "Electrostatic images"
    @test admittance_module.formulas() == (
        :ElectrostaticImages, :IdealGround, :Papadopoulos2010, :Pollaczek)
    for (formulation, heights, layers) in admittance_cases
        functor=@inferred formulation(rho, epsilon, mu, s, nothing)
        pair=EarthPair(1, 2, heights, 0.25, layers)
        reciprocal_pair=EarthPair(2, 1, reverse(heights), 0.25, reverse(layers))
        mutual=@inferred functor(Val(:mutual), pair)
        reciprocal=@inferred functor(Val(:mutual), reciprocal_pair)
        @test isfinite(mutual)
        @test mutual ≈ reciprocal
        @test functor(Val(:self), pair) == mutual
    end

    ideal_ground=admittance_module.Formula(:IdealGround)
    @test description(ideal_ground) == "Ideal ground reference"
    ideal_functor=@inferred ideal_ground(rho, epsilon, mu, s, nothing)
    underground_pair=EarthPair(1, 2, (-1.0, -1.2), 0.25, (2, 2))
    @test iszero(ideal_functor(Val(:self), underground_pair))
    @test iszero(ideal_functor(Val(:mutual), underground_pair))

    underground_impedance=impedance_module.Formula(:Papadopoulos2010)
    underground_admittance=admittance_module.Formula(:Papadopoulos2010)
    for formulation in (underground_impedance, underground_admittance)
        @test_throws DimensionMismatch formulation(
            rho, epsilon[1:1], mu, s, nothing)
        functor=formulation(rho, epsilon, mu, s, nothing)
        wrong_layer=EarthPair(1, 2, (1.0, 1.2), 0.25, (1, 1))
        @test_throws ArgumentError functor(Val(:mutual), wrong_layer)
    end
    @test_throws ArgumentError impedance_module.Formula(:Papadopoulos2010; bad = +)
    @test_throws ArgumentError admittance_module.Formula(:Papadopoulos2010; bad = +)
    @test_throws ArgumentError impedance_module.Formula(:Pollaczek)(
        rho, epsilon, mu, s, ComplexF64(1))
    @test_throws ArgumentError admittance_module.Formula(:Pollaczek)(
        rho, epsilon, mu, s, ComplexF64(1))
end
