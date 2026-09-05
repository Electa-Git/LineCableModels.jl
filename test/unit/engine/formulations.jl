@testitem "Engine / formula selector / owner-local lowering" tags=[:unit] setup=[
    UseEngineSupport,
] begin
    const EP=LineCableModels.EarthProps
    const EH=EP.EHEM
    const EN=LineCableModels.Engine

    inner=state->oftype(state.jω, 7)
    insulation_impedance_route=(r_in, r_ex, mu_r, s, values)->zero(s)
    insulation_admittance_route=(material, frequency, temperature, values)->complex(
        inv(material.rho), frequency*material.eps_r)
    semicon_admittance_route=(material, frequency, temperature, values)->complex(
        inv(material.rho), frequency*material.eps_r)
    formulation=@inferred Formulation(
        internal_impedance = formula(:Schelkunoff1934; inner),
        insulation_impedance = formula(
            :Ametani1980;
            route = insulation_impedance_route
        ),
        earth_impedance = formula(:Papadopoulos2010),
        insulation_admittance = formula(
            :Ametani2004;
            route = insulation_admittance_route
        ),
        semicon_admittance = formula(
            :Ametani2004;
            route = semicon_admittance_route
        ),
        earth_admittance = formula(:Papadopoulos2010),
        earth_properties = formula(:CIGRE2019; epsilon_infinity = 10.0),
        equivalent_earth = formula(:Xue2021; order = :before, layer = 2)
    )

    @test formula_id(formulation.methods.internal_impedance) === :Schelkunoff1934
    @test EN.InternalImpedance.routes(
        formulation.methods.internal_impedance
    ).inner === inner
    @test formula_id(formulation.methods.insulation_impedance) === :Ametani1980
    @test formulation.methods.insulation_impedance.route ===
          insulation_impedance_route
    @test formula_id(formulation.methods.earth_impedance) === :Papadopoulos2010
    @test formula_id(formulation.methods.insulation_admittance) === :Ametani2004
    @test formulation.methods.insulation_admittance.route ===
          insulation_admittance_route
    @test formula_id(formulation.methods.semicon_admittance) === :Ametani2004
    @test formulation.methods.semicon_admittance.route === semicon_admittance_route
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
    @test defaults.methods.earth_properties === nothing
    @test defaults.methods.equivalent_earth isa EH.AfterFD
    @test EH.rule(defaults.methods.equivalent_earth) == EH.Layer(-1)
    @test Formulation(earth_properties = :default).methods.earth_properties ===
          nothing
    explicit_default=Formulation(equivalent_earth = :default)
    @test explicit_default.methods.equivalent_earth isa EH.AfterFD
    @test EH.rule(explicit_default.methods.equivalent_earth) == EH.Layer(-1)
    bare=Formulation(
        insulation_impedance = :Ametani1980,
        insulation_admittance = :Gustavsen2013,
        semicon_admittance = :Ametani2004
    )
    @test formula_id(defaults.methods.insulation_impedance) === :Ametani1980
    @test formula_id(defaults.methods.insulation_admittance) === :Ametani2004
    @test formula_id(defaults.methods.semicon_admittance) === :Ametani2004
    for (owner, identifier) in (
        EN.InternalImpedance=>:Schelkunoff1934,
        EN.InsulationImpedance=>:Ametani1980,
        EN.EarthImpedance=>:Papadopoulos2010,
        EN.InsulationAdmittance=>:Ametani2004,
        EN.SemiconAdmittance=>:Ametani2004,
        EN.EarthAdmittance=>:Papadopoulos2010
    )
        @test owner.DEFAULT === identifier
        @test formula_id(owner.Formula(:default)) === identifier
        @test :default ∉ owner.formulas()
    end
    @test_throws ArgumentError EN.InsulationAdmittance.Formula(:Marti2001)
    @test formula_id(bare.methods.insulation_impedance) === :Ametani1980
    @test formula_id(bare.methods.insulation_admittance) === :Gustavsen2013
    @test formula_id(bare.methods.semicon_admittance) === :Ametani2004

    constants_formulation=@inferred CableConstantsFormulation()
    @test keys(constants_formulation.methods) == (
        :internal_impedance,
        :insulation_impedance,
        :insulation_admittance,
        :semicon_admittance
    )
    @test formula_id(constants_formulation.methods.internal_impedance) ===
          :Schelkunoff1934
    @test formula_id(constants_formulation.methods.insulation_impedance) ===
          :Ametani1980
    @test formula_id(constants_formulation.methods.insulation_admittance) ===
          :Ametani2004
    @test formula_id(constants_formulation.methods.semicon_admittance) ===
          :Ametani2004
    @test constants_formulation.options == (temperature_correction = true,)
    @test CableConstantsFormulation(
        options = (temperature_correction = false,)
    ).options.temperature_correction === false
    @test_throws ArgumentError CableConstantsFormulation(
        options = (kron_reduction = true,)
    )

    after_default=@inferred Formulation(equivalent_earth = formula(:Xue2021))
    after_explicit=@inferred Formulation(
        equivalent_earth = formula(:Xue2021; order = :after)
    )
    layer=@inferred Formulation(
        equivalent_earth = formula(:Layer; order = :before, layer = 2)
    )
    @test after_default.methods.equivalent_earth isa EH.AfterFD
    @test after_explicit.methods.equivalent_earth isa EH.AfterFD
    @test layer.methods.equivalent_earth isa EH.BeforeFD
    @test EH.rule(layer.methods.equivalent_earth) == EH.Layer(2)

    selection=formula(:Xue2021; order = :before)
    @test formula_id(selection) === :Xue2021
    @test isconcretetype(typeof(selection))
    @test_throws ArgumentError formula(:Xue2021; order = :sideways)
    @test_throws ArgumentError Formulation(
        earth_properties = formula(:CIGRE2019; order = :before)
    )
    @test_throws ArgumentError Formulation(
        insulation_impedance = formula(:Ametani1980; order = :before)
    )
    @test_throws ArgumentError Formulation(
        insulation_admittance = formula(:Ametani2004; order = :after)
    )
    @test_throws ArgumentError Formulation(
        semicon_admittance = formula(:Ametani2004; order = :after)
    )
    @test_throws ArgumentError Formulation(
        equivalent_earth = formula(:Layer; unknown = true)
    )
    @test_throws ArgumentError Formulation(equivalent_earth = formula(:Layer))
end

@testitem "Engine / external earth-impedance vocabulary" tags=[:unit] setup=[
    UseEngineSupport,
    TestFixtures
] begin
    earth_impedance=LineCableModels.Engine.EarthImpedance
    formulations=(
        earth_impedance.Formula(:DeriSemlyen1981)=>
            "Deri-Semlyen complex ground-return-plane approximation (1981)",
        earth_impedance.Formula(:Ametani2009)=>
            "Ametani pair-complete homogeneous-earth impedance (2009)",
        earth_impedance.Formula(:Lucca1994)=>
            "Lucca pair-complete homogeneous-earth impedance (1994)"
    )
    for (formulation, label) in formulations
        @test supertype(typeof(formulation)) ===
              LineCableModels.Engine.EarthImpedanceFormulation
        @test description(formulation) == label
    end
    @test !isdefined(earth_impedance, :ReferenceEarthImpedance)
    @test !isdefined(earth_impedance, :DirectNumericalIntegration)
end

@testitem "Engine / insulation formulations / analytical limits across precision" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestNumerics
] begin
    using TOML

    impedance_formulation=InsulationImpedance.Formula(:Ametani1980)
    admittance_formulation=InsulationAdmittance.Formula(:Gustavsen2013)
    @test description(impedance_formulation) ==
          "Ametani coaxial-insulation magnetic impedance (1980)"
    @test occursin("lossless", lowercase(description(admittance_formulation)))
    @test formula_id(impedance_formulation) === :Ametani1980
    @test formula_id(admittance_formulation) === :Gustavsen2013
    @test InsulationImpedance.formulas() == (:Ametani1980,)
    @test InsulationAdmittance.formulas() == (:Ametani2004, :Gustavsen2013)
    @test SemiconAdmittance.formulas() == (:Ametani2004,)
    @test !isdefined(InsulationAdmittance, :Gustavsen2013)
    @test !isdefined(InsulationAdmittance, :Ametani2004)
    @test !isdefined(SemiconAdmittance, :Ametani2004)
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
            material=Material(
                :insulator, T(1.0e12), relative_permittivity, one(T),
                T(20), zero(T)
            )
            evaluated=admittance_formulation(material, T(50), T(20))
            potential=LineCableModels.Engine.potential_coefficient(
                r_in, r_ex, evaluated, s
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
            @test iszero(LineCableModels.Engine.potential_coefficient(
                zero(T), r_ex, evaluated, s
            ))
            @test iszero(LineCableModels.Engine.potential_coefficient(
                r_ex, r_ex, evaluated, s
            ))
        end
    end

    impedance_route=(r_in, r_ex, mu_r, s, values)->values.scale*s
    admittance_route=(material, frequency, temperature, values)->complex(
        inv(material.rho), values.scale*material.eps_r)
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
    material=Material(:insulator, 1.0e12, 2.3, 1.0, 20.0, 0.0)
    @test imag(@inferred(experimental_admittance(material, 50.0, 20.0))) ≈ 6.9
    @test_throws ArgumentError InsulationImpedance.Formula(:Ametani1980; bad = true)
    @test_throws ArgumentError InsulationAdmittance.Formula(:Ametani2004; bad = true)
end

@testitem "Engine / internal impedance / passivity and solid-conductor limits" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    formulation=InternalImpedance.Formula(:Schelkunoff1934)
    @test occursin("Schelkunoff", description(formulation))
    @test InternalImpedance.formula_id(formulation) === :Schelkunoff1934
    @test :Schelkunoff1934 in InternalImpedance.formulas()

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
        :Schelkunoff1934; inner = custom_inner)
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
        (impedance_module.Formula(:Pollaczek1926), (-1.0, -1.2), (2, 2)),
        (impedance_module.Formula(:Pollaczek1926), (10.0, 12.0), (1, 1))
    )
    @test occursin("Papadopoulos", description(first(impedance_cases)[1]))
    @test occursin("Pollaczek", description(impedance_cases[2][1]))
    @test occursin("overhead", description(impedance_cases[3][1]))
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
    @test all(in(impedance_module.formulas()),
        (
            :Papadopoulos2010,
            :Pollaczek1926,
            :Theodoulidis2015,
            :Tsiamitros2008
        ))
    for (formulation, heights, layers) in impedance_cases
        functor=formulation(rho, epsilon, mu, s, nothing)
        pair=EarthPair(1, 2, heights, 0.25, layers)
        reciprocal_pair=EarthPair(2, 1, reverse(heights), 0.25, reverse(layers))
        self_pair=EarthPair(1, 1, (heights[1], heights[1]), 0.01, layers)
        mutual=functor(Val(:mutual), pair)
        reciprocal=functor(Val(:mutual), reciprocal_pair)
        self=functor(Val(:self), self_pair)
        @test isfinite(mutual)
        @test isfinite(self)
        @test mutual ≈ reciprocal
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
        @test impedance_module.formula_id(regression_formula) ===
              :Papadopoulos2010
        @test impedance_module.propagation(regression_formula) === Val(:explicit)
        pair=EarthPair(1, 2, (-1.0, -1.2), 0.25, (2, 2))
        reference=regression_functor(Val(:mutual), pair)
        @test isfinite(reference)

        explicit=impedance_module.Γ(regression_functor)
        explicit_functor=@inferred regression_formula(
            regression_rho,
            regression_epsilon,
            regression_mu,
            regression_frequency,
            explicit
        )
        @test impedance_module.Γ(explicit_functor) === explicit
        @test explicit_functor(Val(:mutual), pair) ≈ reference rtol=1.0e-12
    end

    admittance_cases=(
        (admittance_module.Formula(:Papadopoulos2010), (-1.0, -1.2), (2, 2)),
        (admittance_module.Formula(:Xue2021), (-1.0, -1.2), (2, 2)),
        (admittance_module.Formula(:Ametani2021), (10.0, 12.0), (1, 1))
    )
    @test occursin("Papadopoulos", description(first(admittance_cases)[1]))
    @test occursin("Xue", description(admittance_cases[2][1]))
    @test occursin("space", description(admittance_cases[3][1]))
    @test all(in(admittance_module.formulas()),
        (
            :Ametani2021,
            :Papadopoulos2010,
            :Theethayi2007,
            :Xue2021
        ))
    for (formulation, heights, layers) in admittance_cases
        functor=formulation(rho, epsilon, mu, s, nothing)
        pair=EarthPair(1, 2, heights, 0.25, layers)
        reciprocal_pair=EarthPair(2, 1, reverse(heights), 0.25, reverse(layers))
        self_pair=EarthPair(1, 1, (heights[1], heights[1]), 0.01, layers)
        mutual=functor(Val(:mutual), pair)
        reciprocal=functor(Val(:mutual), reciprocal_pair)
        @test isfinite(mutual)
        @test isfinite(functor(Val(:self), self_pair))
        @test mutual ≈ reciprocal
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
    @test_throws ArgumentError impedance_module.Formula(:Pollaczek1926)(
        rho, epsilon, mu, s, ComplexF64(1))
    @test_throws ArgumentError admittance_module.Formula(:Xue2021)(
        rho, epsilon, mu, s, ComplexF64(1))
end
