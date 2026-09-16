@testitem "Engine / insulation formulations / analytical limits across precision" tags=[:unit] setup=[
    UseEngineSupport,
    TestNumerics,
    FormulaContractModels
] begin

    impedance_formulation=InsulationImpedance.Formula(:default)
    admittance_formulation=InsulationAdmittance.Formula(:default)
    @test description(impedance_formulation) ==
          "Ametani coaxial-insulation magnetic impedance (1980)"
    @test occursin("lossless", lowercase(description(admittance_formulation)))
    @test formula_id(impedance_formulation) === :ametani1980
    @test formula_id(admittance_formulation) === :lossless
    @test :default in InsulationImpedance.formulas()
    @test all(in(InsulationAdmittance.formulas()), (:lossless, :lossy, :default))
    @test all(in(SemiconAdmittance.formulas()), (:lossless, :lossy, :default))
    # Current lossless radial fields: H=I/(2pi*r) and
    # E=V/(r*log(b/a)); integrate their magnetic/electric energy.
    for T in (Float32, Float64, BigFloat)
        setprecision(BigFloat, 128) do
            r_in=parse(T, "0.005")
            r_ex=parse(T, "0.01")
            relative_permeability=one(T)
            relative_permittivity=T(3)
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
                    angular_frequency*T(4)*T(pi)*parse(T,"1e-7")*relative_permeability*log(r_ex/r_in)/(T(2)*T(pi))
                )
            )
            @test TestNumerics.isapprox_scaled(
                potential,
                Complex{T}(log(r_ex/r_in)/(T(2)*T(pi)*parse(T,"8.8541878128e-12")*relative_permittivity))
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

    experimental_impedance=FormulaContractModels.InsulationReactance()
    experimental_admittance=FormulaContractModels.InsulationLaw()
    @test @inferred(experimental_impedance(0.01, 0.02, 1.0, 2.0im)) == 4.0im
    @test_throws DomainError experimental_impedance(-0.01, 0.02, 1.0, 2.0im)
    @test_throws DomainError experimental_impedance(0.01, 0.02, -1.0, 2.0im)
    material=Material(:insulator, 1.0e12, 2.3, 1.0, 20.0, 0.0)
    @test imag(@inferred(experimental_admittance(material, 50.0, 20.0))) ≈ 135.0
    @test_throws ArgumentError InsulationImpedance.Formula(:default; parameters = (bad = true,))
    @test_throws ArgumentError InsulationAdmittance.Formula(:lossy; parameters = (bad = true,))
end

@testitem "Engine / internal impedance / passivity and solid-conductor limits" tags=[:unit] setup=[
    UseEngineSupport,FormulaContractModels
] begin
    formulation=InternalImpedance.Formula(:default)
    @test occursin("Schelkunoff", description(formulation))
    @test InternalImpedance.formula_id(formulation) === :schelkunoff1934
    @test :default in InternalImpedance.formulas()

    r_in=0.005
    r_ex=0.01
    rho=1.7241e-8
    relative_permeability=1.0
    s=ComplexF64(0.0, 2π*50.0)
    interaction=@inferred formulation(
        r_in, r_ex, rho, relative_permeability, s)
    inner=@inferred interaction(Val(:inner))
    outer=@inferred interaction(Val(:outer))
    transfer=@inferred interaction(Val(:transfer))
    @test all(isfinite, (inner, outer, transfer))
    @test real(inner) > 0
    @test real(outer) > 0
    @test real(transfer) > 0
    @test imag(inner) >= 0
    @test imag(outer) >= 0

    solid=@inferred formulation(
        0.0, r_ex, rho, relative_permeability, s)
    @test iszero(solid(Val(:inner)))
    @test iszero(solid(Val(:transfer)))
    solid_outer=solid(Val(:outer))
    @test isfinite(solid_outer)
    @test real(solid_outer) > 0
    @test_throws ArgumentError interaction(:unsupported)

    custom_inner=FormulaContractModels.SurfaceLaw(kinds=(:inner,),coefficients=(inner=7+0im,))
    experiment=(inner=custom_inner,outer=formulation,transfer=formulation)
    experimental=InternalImpedance.surface_impedances(experiment,
        r_in,r_ex,rho,relative_permeability,s)
    @test experimental.inner==7
    @test experimental.outer==outer
    @test experimental.transfer==transfer

end
