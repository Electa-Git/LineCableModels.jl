@testitem "Engine / numerical declarations follow selected types and indexed equations" tags=[:unit] setup=[FormulaContractModels] begin
    const E=LineCableModels.Engine
    const II=E.InternalImpedance
    const EI=E.EarthImpedance
    const FM=LineCableModels.FormulaMethod
    const M=FormulaContractModels
    internal=FM(II.Formula(:default),II.internal_impedance,Val(:outer))
    external=FM(EI.Formula(:default),EI.earth_impedance,Val(:self),Val(1),Val(1))
    @test formulation_options(internal)==FormulationOptions()
    @test formulation_options(external).data.integration.method === :quad
    @test_throws ArgumentError formulation_options(internal,FormulationOptions(integration=(method=:quad,)))
    @test_throws ArgumentError computation_options(LineCableModelsCoaxial, ComputationOptions((integration_method=:quad,)))
    # The custom selection uses existing admission and equation generics. It is
    # not a changed implementation of the built-in's claimed scientific identity.
    custom=M.selection(EI;layers=2:2)
    pair=E.EarthPair(1,2,(-1.0,-1.0),1.0,(2,2))
    bound=validate(custom,pair)
    @test bound.equation.selection === custom
    @test isempty(bound.options.data)
    @test_throws ArgumentError validate(M.selection(EI;layers=2:2,
        options=(integration=(method=:quad,),)),pair)
    rho=[Inf,100.0]
    epsilon=8.8541878128e-12 .* [1,10]
    mu=fill(4pi*1e-7,2)
    value=custom(rho,epsilon,mu,100.0im,pair;thickness=[Inf,Inf])()
    @test isfinite(value)
    @test formula_id(custom) !== formula_id(external.selection)
    @test formulation_options(external).data.integration.method === :quad
end

@testitem "Engine / public internal surfaces consume spectral options and retained resources" tags=[:unit] setup=[TestFixtures,FormulaContractModels] begin
    using QuadGK
    const E=LineCableModels.Engine
    const II=E.InternalImpedance
    const M=FormulaContractModels
    @test_throws ArgumentError II.Formula(:default;options=(integration=(method=:quad,),))
    resources=(segments=alloc_segbuf(Float64,ComplexF64,Float64;size=32),
        images=ComplexF64[],exponents=ComplexF64[])
    base=II.Formula(:default)
    args=(0.008,0.01,1.7241e-8,1.0,100.0im)
    reference=II.surface_impedances(base,args...)
    for method in (:quad,:trapz,:cim)
        selected=M.SpectralSurface(method)
        surfaces=(inner=base,outer=selected,transfer=base)
        values=II.surface_impedances(surfaces,args...;workspace=resources)
        @test values.outer ≈ 5e-5 rtol=3e-6
        @test values.inner == reference.inner
        @test values.transfer == reference.transfer
        @test last(selected.seen) === (Val(method),resources)
    end
    problem=TestFixtures.line_parameters_problem(frequencies=[50.0])
    results=map((:quad,:trapz,:cim)) do method
        selected=M.SpectralSurface(method)
        result=compute(problem,Formulation(internal_impedance=
            (inner=base,outer=selected,transfer=base)))
        @test !isempty(selected.seen)
        @test all(record -> record[1] === Val(method),selected.seen)
        result
    end
    @test results[2].Z.values ≈ results[1].Z.values rtol=3e-6
    @test results[3].Z.values ≈ results[1].Z.values rtol=3e-6
end

@testitem "Engine / surface current basis reproduces concentric wall contributions" tags=[:unit] begin
    using LinearAlgebra
    const II = LineCableModels.Engine.InternalImpedance
    selected = II.Formula(:default)
    for frequency in (1e-4, 50.0, 1e5)
        coefficients = II.surface_impedances(selected, 0.008, 0.01, 1.7241e-8, 1.0,
            complex(0.0, 2pi*frequency))
        W = [coefficients.inner coefficients.transfer; coefficients.transfer coefficients.outer]
        # Terminal currents are (contained metal, enclosing wall). Surface
        # currents are (-contained, contained + wall), exactly the current map.
        B = [-1.0 1.0; 0.0 1.0]
        lifted = B * W * transpose(B)
        @test lifted[1, 1] ≈ coefficients.inner - 2coefficients.transfer + coefficients.outer
        @test lifted[1, 2] ≈ coefficients.outer - coefficients.transfer
        @test lifted[2, 2] == coefficients.outer
        if frequency == 1e-4
            resistance = 1.7241e-8 / (pi * (0.01^2 - 0.008^2))
            @test real(lifted[2, 2]) ≈ resistance rtol=1e-10
            @test abs(real(lifted[1, 1])) < 1e-10resistance
            @test abs(real(lifted[1, 2])) < 1e-10resistance
        end
    end
end

@testitem "Earth / artificial material values are distinct from source restrictions" tags=[:unit] setup=[FormulaContractModels] begin
    const M=FormulaContractModels
    material=M.EP.EarthMaterial(100.0, -10.0, 1.0)
    @test material.eps_r == -10
    @test M.EP.EarthLayer <: M.EP.AbstractEarthLayer <: M.EP.AbstractEarthModel
    @test M.EP.EarthMaterial <: M.EP.AbstractEarthMaterial <:
          LineCableModels.AbstractMaterial
    rho=[Inf, 100.0, 200.0]
    epsilon=8.8541878128e-12 .* [1, -10, -20]
    mu=fill(4pi*1e-7, 3)
    pair=M.E.EarthPair(1, 2, (-0.25, -1.5), 1.0, (2, 3))
    @test isfinite(M.selection(M.EI)(rho, epsilon, mu, 100.0im, pair;
        thickness = [Inf, 0.5, Inf])())
    @test_throws DomainError M.EI.Formula(:default)(
        rho[1:2], epsilon[1:2], mu[1:2], 100.0im,
        M.E.EarthPair(1, 2, (-0.25, -1.5), 1.0, (2, 2)))
end

@testitem "Engine / internal consumers request only their actual surface kinds" tags=[:unit] setup=[FormulaContractModels] begin
    const II=LineCableModels.Engine.InternalImpedance
    selected=FormulaContractModels.SurfaceLaw(kinds=(:outer,))
    @test validate(selected,(:outer,)) === selected
    @test_throws ArgumentError validate(selected,(:inner,:outer,:transfer))
    @test_throws ArgumentError II.surface_impedances(selected,0.0,0.01,1.7e-8,1.0,100im)
    copper=Material(kind=:conductor,rho=1.7e-8)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3)
    design=build(CableDesign,"outer-only",
        terminal(:core,solid(copper,Disk(0.01)),insulation(dielectric;t=0.002)))
    system=build(LineCableSystem,design,Pose2(0.0,5.0);connections=Dict(:core=>1))
    problem=LineParametersProblem(system;frequencies=[50.0],earth_props=homogeneous(rho=100.0))
    result=compute(problem,Formulation(internal_impedance=selected))
    @test all(isfinite,result.Z.values)
    @test keys(details(result).data.formulations.numerical.internal_impedance)==(:outer,)
    @test length(selected.preparations)==1
    @test only(selected.evaluations)[2] === Val(:outer)
end

@testitem "Engine / each metal prepares shared surface state once per frequency" tags=[:unit] setup=[TestFixtures,FormulaContractModels] begin
    selected=FormulaContractModels.SurfaceLaw()
    problem=TestFixtures.line_parameters_problem(frequencies=[50.0,100.0])
    result=compute(problem,Formulation(internal_impedance=selected))
    expected=sum(length(design.terminal_order) for design in problem.system.designs)*length(problem.frequencies)
    @test length(selected.preparations)==expected
    @test allunique(selected.evaluations)
    @test count(record -> record[2] === Val(:outer),selected.evaluations)==expected
    @test any(record -> record[2] === Val(:inner),selected.evaluations)
    @test any(record -> record[2] === Val(:transfer),selected.evaluations)
    @test all(isfinite,result.Z.values)
    empty!(selected.preparations)
    empty!(selected.evaluations)
    composite=compute(problem,Formulation(internal_impedance=
        (inner=selected,outer=selected,transfer=selected)))
    @test length(selected.preparations)==expected
    @test allunique(selected.evaluations)
    @test Z(composite)==Z(result) && Y(composite)==Y(result)
end
