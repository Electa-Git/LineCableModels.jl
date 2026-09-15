@testitem "Engine / indexed formula declarations and immutable source domains" tags=[:unit] begin
    const E = LineCableModels.Engine
    const EI = E.EarthImpedance
    const EA = E.EarthAdmittance
    const FM = LineCableModels.FormulaMethod
    air = E.EarthPair(1, 2, (10.0, 12.0), 1.0, (1, 1))
    soil = E.EarthPair(1, 2, (-1.0, -2.0), 1.0, (2, 2))
    mixed = E.EarthPair(1, 2, (10.0, -2.0), 1.0, (1, 2))
    self = E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2); radius = 0.01)
    @test validate(self) === self
    @test self.radius == 0.01 && iszero(self.separation)
    @test_throws DomainError validate(E.EarthPair(1, 1, (-1.0, -1.0), 0.0, (2, 2)))
    @test_throws ArgumentError validate(E.EarthPair(
        1, 1, (-1.0, -1.0), 0.01, (2, 2); radius = 0.01))
    @test_throws ArgumentError validate(E.EarthPair(
        1, 1, (-1.0, -2.0), 0.0, (2, 2); radius = 0.01))
    @test_throws ArgumentError validate(E.EarthPair(1, 2, (1.0, 2.0), 1.0, (2, 2)))
    for owner in (EI, EA)
        selected = owner.Formula(:default)
        @test validate(selected, air).equation.arguments[2:3] == (Val(1), Val(1))
        @test validate(selected, soil).equation.arguments[2:3] == (Val(2), Val(2))
        @test validate(selected, self).kind === :self
        @test validate(selected, soil).kind === :mutual
        @test validate(selected, mixed).kind === :mutual
        @test_throws ArgumentError validate(selected, E.EarthPair(
            1, 2, (-1.0, -2.0), 1.0, (2, 3)))
        @test_throws DimensionMismatch validate(selected, 3)
        @test validate(selected, 2) === selected
        author=owner.Formula(:xue2018)
        @test_throws ArgumentError validate(author, mixed)
        @test validate(selected, air).equation isa FM
        @test validate(selected,air).equation.selection === selected
    end
    for id in (:ametani2009, :lucca1994)
        selected = EI.Formula(id)
        @test validate(selected, mixed).kind === :mutual
        @test_throws ArgumentError validate(selected, air)
        @test_throws ArgumentError validate(selected, soil)
        @test_throws ArgumentError validate(selected, self)
    end
    vertical=E.EarthPair(1, 2, (-1.0, -2.0), 0.0, (2, 2))
    @test_throws DomainError validate(EI.Formula(:saad1996), vertical)
    @test validate(EI.Formula(:saad1996), self).kind === :self
    @test validate(EI.Formula(:wedepohl1973), vertical).kind === :mutual
    @test validate(EI.Formula(:wedepohl1973), self).kind === :self
    @test_throws ArgumentError validate(EI.Formula(:pollaczek1926), air)
    @test_throws ArgumentError validate(EI.Formula(:carson1926), soil)
end

@testitem "Engine / Wedepohl mutual distance includes the vertical displacement" tags=[:unit] begin
    const E = LineCableModels.Engine
    selected = E.EarthImpedance.Formula(:wedepohl1973)
    # Eq. (8) depends on axis distance and the sum of depths. Rotating the
    # displacement at fixed midpoint preserves both, including vertical axes.
    for T in (Float32, Float64, BigFloat), frequency in (50, 1000)

        rho = T[Inf, 100]
        epsilon = T(8.8541878128e-12) .* T[1, 10]
        mu = fill(4T(pi)*T(1e-7), 2)
        jω = complex(zero(T), 2T(pi)*frequency)
        horizontal = E.EarthPair(1, 2, (T(-2), T(-2)), one(T), (2, 2))
        reference = selected(rho, epsilon, mu, jω, horizontal)()
        for (h1, h2, separation) in ((-8//5, -12//5, 3//5), (-3//2, -5//2, 0//1))
            pair = E.EarthPair(1, 2, (T(h1), T(h2)), T(separation), (2, 2))
            actual = selected(rho, epsilon, mu, jω, pair)()
            @test isfinite(actual)
            @test actual ≈ reference rtol = 32eps(T)
            reversed = E.EarthPair(2, 1, reverse(pair.heights), pair.separation, (2, 2))
            @test selected(rho, epsilon, mu, jω, reversed)() == actual
        end
    end
    coincident = E.EarthPair(1, 2, (-2.0, -2.0), 0.0, (2, 2))
    @test_throws DomainError validate(selected, coincident)
end

@testitem "Engine / selected media retain one problem-owned scalar Γ contract" tags=[:unit] begin
    const E=LineCableModels.Engine
    μ0,ε0=4pi*1e-7,8.8541878128e-12
    rho,epsilon,mu=[Inf,100.0],[ε0,10ε0],[μ0,μ0]
    jω=complex(0.0,2pi*50)
    pair=E.EarthPair(1,2,(-1.0,-2.0),0.75,(2,2))
    for owner in (E.EarthImpedance,E.EarthAdmittance)
        selected=owner.Formula(:default)
        functor=selected(rho,epsilon,mu,jω,pair)
        @test functor.binding.equation.selection === selected
        for name in (:segments,:tolerance,:formula,:gamma_squared)
            @test !hasproperty(functor.state,name)
        end
        @test !hasfield(typeof(functor),:hooks)
        @test owner.Γ(functor)==0
        prescribed=selected(rho,epsilon,mu,jω,pair;Γ=1e-4im)
        @test owner.Γ(prescribed)==1e-4im
        @test prescribed.state.gamma_medium_squared==functor.state.gamma_medium_squared
        @test prescribed.state.gamma==sqrt.(jω .* mu .* (inv.(rho) .+ jω .* epsilon))
        for wrong in ([1.0],(Γ=0,squared=9),NaN,()->0)
            @test_throws ArgumentError selected(rho,epsilon,mu,jω,pair;Γ=wrong)
        end
        @test_throws ArgumentError owner.Formula(:default;parameters=(unknown=1,))
        @test_throws DimensionMismatch selected(
            [Inf,100.0,999.0],[ε0,10ε0,20ε0],[μ0,μ0,μ0],jω,pair)
    end
    carson=E.EarthImpedance.Formula(:carson1926)
    air=E.EarthPair(1,2,(1.0,2.0),1.0,(1,1))
    @test_throws ArgumentError carson(rho,epsilon,mu,jω,air;Γ=1e-4im)
    @test LineCableModels.ComputationOptions === NamedTuple
end

@testitem "Engine / indexed finite-layer methods determine admission and ordered assembly" tags=[:unit] setup=[FormulaContractModels] begin
    using LinearAlgebra
    const M=FormulaContractModels
    const E=M.E
    const EI=M.EI
    const EA=M.EA
    const EP=M.EP
    empty!(M.calls)
    inventories=(EI.formulas(), EA.formulas())
    selected=M.selection(EI; options = (integration = (method = :trapz, options = (;)),))
    heights=(2.0, -0.25, -1.5)
    pairs=[E.EarthPair(t, s, (heights[s], heights[t]), s==t ? 0.0 : 1.0,
               (s, t); radius = s==t ? 0.01 : nothing) for s in 1:3 for t in 1:3]
    bound=validate(selected, pairs)
    @test isempty(M.calls) # Structural preflight executes no kernels.
    @test count(case -> haskey(case.options, :integration), bound) == 1
    @test bound[1].options.integration.method === Val(:trapz)
    @test_throws ArgumentError validate(selected, pairs[2:end]) # Unconsumed integral controls.
    absent=E.EarthPair(1, 2, (-2.0, -3.0), 1.0, (3, 4))
    @test_throws ArgumentError validate(selected, absent)
    # A numerical specialization alone cannot admit a missing canonical case.
    EI.earth_impedance(::M.LayerImpedance, ::Val{:mutual}, ::Val{3}, ::Val{4},
        f::Float64, pair, workspace) = 0
    @test_throws ArgumentError validate(selected, absent)
    @test isempty(M.calls)

    material=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(CableDesign, "indexed-contract",
        Stack(Group(:phase, Region(:core, Disk(0.01), material))))
    system=build(LineCableSystem, fill(design, 3),
        [(Float64(i), h) for (i, h) in enumerate(heights)];
        connections = [(phase = i,) for i in 1:3])
    earth=build(EP.EarthModel, (
        EP.EarthLayer(100.0, 10.0, 1.0, 0.5), EP.EarthLayer(200.0, 20.0, 1.0)))
    problem=LineParametersProblem(system; earth_props = earth, frequencies = [50.0, 500.0])
    fd=M.DispersiveEarth(exponent=0)
    fd_calls=fd.seen
    formulation=Formulation(earth_impedance = selected, earth_admittance = M.selection(EA),
        earth_properties = fd,
        options = (ideal_transposition = false,))
    result=compute(problem, formulation; options = (trace = true,))
    @test length(fd_calls) == 4
    execution=computation_options(LineCableModelsCoaxial,(;))
    T=eltype(problem)
    blueprints=E.CableBlueprint{T}[E.flatten(LineCableModelsCoaxial(),d,T)
        for d in problem.system.designs]
    workspace=E.LineParametersWorkspace(problem,formulation,execution,blueprints)
    @test !workspace.buffers.uses_earth_systems
    @test length(M.calls) == 36
    for (_, _, row, column, layers, geometry, rho, physical) in M.calls
        @test layers == (column, row)
        @test geometry == (heights[column], heights[row])
        @test rho == [Inf, 50.0, 100.0]
        @test physical.layers == layers && physical.heights == geometry
    end
    for k in eachindex(problem.frequencies), t in 1:3, s in 1:3
        # Declared manufactured rule, evaluated independently of the probe method.
        coefficient=11s+17t+3t+5s+problem.frequencies[k]/100+(s==t ? 101 : 0)
        @test details(result).trace.Zg[t,s,k] ≈ coefficient*(1e-4+1e-3im) rtol=3e-6
        @test details(result).trace.Pg[t,s,k] ≈ coefficient*1e9
    end
    @test result.Z.values[1,2,1] ≈ (11*2+17+3+5*2+.5)*(1e-4+1e-3im)
    @test result.Z.values[2,1,1] ≈ (11+17*2+3*2+5+.5)*(1e-4+1e-3im)
    @test result.Y.values[:, :, 1] ≈ 2pi * 50im * inv(details(result).trace.P[:, :, 1])
    # A full-layer impedance and a separately reduced homogeneous potential
    # consume different material inventories in the same solve.
    buried=build(LineCableSystem, fill(design, 3),
        [(1.0, -0.1), (2.0, -0.25), (3.0, -1.5)];
        connections = [(phase = i,) for i in 1:3])
    buried_problem=LineParametersProblem(buried; earth_props = earth, frequencies = [50.0])
    for order in (:before, :after)
        empty!(M.calls)
        empty!(fd_calls)
        hybrid=Formulation(earth_impedance = M.selection(EI),
            earth_admittance = formula(:xue2018; equivalent_earth = formula(:default; order)),
            earth_properties = fd,
            options = (ideal_transposition = false,))
        value=compute(buried_problem, hybrid)
        @test all(isfinite, value.Z.values) && all(isfinite, value.Y.values)
        @test length(M.calls) == 9
        @test all(record -> length(record[7]) == 3, M.calls)
        @test Set(record[5] for record in M.calls) == Set(((2, 2), (2, 3), (3, 2), (3, 3)))
        @test length(fd_calls) == (order === :before ? 11 : 2)
        @test details(value).formulations.equivalent_earth.earth_impedance === nothing
        @test details(value).formulations.equivalent_earth.earth_admittance.order === order
    end
    @test (EI.formulas(), EA.formulas()) === inventories
end

@testitem "Engine / material laws and numerical options do not change physical routing" tags=[:unit] setup=[TestFixtures,FormulaContractModels] begin
    const E=LineCableModels.Engine
    base=TestFixtures.line_parameters_problem(frequencies=[50.0,1e5])
    fd=FormulaContractModels.DispersiveEarth(exponent=0)
    @test_throws ArgumentError compute(base,
        Formulation(earth_impedance=:carson1926,earth_properties=fd))
    @test isempty(fd.seen)
    selected=Formulation(earth_properties=fd)
    changed=compute(base,selected)
    material=base.earth_props.layers[2]
    equivalent=LineParametersProblem(base.system;temperature=base.temperature,
        earth_props=homogeneous(rho=material.rho/2,eps_r=material.eps_r,mu_r=material.mu_r),
        frequencies=base.frequencies)
    expected=compute(equivalent)
    @test changed.Z.values ≈ expected.Z.values rtol=1e-10
    @test changed.Y.values ≈ expected.Y.values rtol=1e-10
    @test Set(record[2] for record in fd.seen)==Set(base.frequencies)
    @test details(changed).formulations.effective.earth_properties === :DispersiveEarth
    reference=compute(base)
    for method in (:trapz,:cim)
        result=compute(base,Formulation(
            earth_impedance=formula(:default;options=(integration=(method=method,options=(;)),)),
            earth_admittance=formula(:default;options=(integration=(method=method,options=(;)),))))
        @test result.Z.values ≈ reference.Z.values rtol=3e-6
        @test result.Y.values ≈ reference.Y.values rtol=3e-6
        @test all(case->case.options.integration.method === Val(method),
            details(result).formulations.numerical.earth_impedance)
    end
end

@testitem "Engine / earth state and quadrature preserve scalar types and uncertainty" tags=[:unit] begin
    using Measurements
    const E=LineCableModels.Engine
    for T in (Float32, Float64, BigFloat)
        rho=T[Inf, 100]
        epsilon=T(8.8541878128e-12) .* T[1, 10]
        mu=fill(T(4pi*1e-7), 2)
        pair=E.EarthPair(1, 2, (-one(T), -T(2)), T(0.75), (2, 2))
        s=complex(zero(T), T(2)*T(pi)*T(50))
        for owner in (E.EarthImpedance, E.EarthAdmittance), method in (:quad, :trapz)

            options=(integration = (method = method,
                options = (rtol = T===Float32 ? 1e-4 : 1e-6,)),)
            functor=owner.Formula(:xue2018; options)(rho, epsilon, mu, s, pair)
            @test functor.state.Γ isa Complex{T}
            result=functor()
            @test result isa Complex{T}
            @test isfinite(result)
        end
    end
    rho=measurement.([Inf, 100.0], [0.0, 1.0])
    epsilon=measurement.(8.8541878128e-12 .* [1, 10], 0.0)
    mu=measurement.(fill(4pi*1e-7, 2), 0.0)
    s=complex(measurement(0.0, 0.0), measurement(2pi*50, 0.0))
    pair=E.EarthPair(1, 2, (-1.0, -2.0), 0.75, (2, 2))
    for owner in (E.EarthImpedance, E.EarthAdmittance)
        value=owner.Formula(:xue2018)(rho, epsilon, mu, s, pair)()
        @test value isa Complex{Measurement{Float64}}
        @test isfinite(value)
        @test uncertainty(real(value)) > 0
        @test uncertainty(imag(value)) > 0
    end
end
