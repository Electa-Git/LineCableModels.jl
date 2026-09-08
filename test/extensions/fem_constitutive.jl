@testitem "Gmsh FEM / constitutive selection, evaluated state and rejection" tags=[:extension] begin
    using Gmsh, Measurements
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const FD = LineCableModels.Earth.FrequencyDependent
    const TD = LineCableModels.Materials.TemperatureDependent
    copper = Material(:conductor, 1.72e-8, 1, 1, 20, 0.004)
    dielectric = Material(:insulator, 1e8, 2.3, 1, 20, -0.003; tan_delta=0.025)
    design = build(CableDesign, "constitutive", terminal(:core,
        core(copper; r=0.005), insulation(dielectric; t=0.005)))
    system = build(LineCableSystem, design, (0.0,-0.1); connections=Dict(:core=>1))
    problem = LineParametersProblem(system; temperature=80.0, frequencies=[50.0,1000.0],
        earth_props=homogeneous(rho=100.0, eps_r=10.0))
    calls = Float64[]
    soil_law = (m,f,p,o,w) -> begin
        push!(calls, f)
        EarthMaterial(m.rho / (1+f/1000), m.eps_r*(1+f/2000), m.mu_r*(1+f/10000))
    end
    temperature_law = (m,t,p,o,w) -> m.rho * exp((t-m.T0)/1000)
    dielectric_calls = Ref(0)
    dielectric_law = (m,f,t,p,o,w) -> begin
        dielectric_calls[] += 1
        @test m.T0 == t
        complex(inv(m.rho), 2pi*f*8.8541878128e-12*m.eps_r)
    end
    for (method, law) in ((FD.earth_material, soil_law),
            (TD.temperature_resistivity, temperature_law),
            (LineCableModels.Engine.InsulationAdmittance.insulation_material, dielectric_law))
        @eval LineCableModels.computation_options(
            ::LineCableModels.FormulaMethod{:default, typeof($method)}, ::$(typeof(law))) = (;)
    end
    formulation = LineCableModelsFEM(
        earth_properties=formula(:default; hooks=(contribution=soil_law,)),
        temperature_dependence=formula(:default; hooks=(contribution=temperature_law,)),
        insulation_admittance=formula(:default; hooks=(contribution=dielectric_law,)))
    model = FEM._resolved_fem_model(problem, formulation)
    @test calls == problem.frequencies
    @test dielectric_calls[] == 2
    @test keys(formulation.methods) == keys(formulation.definitions) ==
        (:insulation_admittance, :semicon_admittance, :earth_properties, :temperature_dependence)
    metal, passive = model.material_plans
    @test real.(metal.admittivity) ≈ fill(inv(copper.rho*exp(0.06)), 2)
    @test real.(passive.admittivity) ≈ fill(inv(dielectric.rho*exp(0.06)), 2)
    @test problem.system.designs[1].geometry.regions[1].source.material.rho == copper.rho
    for (index,f) in pairs(problem.frequencies)
        earth = model.earth_materials[index]
        @test (earth.rho,earth.eps_r,earth.mu_r) == (100/(1+f/1000),10*(1+f/2000),1+f/10000)
        @test model.mesh_plans[index].domain_radius ≈ sqrt(earth.rho/(pi*f*earth.mu_r*4pi*1e-7))
    end
    record = FEM.formulation_record(formulation)
    @test record.selections.temperature_dependence.identifier === :default
    @test !record.selections.temperature_dependence.replayable
    @test record.selections.temperature_dependence.hooks.contribution.type == string(typeof(temperature_law))
    @test record.selections.semicon_admittance.replayable
    for name in (:internal_impedance,:insulation_impedance,:earth_impedance,:earth_admittance,:pipe_impedance)
        @test_throws MethodError Formulation(:LineCableModelsFEM; NamedTuple{(name,)}((formula(:default),))...)
    end
    @test_throws ArgumentError LineCableModelsFEM(earth_properties=formula(:default; equivalent_earth=:default))
    @test_throws MethodError LineCableModelsFEM(earth_properties=LineCableModels.Earth.EquivalentHomogeneous.Formula(:default))
    finite_soil = EarthModel(100.0,10.0,1.0; thickness=5.0)
    finite_problem = LineParametersProblem(system; frequencies=[50.0], earth_props=finite_soil)
    @test_throws LineCableModelsFEMError FEM._resolved_fem_model(finite_problem, LineCableModelsFEM())
    hot = LineParametersProblem(system; temperature=250.0, frequencies=[50.0],
        earth_props=problem.earth_props)
    @test_throws LineCableModelsFEMError FEM._resolved_fem_model(hot, LineCableModelsFEM())
    @test FEM._resolved_fem_model(hot, LineCableModelsFEM(temperature_dependence=nothing)) isa FEM.FEMResolvedModel
    @test FEM._fem_float64_scalar(measurement(big"1.25",big"0.01")) === 1.25
    @test_throws OverflowError FEM._fem_float64_scalar(big"1e400")
    conducting_limit = LineParametersProblem(system;frequencies=[50.0],
        earth_props=homogeneous(rho=1e-309))
    @test_throws LineCableModelsFEMError FEM._resolved_fem_model(conducting_limit,LineCableModelsFEM())
end

@testitem "Gmsh FEM / metallic enclosure solves and reduces without a pipe formula" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh, LinearAlgebra
    copper = Material(:conductor,1.72e-8,1,1,20,0.004)
    dielectric = Material(:insulator,1e14,2.3)
    air = Material(:insulator,Inf,1.0)
    function enclosed_system(metal)
        first_core = terminal(:a,core(metal;r=0.005),insulation(dielectric;t=0.002))
        second_core = terminal(:b,core(metal;r=0.005),insulation(dielectric;t=0.002))
        wall = terminal(:pipe,sheath(metal;t=0.001),insulation(dielectric;t=0.002))
        design = build(CableDesign,"constitutive-enclosure",
            pipe(at(first_core,-0.01,0),at(second_core,0.01,0);
                shape=Disk(0.025),fill=air,wall))
        build(LineCableSystem,design,Pose2(0.0,-1.0);
            connections=Dict(:a=>1,:b=>2,:pipe=>0))
    end
    problem = LineParametersProblem(enclosed_system(copper);temperature=80.0,
        frequencies=[50.0],earth_props=homogeneous(rho=100.0))
    selected = LineCableModelsFEM(options=(ideal_transposition=false,),
        fem_options=(gmsh_verbosity=0,getdp_verbosity=0))
    result = compute(problem,selected;options=(trace=true,))
    primitive = result.details.fem.primitive
    @test primitive.phase_map == [1,2,0]
    @test size(primitive.Z_primitive) == (3,3,1)
    @test result.details.fem.reduced_phase_map == [1,2]
    @test all(isfinite,result.Z.values) && all(isfinite,result.Y.values)
    z = primitive.Z_primitive[:,:,1]
    p = primitive.P_primitive[:,:,1]
    # Independently eliminate the grounded wall's current/charge unknowns.
    @test result.Z.values[:,:,1] ≈ z[1:2,1:2]-z[1:2,3:3]*(z[3:3,3:3]\z[3:3,1:2])
    @test result.Y.values[:,:,1] ≈ inv(p[1:2,1:2]-p[1:2,3:3]*(p[3:3,3:3]\p[3:3,1:2]))
    @test z ≈ transpose(z) rtol=1e-10
    @test p ≈ transpose(p) rtol=1e-8
    @test eigmin(Symmetric(real(z))) > 0
    @test eigmin(Symmetric(imag(result.Y.values[:,:,1]))) > 0
    corrected = Material(:conductor,1.72e-8*(1+0.004*60))
    reference_problem = LineParametersProblem(enclosed_system(corrected);
        temperature=80.0,frequencies=[50.0],earth_props=problem.earth_props)
    reference = compute(reference_problem,LineCableModelsFEM(temperature_dependence=nothing,
        options=selected.options,fem_options=selected.execution);options=(trace=true,))
    @test z ≈ reference.details.fem.primitive.Z_primitive[:,:,1] rtol=2e-9
    @test result.Y.values ≈ reference.Y.values rtol=2e-9
end

@testitem "Gmsh FEM / real constitutive laws match independent static material solves" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const FD = LineCableModels.Earth.FrequencyDependent
    copper = Material(:conductor, 1.72e-8, 1, 1, 20, 0.004)
    dielectric = Material(:insulator, 1e8, 2.3, 1, 20, -0.003; tan_delta=0.025)
    function system_for(metal, passive)
        design = build(CableDesign,"constitutive-numerical",terminal(:core,
            core(metal; r=0.005), insulation(passive; t=0.005)))
        build(LineCableSystem,[design,design],[(0.0,-0.1),(0.1,-0.1)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
    end
    system = system_for(copper,dielectric)
    options = (reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    execution = (gmsh_verbosity=0,getdp_verbosity=0)
    for law in (:default,:Ametani2004), temperature in (20.0,80.0)
        problem = LineParametersProblem(system; temperature, frequencies=[50.0],
            earth_props=homogeneous(rho=100.0,eps_r=10.0))
        selected = LineCableModelsFEM(insulation_admittance=law; options, fem_options=execution)
        actual = compute(problem,selected)
        corrected = map((copper,dielectric)) do material
            rho = material.rho*(1+material.alpha*(temperature-material.T0))
            Material(material.kind,rho,material.eps_r,material.mu_r,material.T0,0;
                tan_delta=material.tan_delta)
        end
        reference_problem = LineParametersProblem(system_for(corrected...);
            temperature, frequencies=[50.0],earth_props=problem.earth_props)
        reference = compute(reference_problem,LineCableModelsFEM(insulation_admittance=law,
            temperature_dependence=nothing; options,fem_options=execution))
        @test actual.Z.values ≈ reference.Z.values rtol=2e-9
        @test actual.Y.values ≈ reference.Y.values rtol=2e-9
    end
    soil_law = (m,f,p,o,w) -> EarthMaterial(m.rho/(1+f/1000),
        m.eps_r*(1+f/2000),m.mu_r*(1+f/10000))
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{:default,typeof(FD.earth_material)}, ::$(typeof(soil_law))) = (;)
    air = EarthLayer(Inf,1.5,1.2,Inf)
    earth = EarthModel(100.0,10.0,1.0; air_layer=air)
    problem = LineParametersProblem(system; frequencies=[50.0,1000.0],earth_props=earth)
    selected = LineCableModelsFEM(earth_properties=formula(:default;hooks=(contribution=soil_law,));
        options,fem_options=execution)
    actual = compute(problem,selected;options=(trace=true,))
    for (index,f) in pairs(problem.frequencies)
        # Explicit published inputs form the reference; no call to the fitted/selected law.
        static = EarthModel(100/(1+f/1000),10*(1+f/2000),1+f/10000; air_layer=air)
        reference_problem = LineParametersProblem(system;frequencies=[f],earth_props=static)
        reference = compute(reference_problem,LineCableModelsFEM(earth_properties=nothing;
            options,fem_options=execution);options=(trace=true,))
        @test actual.Z.values[:,:,index] ≈ reference.Z.values[:,:,1] rtol=2e-9
        @test actual.Y.values[:,:,index] ≈ reference.Y.values[:,:,1] rtol=2e-9
        @test actual.details.fem.inputs.mesh_plans[index].domain_radius ==
            only(reference.details.fem.inputs.mesh_plans).domain_radius
    end
    vacuum = LineParametersProblem(system;frequencies=problem.frequencies,
        earth_props=homogeneous(rho=100.0,eps_r=10.0))
    other = compute(vacuum,selected)
    @test !isapprox(actual.Z.values,other.Z.values;rtol=1e-5)
    # Buried cable capacitance is dominated by insulation. Exercise the air
    # permittivity response with the same cables placed above the interface.
    overhead = build(LineCableSystem,system.designs,[(0.0,1.0),(0.1,1.0)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    declared_air = compute(LineParametersProblem(overhead;frequencies=[50.0],earth_props=earth),selected)
    vacuum_air = compute(LineParametersProblem(overhead;frequencies=[50.0],earth_props=vacuum.earth_props),selected)
    @test !isapprox(declared_air.Y.values,vacuum_air.Y.values;rtol=1e-5)
end
