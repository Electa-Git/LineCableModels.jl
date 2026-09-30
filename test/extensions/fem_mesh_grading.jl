@testitem "Gmsh FEM / exterior grading preserves local controls and material interfaces" tags=[:extension] begin
    using Gmsh
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    copper = Material(kind=:conductor, rho=1/5.8e7)
    wire = build(CableDesign, "graded-wire",
        Stack(Group(:core, Region(:metal, Disk(0.0425), copper))))
    system = build(LineCableSystem, [wire, wire], [Pose2(0.,1.), Pose2(1.,1.)];
        connections=[Dict(:core=>1), Dict(:core=>2)], line_length=1.)
    problem = LineParametersProblem(system; frequencies=[1e6],
        earth_props=homogeneous(rho=100., eps_r=1., mu_r=1.))
    formulation = Formulation(:LineCableModelsFEM; options=(physics=:quasi_fw,))
    controls = ComputationOptions(domain_skin_depths=24., pml_layers=8, mesh_size_factor=3.)
    baseline = extension._resolved_fem_model(problem, formulation,
        computation_options(LineCableModelsFEM, controls))
    graded = extension._resolved_fem_model(problem, formulation,
        computation_options(LineCableModelsFEM,
            ComputationOptions(; controls.data..., exterior_mesh_size_factor=8.)))
    a, b = only(baseline.mesh_plans), only(graded.mesh_plans)
    for name in (:frequency, :domain_halfwidth, :pml_thickness, :pml_strength,
            :pml_layers, :pml_grading, :volume_quadrature, :domain_mesh_size, :interface_mesh_size,
            :cable_interface_mesh_sizes, :wave_mesh_sizes, :wave_decay_radii)
        @test getproperty(a, name) == getproperty(b, name)
    end
    @test baseline.cable_outer_mesh_sizes == graded.cable_outer_mesh_sizes
    @test getproperty.(baseline.region_plans, :mesh_size) ==
        getproperty.(graded.region_plans, :mesh_size)
    @test all(>(a.domain_mesh_size), b.exterior_mesh_sizes)
    @test extension._mesh_fingerprint(baseline, "test") !=
        extension._mesh_fingerprint(graded, "test")
    for (length, near, far) in ((24.,0.15,1.2), (1.,0.2,0.2), (3.,1e-4,2.),
            (1.,0.17,nextfloat(nextfloat(0.17))))
        count, ratio = extension._exterior_edge_grading(length, near, far)
        widths = ratio .^ (0:count-1)
        widths .*= length / sum(widths)
        @test sum(widths) ≈ length
        @test first(widths) <= near*(1+1e-12)
        @test last(widths) <= far*(1+1e-12)
        @test issorted(widths)
    end
    counts = Int[]
    contour_counts = Vector{Int}[]
    pml_counts = Int[]
    mktempdir() do directory
        for model in (baseline, graded)
            session = extension._start_gmsh(0)
            try
                geometry = extension._build_geometry!(model, "grading-test")
                extension._configure_mesh!(model, geometry)
                Gmsh.gmsh.model.mesh.generate(2)
                path = joinpath(directory, "mesh$(length(counts)).msh")
                Gmsh.gmsh.write(path)
                extension._inspect_loaded_mesh(model, path)
                push!(counts, length(first(Gmsh.gmsh.model.mesh.get_nodes())))
                push!(pml_counts,sum(length(block) for surface in geometry.pml_surfaces
                    for block in Gmsh.gmsh.model.mesh.get_elements(2,surface)[2]))
                push!(contour_counts, [sum(length(first(Gmsh.gmsh.model.mesh.get_nodes(1,c)))
                    for c in curves) for curves in geometry.terminal_curves])
                for i in eachindex(model.terminal_ids)
                    curves = Gmsh.gmsh.model.get_entities_for_physical_group(1,
                        model.tags.voltage_path_base+i)
                    @test !isempty(curves)
                    @test sum(length(block) for c in curves for block in
                        Gmsh.gmsh.model.mesh.get_elements(1,c)[2]) > 0
                end
            finally
                extension._finish_gmsh(session)
            end
        end
        entry = export_data(:onelab, problem, formulation;
            file_name=joinpath(directory,"graded","study.pro"),
            mesh_options=(; controls.data..., exterior_mesh_size_factor=8.))
        session = extension._start_gmsh(0)
        try
            Gmsh.gmsh.open(replace(entry, r"\.pro$"=>".geo"))
            Gmsh.gmsh.model.mesh.generate(2)
            @test extension._inspect_loaded_mesh(graded,"native-grading") === nothing
            surfaces = Gmsh.gmsh.model.get_entities_for_physical_group(2,graded.tags.pml)
            # Preserve the structured exterior exactly. Native arithmetic can
            # perturb an unstructured conductor triangulation by one node.
            @test sum(length(block) for s in surfaces for block in
                Gmsh.gmsh.model.mesh.get_elements(2,s)[2]) == pml_counts[2]
        finally
            extension._finish_gmsh(session)
        end
    end
    @test counts[2] < counts[1]
    @test contour_counts[1] == contour_counts[2]
end

@testitem "Gmsh FEM / physical PML native export nodes and buried path" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"physical-pml-export",terminal(:core,
        core(Material(kind=:conductor,rho=1/5.8e7);r=.01)))
    system = build(LineCableSystem,[wire,wire],[(0.,1.),(1.,-1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[1e6],
        earth_props=homogeneous(rho=.1,eps_r=1.))
    form = Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,))
    controls = (domain_skin_depths=24.,mesh_size_factor=3.,exterior_mesh_size_factor=8.,
        pml_resolution=(interpolation_cells=72,coefficient_change=.12))
    model = FEM._resolved_fem_model(problem,form,
        computation_options(LineCableModelsFEM,ComputationOptions(;controls...)))
    plan = only(model.mesh_plans)
    coordinates = Dict{Int,Matrix{Float64}}()
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"physical-pml-nodes")
        FEM._configure_mesh!(model,geometry)
        gmsh.model.mesh.generate(2)
        for (curve,(count,_)) in geometry.transfinite_curves
            tags,xyz,_ = gmsh.model.mesh.get_nodes(1,curve,true)
            @test length(tags) == count
            coordinates[curve] = reshape(xyz,3,:)
            for adjacent in first(gmsh.model.get_adjacencies(1,curve))
                nodes = Set(Iterators.flatten(gmsh.model.mesh.get_elements(2,adjacent)[3]))
                @test all(in(nodes),tags)
            end
        end
        path = gmsh.model.get_entities_for_physical_group(1,model.tags.voltage_path_base+2)
        bottom = filter(path) do curve
            lo,hi = gmsh.model.get_parametrization_bounds(1,curve)
            maximum(gmsh.model.get_value(1,curve,[only(t)])[2] for t in (lo,hi)) <=
                -plan.domain_halfwidth+1e-10
        end
        @test length(bottom) == length(plan.pml_strips[3])
        @test all(c -> length(first(gmsh.model.get_adjacencies(1,c)))==2,bottom)
    finally
        FEM._finish_gmsh(session)
    end
    mktempdir() do directory
        entry = export_data(:onelab,problem,form;
            file_name=joinpath(directory,"study.pro"),mesh_options=controls)
        session = FEM._start_gmsh(0)
        try
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            gmsh.model.mesh.generate(2)
            for (curve,expected) in coordinates
                _,xyz,_ = gmsh.model.mesh.get_nodes(1,curve,true)
                actual = reshape(xyz,3,:)
                axis = maximum(expected[1,:])-minimum(expected[1,:])>0 ? 1 : 2
                @test actual[:,sortperm(vec(actual[axis,:]))] ≈
                    expected[:,sortperm(vec(expected[axis,:]))] rtol=1e-12 atol=1e-10
            end
            @test FEM._inspect_loaded_mesh(model,"physical-pml-export") === nothing
        finally
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / exterior grading retains the bulk cap inside large cores" tags=[:extension] begin
    using Gmsh
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    wire = build(CableDesign, "large-core", Stack(Group(:core,
        Region(:metal, Disk(3.), Material(kind=:conductor, rho=1/5.8e7)))))
    system = build(LineCableSystem, wire, Pose2(0.,4.); connections=Dict(:core=>1))
    problem = LineParametersProblem(system; frequencies=[1e6],
        earth_props=homogeneous(rho=4pi^2*0.1, eps_r=1.))
    caps = Float64[]
    for factor in (1.,8.)
        options = computation_options(LineCableModelsFEM, ComputationOptions(
            domain_skin_depths=24., pml_layers=4, mesh_size_factor=3.,
            exterior_mesh_size_factor=factor))
        model = extension._resolved_fem_model(problem, LineCableModelsFEM(), options)
        @test model.fine_mesh_size > only(model.mesh_plans).domain_mesh_size
        session = extension._start_gmsh(0)
        try
            geometry = extension._build_geometry!(model, "large-core-test")
            extension._configure_mesh!(model, geometry)
            # Inspect the actual native size constraints. Boundary-layer
            # triangulation can differ by a node for identical targets; a
            # node count does not establish that the bulk cap was retained.
            surface = only(only(geometry.material_surfaces))
            constants = Float64[]
            for field in Gmsh.gmsh.model.mesh.field.list()
                Gmsh.gmsh.model.mesh.field.get_type(field) == "Restrict" || continue
                surface in Gmsh.gmsh.model.mesh.field.get_numbers(field,"SurfacesList") || continue
                source = Int(Gmsh.gmsh.model.mesh.field.get_number(field,"InField"))
                Gmsh.gmsh.model.mesh.field.get_type(source) == "MathEval" || continue
                value = tryparse(Float64,Gmsh.gmsh.model.mesh.field.get_string(source,"F"))
                value === nothing || push!(constants,value)
            end
            cap = min(Gmsh.gmsh.option.get_number("Mesh.MeshSizeMax"),minimum(constants))
            @test cap <= only(model.mesh_plans).domain_mesh_size
            push!(caps,cap)
        finally
            extension._finish_gmsh(session)
        end
    end
    @test caps[1] == caps[2]
end
@testitem "Gmsh FEM / fixed PML progression and native node positions" tags=[:extension] begin
    using Gmsh
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    g0 = (192/191)*log(1536)
    @test FEM._pml_progression(192,g0,1.,2.) ≈ (8*192)^(1/191) rtol=2eps()
    @test FEM._pml_progression(1,1e300,1.,2.) == 1.
    @test FEM._pml_progression(8,1e-20,1.,2.) == 1.
    for args in ((8,1e300,1.,2.), (8,100.,1.,2.), (8,0.,1.,nextfloat(1.)),
            (8,0.,1.,1.), (8,0.,1.,Inf))
        @test_throws ArgumentError FEM._pml_progression(args...)
    end
    session = FEM._start_gmsh(0)
    try
        gmsh.model.add("normal-progression")
        # Independent analytic coordinates, both curve orientations and counts.
        for (k,(n,g)) in enumerate(((1,g0),(8,0.),(8,1e-20),(5,2.),(128,g0),(192,g0)))
            a = gmsh.model.geo.add_point(1.,Float64(k),0.)
            b = gmsh.model.geo.add_point(2.,Float64(k),0.)
            reversed = iseven(k)
            curve = gmsh.model.geo.add_line(reversed ? b : a,reversed ? a : b)
            r = FEM._pml_progression(n,g,1.,2.)
            gmsh.model.geo.mesh.set_transfinite_curve(curve,n+1,"Progression",reversed ? inv(r) : r)
        end
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(1)
        for (curve,(n,g)) in enumerate(((1,g0),(8,0.),(8,1e-20),(5,2.),(128,g0),(192,g0)))
            xyz = reshape(gmsh.model.mesh.get_nodes(1,curve,true)[2],3,:)
            actual = sort(vec(xyz[1,:])) .- 1
            expected = g == 0 ? collect(0:n)./n : expm1.(g.*(0:n)./n)./expm1(g)
            @test length(actual) == n+1
            # Bound each normalized node error, independent of node count.
            @test maximum(abs,actual.-expected) <= 2e-7
        end
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / directional PML topology export and mesh identity" tags=[:extension] begin
    using Gmsh
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"pml-controls",terminal(:core,
        core(Material(kind=:conductor,rho=1.72e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[1e4],earth_props=homogeneous(rho=100.))
    form = Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,))
    resolve(controls) = FEM._resolved_fem_model(problem,form,
        computation_options(LineCableModelsFEM,ComputationOptions(controls)))
    controls = (pml_layers=(5,4,3),pml_grading=(2.,0.,1.))
    model = resolve(controls)
    plan = only(model.mesh_plans)
    @test plan.pml_layers == (5,4,3)
    @test plan.pml_grading == (2.,0.,1.)
    scalar = resolve((pml_layers=5,pml_grading=2.))
    tupled = resolve((pml_layers=(5,5,5),pml_grading=(2.,2.,2.)))
    @test FEM._mesh_fingerprint(scalar,"test") == FEM._mesh_fingerprint(tupled,"test")
    key = FEM._mesh_fingerprint(model,"test")
    for direction in 1:3, field in (:pml_layers,:pml_grading)
        value = getproperty(controls,field)
        changed = merge(controls,NamedTuple{(field,)}((Base.setindex(value,value[direction]+1,direction),)))
        @test FEM._mesh_fingerprint(resolve(changed),"test") != key
    end
    metadata = FEM._mesh_metadata(model,key,"test",:generated)
    serialized = FEM.JSON3.read(FEM.JSON3.write(metadata))
    @test Tuple(serialized.pml_layers) == (5,4,3)
    @test Tuple(serialized.pml_grading) == (2.,0.,1.)
    for property in (:domain_halfwidth,:pml_thickness,:pml_strength,:domain_mesh_size,
            :exterior_mesh_sizes,:interface_mesh_size,:volume_quadrature)
        @test getproperty(plan,property) == getproperty(only(scalar.mesh_plans),property)
    end
    @test model.conductor_mesh == scalar.conductor_mesh
    @test getproperty.(model.region_plans,:mesh_size) == getproperty.(scalar.region_plans,:mesh_size)
    coordinates = Dict{Int,Matrix{Float64}}()
    constraints = Dict{Int,Tuple{Int,Float64}}()
    mktempdir() do root
        session = FEM._start_gmsh(0)
        try
            geometry = FEM._build_geometry!(model,"directional-pml")
            FEM._configure_mesh!(model,geometry)
            gmsh.model.mesh.generate(2)
            merge!(constraints,geometry.transfinite_curves)
            corners = 0
            for surface in geometry.pml_surfaces
                _, vertices = geometry.transfinite_surfaces[surface]
                points = [gmsh.model.get_value(0,p,Float64[]) for p in vertices]
                xmid = sum(p[1] for p in points)/4 - model.centre[1]
                ymid = sum(p[2] for p in points)/4
                if abs(xmid) > plan.domain_halfwidth && abs(ymid) > plan.domain_halfwidth
                    elements = gmsh.model.mesh.get_elements(2,surface)
                    @test elements[1] == [2] # first-order triangles
                    corners += length(only(elements[2]))
                end
            end
            @test corners == 4*5*(4+3)
            directions = zeros(Int,3)
            for (curve,(count,ratio)) in constraints
                lo,hi = gmsh.model.get_parametrization_bounds(1,curve)
                a,b = [gmsh.model.get_value(1,curve,[only(t)]) for t in (lo,hi)]
                xmid,ymid = (a[1]+b[1])/2-model.centre[1],(a[2]+b[2])/2
                axis,direction = if abs(xmid)>plan.domain_halfwidth && abs(a[1]-b[1])>0
                    (1,1)
                elseif abs(ymid)>plan.domain_halfwidth && abs(a[2]-b[2])>0
                    (2,ymid>0 ? 2 : 3)
                else
                    (0,0)
                end
                tags,xyz,_ = gmsh.model.mesh.get_nodes(1,curve,true)
                points = reshape(xyz,3,:)
                coordinates[curve] = points[:,sortperm(vec(points[axis==0 ? 1 : axis,:]))]
                @test length(tags) == count
                axis==0 && continue
                directions[direction] += 1
                n,g = controls.pml_layers[direction],controls.pml_grading[direction]
                @test count == n+1
                onset = axis==1 ? model.centre[1]+sign(xmid)*plan.domain_halfwidth : sign(ymid)*plan.domain_halfwidth
                actual = sort(abs.(vec(points[axis,:]).-onset)./plan.pml_thickness[direction])
                expected = g==0 ? collect(0:n)./n : expm1.(g.*(0:n)./n)./expm1(g)
                @test maximum(abs,actual.-expected) <= 2e-7
                # Every adjacent patch uses the same curve nodes.
                for adjacent in first(gmsh.model.get_adjacencies(1,curve))
                    field_nodes = Set(Iterators.flatten(gmsh.model.mesh.get_elements(2,adjacent)[3]))
                    @test all(in(field_nodes),tags)
                end
            end
            @test all(>(0),directions)
            path = gmsh.model.get_entities_for_physical_group(1,model.tags.voltage_path_base+2)
            bottom = filter(path) do curve
                lo,hi = gmsh.model.get_parametrization_bounds(1,curve)
                y = [gmsh.model.get_value(1,curve,[only(t)])[2] for t in (lo,hi)]
                maximum(y) <= -plan.domain_halfwidth + 1e-10
            end
            @test length(bottom) == 1
            @test constraints[only(bottom)][1] == controls.pml_layers[3]+1
            @test length(first(gmsh.model.get_adjacencies(1,only(bottom)))) == 2
            # Exercise actual cache reuse with the resolved directional controls.
            gmsh.write(joinpath(root,"model.msh"))
            cached = FEM._mesh_fingerprint(model,FEM._gmsh_version())
            FEM._cache_mesh!(joinpath(root,"model.msh"),joinpath(root,"meshes",cached,"model.msh"),
                joinpath(root,"meshes",cached,"mesh.json"),FEM._mesh_metadata(model,cached,FEM._gmsh_version(),:generated))
            run = FEM._create_run(root)
            execution = computation_options(LineCableModelsFEM,ComputationOptions(controls))
            FEM._select_mesh!(run,model,geometry,execution,root)
            @test run.mesh_source == :cache
        finally
            FEM._finish_gmsh(session)
        end
        entry = export_data(:onelab,problem,form;file_name=joinpath(root,"detached","study.pro"),mesh_options=controls)
        session = FEM._start_gmsh(0)
        try
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            gmsh.model.mesh.generate(2)
            for (curve,expected) in coordinates
                _,xyz,_ = gmsh.model.mesh.get_nodes(1,curve,true)
                points = reshape(xyz,3,:)
                axis = maximum(expected[1,:])-minimum(expected[1,:])>0 ? 1 : 2
                actual = points[:,sortperm(vec(points[axis,:]))]
                # Sort both, including vertical tangential edges.
                wanted = expected[:,sortperm(vec(expected[axis,:]))]
                @test actual ≈ wanted rtol=1e-12 atol=1e-10
            end
            @test FEM._inspect_loaded_mesh(model,"detached-directional") === nothing
        finally
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / physical PML preserves qualified native strips" tags=[:extension] begin
    using Gmsh, TOML
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    copper = Material(kind=:conductor,rho=1/5.8e7)
    wire = build(CableDesign,"physical-pml-wire",Stack(Group(:core,Region(:metal,Disk(.085),copper))))
    system = build(LineCableSystem,[wire,wire],[Pose2(0.,1.),Pose2(1.,1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)],line_length=1.)
    frozen = TOML.parsefile(joinpath(@__DIR__,"../fixtures/data/fem/pml_physical_strips.toml"))["cases"]
    problem = LineParametersProblem(system;frequencies=[c["frequency_hz"] for c in frozen],
        earth_props=homogeneous(rho=.1,eps_r=1.,mu_r=1.))
    form = Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,))
    controls = (domain_skin_depths=24.,mesh_size_factor=3.,exterior_mesh_size_factor=8.,
        pml_resolution=(interpolation_cells=72,coefficient_change=.12))
    options = computation_options(LineCableModelsFEM,ComputationOptions(;controls...))
    model = FEM._resolved_fem_model(problem,form,options)
    for (plan,expected) in zip(model.mesh_plans,frozen)
        @test plan.pml_grading === nothing
        @test plan.pml_layers == map(s -> sum(x.count for x in s),plan.pml_strips)
        for (direction,strips) in zip(("side","top","bottom"),plan.pml_strips)
            @test length(strips) == length(expected[direction])
            for (strip,row) in zip(strips,expected[direction])
                @test strip.count == row["count"]
                @test strip.start ≈ row["start"] rtol=2e-12 atol=1e-14
                @test strip.stop ≈ row["stop"] rtol=2e-12 atol=1e-14
                @test strip.ratio ≈ row["ratio"] rtol=2e-12 atol=1e-14
            end
        end
    end
    other = FEM._resolved_fem_model(problem,form,computation_options(LineCableModelsFEM,
        ComputationOptions(;controls...,pml_resolution=(interpolation_cells=72,coefficient_change=.1))))
    @test FEM._mesh_fingerprint(model,"test") != FEM._mesh_fingerprint(other,"test")
    # The public export must serialize the same native strip constraints.
    mktempdir() do directory
        entry = export_data(:onelab,problem,form;file_name=joinpath(directory,"study.pro"),mesh_options=controls)
        @test isfile(entry)
        for (index,plan) in enumerate(model.mesh_plans)
            path = joinpath(directory,"geometry","case-"*lpad(index,4,'0')*".geo")
            @test isfile(path)
            text = read(path,String)
            @test occursin("Transfinite Surface",text)
            for strips in plan.pml_strips, strip in strips
                @test occursin(" = $(strip.count+1) Using Progression ",text)
            end
        end
    end
end
