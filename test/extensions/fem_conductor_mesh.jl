@testitem "Gmsh FEM / conductor grading preserves interfaces and native edits" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    # A 3 mm strand and a much larger wire with a different material.
    designs = [build(CableDesign, "wire-$i", terminal(:core,
        core(Material(kind=:conductor, rho=rho); r=radius)))
        for (i, (radius, rho)) in enumerate(((.0015, 2e-8), (.0425, 8e-8)))]
    system = build(LineCableSystem, designs, [(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system; frequencies=[.1,1e6],
        earth_props=homogeneous(rho=100.,eps_r=1.))
    formulation = LineCableModelsFEM()
    controls = (pml_layers=4,mesh_size_factor=3.,exterior_mesh_size_factor=8.)
    resolved(extra=(;)) = FEM._resolved_fem_model(problem,formulation,
        computation_options(LineCableModelsFEM,ComputationOptions(;controls...,extra...)))
    model = resolved()
    @test FEM._conductor_circle_segments(model.conductor_mesh.geometry_tolerance) == 96
    finer = resolved((conductor_geometry_tolerance=2.5e-4,))
    @test FEM._conductor_circle_segments(finer.conductor_mesh.geometry_tolerance) == 192
    @test FEM._mesh_fingerprint(model,"test") != FEM._mesh_fingerprint(finer,"test")
    changed_material = deepcopy(model)
    changed_material.material_plans[1].admittivity .*= 2
    @test FEM._mesh_fingerprint(model,"test") != FEM._mesh_fingerprint(changed_material,"test")
    for region in model.region_plans
        material = model.material_plans[region.material_index]
        low, high = [FEM._conductor_mesh_sizes(model,region,p) for p in model.mesh_plans]
        @test !low.active && high.active
        @test high.delta ≈ sqrt(inv(π*1e6*4π*1e-7*real(material.admittivity[2])))
        @test high.first_size <= high.delta/6
        @test high.extent <= min(5high.delta,.8region.shape.r)
    end
    snapshots = Dict()
    function snapshot(model, geometry)
        FEM._inspect_loaded_mesh(model,"conductor-controls")
        counts = [sum(length(block) for curve in curves
            for block in gmsh.model.mesh.get_elements(1,curve)[2])
            for curves in geometry.terminal_curves]
        @test counts == [96,96]
        # Interface geometry stays fixed across the frequency-dependent mesh.
        return counts
    end
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"conductor-controls",model.mesh_plans[2])
        for index in (2,1,2)
            FEM._update_exterior_mesh!(model,geometry,model.mesh_plans[index])
            FEM._configure_mesh!(model,geometry,model.mesh_plans[index])
            if index == 1
                all_surfaces = Set(last.(gmsh.model.get_entities(2)))
                for (layer,_) in values(geometry.conductor_fields)
                    @test Set(gmsh.model.mesh.field.get_numbers(layer,"ExcludedSurfacesList")) == all_surfaces
                end
            end
            gmsh.model.mesh.generate(2)
            result = snapshot(model,geometry)
            if haskey(snapshots,index)
                @test result == snapshots[index]
            end
            snapshots[index] = result
        end
    finally
        FEM._finish_gmsh(session)
    end
    # Compare exports with fresh construction. Reusing a model changes native
    # entity tags, so its unstructured triangulation need not be identical.
    for index in (1,2)
        fresh_session = FEM._start_gmsh(0)
        try
            geometry = FEM._build_geometry!(model,"fresh-controls",model.mesh_plans[index])
            FEM._configure_mesh!(model,geometry,model.mesh_plans[index])
            gmsh.model.mesh.generate(2)
            snapshots[index] = snapshot(model,geometry)
        finally
            FEM._finish_gmsh(fresh_session)
        end
    end
    mktempdir() do directory
        entry = export_data(:onelab,problem,formulation;
            file_name=joinpath(directory,"study.pro"),mesh_options=controls)
        session = FEM._start_gmsh(0)
        try
            for index in (1,2)
                gmsh.clear(); gmsh.onelab.clear(); gmsh.parser.clear()
                gmsh.parser.set_number("FrequencyIndex",[Float64(index)])
                gmsh.onelab.set_number("Inputs/01Frequency case",[Float64(index)])
                gmsh.open(replace(entry,r"\.pro$"=>".geo"))
                gmsh.model.mesh.generate(2)
                @test FEM._inspect_loaded_mesh(model,"native-conductor-controls") === nothing
                contours = [gmsh.model.get_entities_for_physical_group(1,
                    model.tags.terminal_contour_base+i) for i in 1:2]
                @test [sum(length(block) for c in curves for block in
                    gmsh.model.mesh.get_elements(1,c)[2]) for curves in contours] == snapshots[index]
                for (i,region) in enumerate(model.region_plans)
                    sizes = FEM._conductor_mesh_sizes(model,region,model.mesh_plans[index])
                    @test only(gmsh.parser.get_number("ConductorRegion$(i)First")) ≈ sizes.first_size
                    @test only(gmsh.parser.get_number("ConductorRegion$(i)Bulk")) ≈ sizes.bulk
                end
            end
            # Native frequency/material edits recompute spacing. Global mesh
            # coarsening must not erase the independently prescribed layers.
            edit = joinpath(directory,"edit.geo")
            write(edit,"""
                Include "study_data.pro";
                MeshSizeFactor = $(3controls.mesh_size_factor);
                FrequencyIndex = 2;
                Frequencies(1) = 4e6;
                MaterialSigma_1(1) = MaterialSigma_1(1)*4;
                Include "geometry/case-0002.geo";
                """)
            before = only(gmsh.parser.get_number("ConductorRegion1Delta"))
            gmsh.clear(); gmsh.onelab.clear(); gmsh.parser.clear()
            gmsh.open(edit)
            @test only(gmsh.parser.get_number("ConductorRegion1Delta")) ≈ before/4
            gmsh.model.mesh.generate(2)
            @test FEM._inspect_loaded_mesh(model,"native-edited-conductor-controls") === nothing
            @test only(gmsh.parser.get_number("ConductorCircleSegments")) == 96
        finally
            gmsh.onelab.clear()
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / screened cable voltage paths retain field edges" tags=[:extension] begin
    using Gmsh
    include(joinpath(@__DIR__,"..","manual","calculations","fem_conductor_mesh","fixtures.jl"))
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    design = ConductorMeshFixtures.screened_cable()
    system = build(LineCableSystem,design,Pose2(0.,-1.);
        connections=Dict(name=>i for (i,name) in enumerate(design.terminal_order)))
    problem = LineParametersProblem(system;frequencies=[1e6],
        earth_props=homogeneous(rho=100.,eps_r=1.))
    controls = (domain_skin_depths=24.,pml_layers=4,mesh_size_factor=3.,
        exterior_mesh_size_factor=8.,conductor_geometry_tolerance=2.5e-4,
        conductor_skin_depth_elements=12.,conductor_mesh_growth=1.25^0.25)
    model = FEM._resolved_fem_model(problem,LineCableModelsFEM(),
        computation_options(LineCableModelsFEM,ComputationOptions(;controls...)))
    foil = only(filter(r->r.terminal_index==3,model.region_plans))
    function foil_edge_size()
        curves = gmsh.model.get_entities_for_physical_group(1,model.tags.terminal_contour_base+3)
        return maximum(curves) do curve
            edges = reshape(only(gmsh.model.mesh.get_elements(1,curve)[3]),2,:)
            maximum(eachcol(edges)) do edge
                a,b = (gmsh.model.mesh.get_node(tag)[1] for tag in edge)
                hypot((a.-b)...)
            end
        end
    end
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"screen-paths")
        for reuse in (false,true)
            reuse && FEM._update_exterior_mesh!(model,geometry,only(model.mesh_plans))
            FEM._configure_mesh!(model,geometry)
            gmsh.model.mesh.generate(2)
            @test FEM._inspect_loaded_mesh(model,"screen-paths") === nothing
            @test foil_edge_size() <= foil.mesh_size*(1+1e-9)
        end
    finally
        FEM._finish_gmsh(session)
    end
    mktempdir() do directory
        entry = export_data(:onelab,problem,LineCableModelsFEM();
            file_name=joinpath(directory,"study.pro"),mesh_options=controls)
        session = FEM._start_gmsh(0)
        try
            gmsh.onelab.clear(); gmsh.parser.clear()
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            gmsh.model.mesh.generate(2)
            @test FEM._inspect_loaded_mesh(model,"native-screen-paths") === nothing
            @test foil_edge_size() <= foil.mesh_size*(1+1e-9)
            extension = only(filter(tag -> gmsh.model.mesh.field.get_type(tag)=="Extend",
                gmsh.model.mesh.field.list()))
            distance = gmsh.model.mesh.field.get_number(extension,"DistMax")
            write(joinpath(directory,"edit.geo"),"""
                Include "study_data.pro";
                MeshSizeFactor = $(0.5controls.mesh_size_factor);
                Include "geometry/case-0001.geo";
                """)
            gmsh.clear(); gmsh.onelab.clear(); gmsh.parser.clear()
            gmsh.open(joinpath(directory,"edit.geo"))
            @test gmsh.model.mesh.field.get_number(extension,"DistMax") ≈ distance/2
            gmsh.model.mesh.generate(1)
            @test foil_edge_size() <= foil.mesh_size/2*(1+1e-9)
        finally
            gmsh.onelab.clear(); gmsh.parser.clear()
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / merged arc endpoints do not create a full circle" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    session = FEM._start_gmsh(0)
    try
        Gmsh.gmsh.model.add("merged-arc")
        registry = FEM.FEMLoopRegistry(1e-4)
        centre,radius = (0.00212,1.0),0.00212
        # Independent contact calculations differ in angle, but the physical
        # endpoints are the same within the registry's Float64 tolerance.
        FEM._register_circle_break!(registry,centre,radius,π+1.36e-12)
        curves = FEM._circle_arc_path!(registry,centre,radius,2.6,1.0)
        @test length(curves) == 2
        @test all(c -> begin a,b=registry.curve_points[abs(c)]; a!=b end,curves)
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / sector partitions retain terminals and native controls" tags=[:extension] begin
    using Gmsh
    include(joinpath(@__DIR__,"..","manual","calculations","fem_conductor_mesh","fixtures.jl"))
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    design = ConductorMeshFixtures.sector_cable()
    system = build(LineCableSystem,design,Pose2(0.,1.);
        connections=Dict(name=>i for (i,name) in enumerate(design.terminal_order)))
    problem = LineParametersProblem(system;frequencies=[.1,1e4,1e6],
        earth_props=homogeneous(rho=100.,eps_r=1.))
    formulation = LineCableModelsFEM()
    controls = (pml_layers=4,mesh_size_factor=3.,exterior_mesh_size_factor=8.)
    model = FEM._resolved_fem_model(problem,formulation,
        computation_options(LineCableModelsFEM,ComputationOptions(;controls...)))
    counts = Dict{Int,Vector{Int}}()
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"sector-controls")
        @test length(geometry.sector_partitions) == 3
        for k in (3,1,2,3)
            FEM._update_exterior_mesh!(model,geometry,model.mesh_plans[k])
            FEM._configure_mesh!(model,geometry,model.mesh_plans[k])
            gmsh.model.mesh.generate(2)
            @test FEM._inspect_loaded_mesh(model,"sector-controls") === nothing
            for (index,partition) in geometry.sector_partitions
                terminal = model.region_plans[index].terminal_index
                contour = Set(geometry.terminal_curves[terminal])
                @test contour == Set(first.(partition.curve_pairs))
                @test isempty(intersect(contour,abs.(partition.spokes)))
                @test isempty(intersect(contour,getindex.(partition.curve_pairs,2)))
                @test Set(geometry.terminal_surfaces[terminal]) ==
                    Set([collect(keys(partition.patches));partition.core])
                # Every strip retains the prescribed triangular arrangement.
                @test all(s -> haskey(geometry.transfinite_surfaces,s),keys(partition.patches))
            end
            value = [geometry.transfinite_curves[abs(first(p.spokes))][1]
                for (_,p) in sort!(collect(geometry.sector_partitions);by=first)]
            haskey(counts,k) && @test value == counts[k]
            counts[k] = value
        end
        @test all(counts[3] .> counts[1])
    finally
        FEM._finish_gmsh(session)
    end
    mktempdir() do directory
        entry = export_data(:onelab,problem,formulation;
            file_name=joinpath(directory,"study.pro"),mesh_options=controls)
        session = FEM._start_gmsh(0)
        try
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            gmsh.model.mesh.generate(2)
            @test FEM._inspect_loaded_mesh(model,"native-sector-controls") === nothing
            @test only(gmsh.parser.get_number("ConductorRegion1Segments")) == 192
            @test only(gmsh.parser.get_number("ConductorRegion1Layers"))+1 == counts[1][1]
            edit = joinpath(directory,"edit.geo")
            write(edit,"""
                Include "study_data.pro";
                FrequencyIndex = 3;
                Frequencies(2) = 4e6;
                MaterialSigma_1(2) = MaterialSigma_1(2)*4;
                ConductorMeshGrowth = 1;
                Include "geometry/case-0003.geo";
                """)
            high = FEM._conductor_mesh_sizes(model,first(model.region_plans),last(model.mesh_plans))
            gmsh.clear(); gmsh.onelab.clear(); gmsh.parser.clear()
            gmsh.open(edit)
            @test only(gmsh.parser.get_number("ConductorRegion1Delta")) ≈ high.delta/4
            # Parsing native edits exercises the uniform-spacing branch
            # without launching a deliberately very fine global mesh.
            @test only(gmsh.parser.get_number("ConductorRegion1Layers")) > counts[3][1]
        finally
            gmsh.onelab.clear()
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / annular grading covers both faces of a thin metal wall" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor,rho=2e-8)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    cable = build(CableDesign,"thin-wall",Stack(
        terminal(:core,core(copper;r=.0015)),insulation(dielectric;t=.001),
        terminal(:wall,sheath(copper;t=.00015)),jacket(dielectric;t=.0005)))
    system = build(LineCableSystem,cable,Pose2(0.,-.1);
        connections=Dict(:core=>1,:wall=>2))
    problem = LineParametersProblem(system;frequencies=[1e6],
        earth_props=homogeneous(rho=100.,eps_r=1.))
    options = computation_options(LineCableModelsFEM,ComputationOptions(
        pml_layers=4,conductor_thickness_elements=8,conductor_mesh_growth=1.))
    model = FEM._resolved_fem_model(problem,LineCableModelsFEM(),options)
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"thin-wall")
        FEM._configure_mesh!(model,geometry)
        gmsh.model.mesh.generate(2)
        @test FEM._inspect_loaded_mesh(model,"thin-wall") === nothing
        @test length(geometry.conductor_fields) == 2
        counts = [sum(length(block) for c in curves for block in
            gmsh.model.mesh.get_elements(1,c)[2]) for curves in geometry.terminal_curves]
        @test counts[1] == 96
        @test counts[2] > 192 # The local wall target is finer than the angular bound.
        wall = only(filter(r->r.terminal_index==2,model.region_plans))
        sizes = FEM._conductor_mesh_sizes(model,wall,only(model.mesh_plans))
        @test sizes.active
        @test sizes.first_size <= min(sizes.delta/6,.00015/8)*(1+1e-12)
        @test sizes.bulk ≈ .00015/8
        @test sizes.extent ≈ .45*.00015
    finally
        FEM._finish_gmsh(session)
    end
end
