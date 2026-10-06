@testitem "FEM / intent options and declared native overrides" tags=[:extension] begin
    using Gmsh
    E=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    controls(;kwargs...)=computation_options(E.LineCableModelsFEM,ComputationOptions(;kwargs...)).data
    @test (@inferred computation_options(E.LineCableModelsFEM,ComputationOptions())).data.pml_element_family==1
    @test (@inferred computation_options(E.LineCableModelsFEM,ComputationOptions(overrides=(MumpsOrdering=-1,PetscPrealloc=256)))).data.mumps_ordering==-1
    @test controls(overrides=(PmlSideLayers=48,PmlTopLayers=48,PmlBottomLayers=48,)).pml_layers==(48,48,48)
    @test controls(overrides=(PmlTopLayers=24,)).pml_layers==(16,24,16)
    @test controls(overrides=(PmlSideLayers=8,PmlTopLayers=6,PmlBottomLayers=4,)).pml_layers==(8,6,4)
    @test_throws ArgumentError controls(overrides=(NotANativeParameter=1,))
    @test_throws ArgumentError controls(overrides=Dict(:PmlSideLayers=>48))
    @test_throws ArgumentError controls(mesh_path="mesh.msh",mesh_policy=:remesh)
    @test controls(overrides=(PmlQuadrangles=0,)).pml_element_family==0
    @test controls(overrides=(PetscPrealloc=0,MumpsOrdering=-1)).petsc_prealloc==0
    @test controls(overrides=(MumpsOrdering=-1,)).mumps_ordering==-1
    @test !haskey(controls(),:mumps_forward_error_tolerance)
    for spec in E.FEM_NATIVE_OVERRIDES
        for name in spec.names
            @test getproperty(controls(overrides=(;name=>spec.default)),spec.field)==getproperty(controls(),spec.field)
            for bad in ("1; Error(\"injection\");",true,NaN,Inf)
                @test_throws ArgumentError controls(overrides=(;name=>bad))
            end
        end
        err=try controls(;spec.field=>spec.default);catch e;e;end
        @test err isa ArgumentError
        @test occursin(string(first(spec.names)),sprint(showerror,err))
    end
    @test_throws ArgumentError controls(overrides=(PhysicalVolumeQuadrature=1,))
    @test_throws ArgumentError controls(overrides=(PmlSideLayers=0,PmlTopLayers=0,PmlBottomLayers=0,))
    @test_throws ArgumentError controls(mumps_forward_error_tolerance=.02)
end

@testitem "FEM / one resolution scale coarsens conductor meshes" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    records=[]
    for scale in (1.,2.)
        N.geometry(N.problem(;frequency=1e6,rho=100.,radius=.01);options=(mesh_size_factor=scale,),mesh=true) do g,log
            regions=g.model.get_entities_for_physical_group(2,10001)
            count=sum(length(g.model.mesh.get_nodes(2,region,true)[1]) for region in regions)
            push!(records,(count,only(g.parser.get_number("FEMCircleSegments")),only(g.parser.get_number("FEMFirst")),only(g.parser.get_number("FEMBulk")),only(g.parser.get_number("FEMPmlPPW_0"))))
        end
    end
    @test records[2][1]<records[1][1]
    @test records[2][2]==records[1][2]
    @test records[2][3]>records[1][3]
    @test records[2][4]==2records[1][4]
    @test records[2][5]==records[1][5]/2
end

@testitem "FEM / unified detached options and guarded native inputs" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh
    N=NativeFEMFixtures;E=N.FEM;p=N.problem()
    options=(plot_field_maps=true,solver_threads=2,output_basis=:total,
        linear_solver=:gmres,overrides=(PmlSideLayers=48,PmlTopLayers=48,PmlBottomLayers=48,MumpsOrdering=0))
    N.bundle(p;options) do root,entry
        data=read(joinpath(root,"model_data.pro"),String)
        @test occursin("PlotFieldMaps = {1,",data)
        @test occursin("GetDPThreads = {2,",data)
        @test occursin("OutputTotal = {1,",data)
        @test occursin("LinearSolver = {1,",data)
        @test occursin("FEMMumpsForwardErrorBudget = $(E._pro_number(E.FEM_MUMPS_FORWARD_ERROR_BUDGET));",data)
        for spec in E.FEM_NATIVE_OVERRIDES,name in spec.names
            @test occursin("Expert/$name",data)
            if spec.closed === false
                line=only(filter(l -> startswith(strip(l),"$name ="),split(data,'\n')))
                @test parse(Float64,match(r"Min ([^,]+)",line)[1])>first(spec.range)
            end
        end
    end
    for name in (:frequency_workers,:mesh_policy,:mesh_path,:keep_run_directory,:resume_run_directory,:on_result,:log_file,:trace,:timing,:verbosity,:getdp_executable)
        mktempdir() do root
            @test_throws ArgumentError export_data(:onelab,p,Formulation(:LineCableModelsFEM);file_name=joinpath(root,"bad","model.pro"),options=(;name=>nothing))
            @test !isdir(joinpath(root,"bad"))
        end
    end
    N.bundle(p) do root,entry
        getdp=E._getdp_selection(computation_options(E.LineCableModelsFEM,ComputationOptions())).path
        for (name,value,message) in (("PhysicalVolumeQuadrature",1,"PhysicalVolumeQuadrature"),("PmlTopLayers",0,"PML intervals"),("PmlReflection",1,"PML reflection"),("MeshSizeFactor",0,"physical mesh controls"),("PmlBottomThicknessFactor",0,"Native mesh factors"))
            probe=joinpath(root,"invalid.pro")
            write(probe,"Include \"model_data.pro\";\n$name = $value;\nInclude \"formulations/parameters.pro\";\n")
            Gmsh.gmsh.initialize(String[],false,false)
            try
                Gmsh.gmsh.option.set_number("General.Terminal",0)
                @test_throws Exception Gmsh.gmsh.parser.parse(probe)
            finally
                Gmsh.finalize()
            end
            output=IOBuffer();proc=run(pipeline(ignorestatus(`$getdp $probe -v 3`),stdout=output,stderr=output))
            @test !success(proc)
            @test occursin(message,String(take!(output)))
        end
    end
end
