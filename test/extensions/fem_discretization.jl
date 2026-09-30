@testitem "Gmsh FEM / prescribed discretization managed and detached" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"discretization",terminal(:core,
        core(Material(kind=:conductor,rho=1.72e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[50.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    form = Formulation(:LineCableModelsFEM;options=(Γ=.01+.02im,
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    function inventory(model)
        pml = Set(gmsh.model.get_entities_for_physical_group(2,model.tags.pml))
        surfaces = [s for (_,s) in gmsh.model.get_entities(2)]
        kinds = Dict(s=>first(gmsh.model.mesh.get_elements(2,s)) for s in surfaces)
        counts = Dict(s=>sum(length,gmsh.model.mesh.get_elements(2,s)[2]) for s in surfaces)
        lines = Dict(c=>gmsh.model.mesh.get_elements(1,c)[3] for (_,c) in gmsh.model.get_entities(1))
        physical = Dict(s=>gmsh.model.mesh.get_elements(2,s)[3] for s in setdiff(surfaces,pml))
        (;pml,kinds,counts,lines,physical)
    end
    mktempdir() do root
        inventories = []
        fingerprints = String[]
        for (family,physical,quad) in ((:triangle,nothing,9),(:triangle,3,9),
                (:quadrangle,3,9),(:quadrangle,12,16))
            mesh_controls = (pml_layers=8,mesh_size_factor=3.,
                pml_element_family=family,physical_volume_quadrature=physical,pml_quadrature=quad)
            options = computation_options(LineCableModelsFEM,ComputationOptions(;
                mesh_controls...,keep_run_directory=true,mesh_policy=:remesh,
                frequency_workers=1,solver_threads=1,gmsh_verbosity=0,getdp_verbosity=2))
            model = FEM._resolved_fem_model(problem,form,options)
            plan = only(model.mesh_plans)
            push!(fingerprints,FEM._mesh_fingerprint(model,"test",plan))
            session = FEM._start_gmsh(0)
            mesh = joinpath(root,"$(family)-$(physical)-$(quad).msh")
            try
                geometry = FEM._build_geometry!(model,"discretization",plan)
                FEM._configure_mesh!(model,geometry,plan)
                gmsh.model.mesh.generate(2)
                gmsh.write(mesh)
                @test FEM._inspect_loaded_mesh(model,mesh) === nothing
                inv = inventory(model)
                @test all(inv.kinds[s] == [family===:triangle ? 2 : 3] for s in inv.pml)
                @test all(inv.kinds[s] == [2] for s in keys(inv.physical))
                push!(inventories,inv)
            finally
                FEM._finish_gmsh(session)
            end
            # Reuse this exact mesh for the managed/native comparison: do not
            # turn unrelated remeshing roundoff into an execution-parity test.
            result = compute(problem,form;options=ComputationOptions(;
                mesh_controls...,mesh_policy=:reuse,mesh_path=mesh,
                keep_run_directory=true,frequency_workers=1,solver_threads=1,
                gmsh_verbosity=0,getdp_verbosity=2))
            entry = export_data(:onelab,problem,form;
                file_name=joinpath(root,"bundle-$(family)-$(physical)-$(quad)","study.pro"),
                mesh_options=mesh_controls)
            cmd = `$getdp $entry -msh $mesh -solve LineCableModelsFEM -setnumber PlotFieldMaps 0 -v 2`
            run(pipeline(addenv(cmd,"OPENBLAS_NUM_THREADS"=>"1","OMP_NUM_THREADS"=>"1");stdout=devnull))
            for quantity in (Z,Y)
                native = zeros(ComplexF64,2,2)
                table = joinpath(dirname(entry),"results/f0001-quasi-fw-b0000/matrices/$quantity.tsv")
                for row in split.(readlines(table)[3:end],'\t')
                    native[parse(Int,row[1]),parse(Int,row[2])] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
                end
                @test native ≈ observe(result,quantity)[:,:,1] rtol=1e-10
            end
            data = read(joinpath(dirname(entry),"study_data.pro"),String)
            @test occursin("PmlQuadrangles = {$(Int(family===:quadrangle)),",data)
            @test occursin("PhysicalVolumeQuadrature = {$(something(physical,0)),",data)
            # Detached geometry must construct the requested cells as well.
            session = FEM._start_gmsh(0)
            try
                gmsh.open(replace(entry,r"\.pro$"=>".geo"))
                gmsh.model.mesh.generate(2)
                detached = inventory(model)
                @test detached.kinds == inventories[end].kinds
                @test all(detached.counts[s] == inventories[end].counts[s] for s in detached.pml)
                @test detached.lines == inventories[end].lines
            finally
                gmsh.onelab.clear(); gmsh.parser.clear()
                FEM._finish_gmsh(session)
            end
        end
        @test fingerprints[1] == fingerprints[2]
        @test fingerprints[2] != fingerprints[3]
        @test fingerprints[3] == fingerprints[4]
        @test inventories[1].physical == inventories[3].physical
        @test inventories[1].lines == inventories[3].lines
        @test all(inventories[1].counts[s] == 2inventories[3].counts[s] for s in inventories[1].pml)
    end
end
