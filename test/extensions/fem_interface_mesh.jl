@testitem "Gmsh FEM / cable-local interface grading and native controls" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    wire = build(CableDesign,"interface-control",terminal(:core,
        core(Material(kind=:conductor,rho=2e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,1.),(.5,-1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    f = 1e3
    problem = LineParametersProblem(system;frequencies=[f],earth_props=homogeneous(rho=100.,eps_r=1.))
    controls = (domain_skin_depths=8.,pml_layers=4,mesh_size_factor=3.,exterior_mesh_size_factor=8.)
    gamma = .99sqrt(im*2pi*f*4pi*1e-7*(.01+im*2pi*f*8.8541878128e-12))
    function inventory(model,expected_full_sources)
        interface = Set(gmsh.model.get_entities_for_physical_group(1,model.tags.interface))
        full_sources = count(gmsh.model.mesh.field.list()) do field
            gmsh.model.mesh.field.get_type(field)=="Distance" &&
                !isempty(intersect(interface,gmsh.model.mesh.field.get_numbers(field,"CurvesList")))
        end
        @test full_sources == expected_full_sources
        FEM._inspect_loaded_mesh(model,"interface-control")
        groups(tag) = gmsh.model.get_entities_for_physical_group(1,tag)
        lines = Dict(c=>length(first(gmsh.model.mesh.get_nodes(1,c,true)))
            for i in 1:2 for c in [groups(4000+i);groups(7000+i)])
        pml = sum(length(block) for s in gmsh.model.get_entities_for_physical_group(2,1005)
            for block in gmsh.model.mesh.get_elements(2,s)[2])
        return (;nodes=length(first(gmsh.model.mesh.get_nodes())),lines,pml)
    end
    mktempdir() do dir
        for (index,Γ) in enumerate((0.,gamma))
            form = Formulation(:LineCableModelsFEM;options=(;Γ))
            entry = export_data(:onelab,problem,form;
                file_name=joinpath(dir,"g$index","study.pro"),mesh_options=controls)
            before = read(joinpath(dirname(entry),"study_data.pro"),String)
            records = []
            identities = String[]
            for factor in (1.,1000.)
                model = FEM._resolved_fem_model(problem,form,
                    computation_options(LineCableModelsFEM,ComputationOptions(;
                        controls...,interface_refinement_factor=factor)))
                push!(identities,FEM._mesh_fingerprint(model,"test"))
                session = FEM._start_gmsh(0)
                try
                    geometry = FEM._build_geometry!(model,"interface-control")
                    FEM._configure_mesh!(model,geometry)
                    gmsh.model.mesh.generate(2)
                    managed = inventory(model,0)
                    gmsh.clear(); gmsh.parser.clear(); gmsh.onelab.clear()
                    gmsh.onelab.set_number("Mesh/09Interface footprint factor",[factor])
                    gmsh.open(replace(entry,r"\.pro$"=>".geo"))
                    @test only(gmsh.parser.get_number("InterfaceRefinementFactor")) == factor
                    gmsh.model.mesh.generate(2)
                    native = inventory(model,0)
                    @test native.lines == managed.lines
                    @test native.pml == managed.pml
                    push!(records,(;managed,native))
                finally
                    gmsh.onelab.clear(); gmsh.parser.clear()
                    FEM._finish_gmsh(session)
                end
            end
            @test identities[1] != identities[2]
            @test records[2].managed.nodes > records[1].managed.nodes
            @test records[2].native.nodes > records[1].native.nodes
            @test records[1].managed.lines == records[2].managed.lines
            @test records[1].managed.pml == records[2].managed.pml
            @test read(joinpath(dirname(entry),"study_data.pro"),String) == before
        end
    end
end
