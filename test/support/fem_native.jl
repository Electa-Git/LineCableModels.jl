@testmodule NativeFEMFixtures begin
    using LineCableModels, Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const gmsh = Gmsh.gmsh

    function problem(; frequency=1e6, rho=1000., eps_r=12., radius=.01,
            positions=[(0.,1.)], radii=fill(radius,length(positions)))
        copper=Material(kind=:conductor,rho=1.72e-8,eps_r=1.,mu_r=1.)
        designs=[build(CableDesign,"wire-$i",
            terminal(:core,Region(:metal,Disk(r),copper))) for (i,r) in pairs(radii)]
        system=build(LineCableSystem,designs,Pose2.(first.(positions),last.(positions));
            connections=[Dict(:core=>i) for i in eachindex(designs)],
            system_id="native-mesh",line_length=1.)
        LineParametersProblem(system;temperature=20.,frequencies=[frequency],
            earth_props=homogeneous(;rho,eps_r,mu_r=1.))
    end

    function bundle(f,problem; gamma=0., options=(;))
        mktempdir() do directory
            form=Formulation(:LineCableModelsFEM;options=(Γ=gamma,
                reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
            entry=export_data(:onelab,problem,form;
                file_name=joinpath(directory,"model.pro"),mesh_options=options)
            f(directory,entry)
        end
    end

    function parameters(f,problem;gamma=0.,options=(;),expressions)
        bundle(problem;gamma,options) do directory,entry
            probe=joinpath(directory,"probe.pro")
            open(probe,"w") do io
                println(io,"Include \"model_data.pro\";")
                println(io,"Include \"formulations/parameters.pro\";")
                for (i,(name,expression)) in enumerate(expressions)
                    println(io,"FEMProbe$i = $expression;")
                    println(io,"Printf(\"NativeProbe$i %.17g\",FEMProbe$i);")
                end
            end
            gmsh.initialize(String[],false,false)
            values=try
                gmsh.option.set_number("General.Terminal",0)
                gmsh.parser.parse(probe)
                Dict(name=>only(gmsh.parser.get_number("FEMProbe$i"))
                    for (i,(name,_)) in enumerate(expressions))
            finally
                Gmsh.finalize()
            end
            selection=FEM._getdp_selection(FEM.computation_options(FEM.LineCableModelsFEM,ComputationOptions()))
            output=read(pipeline(`$(selection.path) $probe -v 3`;stderr=stdout),String)
            native=Dict(name=>parse(Float64,only(match(Regex("NativeProbe$i ([^\\r\\n]+)"),output).captures))
                for (i,(name,_)) in enumerate(expressions))
            f(values,native)
        end
    end

    function geometry(f,problem;gamma=0.,options=(;),mesh=false)
        bundle(problem;gamma,options) do directory,entry
            gmsh.initialize(String[],false,false)
            try
                gmsh.option.set_number("General.Terminal",0)
                gmsh.logger.start()
                gmsh.open(replace(entry,r"\.pro$"=>".geo"))
                mesh && gmsh.model.mesh.generate(2)
                f(gmsh,copy(gmsh.logger.get()))
            finally
                gmsh.logger.stop()
                Gmsh.finalize()
            end
        end
    end
end
