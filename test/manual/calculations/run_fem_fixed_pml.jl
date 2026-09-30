# Bounded observation only: three prescribed meshes, two frequencies, quasi-fw.
# Run from the repository root. No reference tolerance or setting selection.
# Restarting uses the backend's ordinary compatible-run resume mechanism.
using LineCableModels, Gmsh, Printf, TOML, JSON3

function fixed_pml_observation(root)
    mkpath(root)
    copper = Material(MaterialsLibrary(add_defaults=true),:copper)
    wire = build(CableDesign,"two_bare_wires",Stack(Group(:core,
        Region(:core_metal,Disk(.0425),copper))))
    system = build(LineCableSystem,[wire,wire],[Pose2(0.,1.),Pose2(1.,-1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)],system_id="two_bare_wires",line_length=1.)
    frequencies = [.1,1e6]
    problem = LineParametersProblem(system;frequencies,temperature=20.,
        earth_props=homogeneous(rho=100.,eps_r=1.,mu_r=1.))
    formulation = Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    cases = (("192-192-192",(192,192,192)),("96-96-96",(96,96,96)),("96-96-48",(96,96,48)))
    controls = (domain_skin_depths=24.,pml_grading=(192/191)*log(1536),
        mesh_size_factor=3.,exterior_mesh_size_factor=8.,
        conductor_geometry_tolerance=1e-3,conductor_skin_depth_elements=6.,
        conductor_mesh_growth=sqrt(1.25),conductor_skin_depths=5.,conductor_thickness_elements=4,
        frequency_workers=4,solver_threads=1,mesh_policy=:reuse,resume_run_directory=:latest,
        keep_run_directory=true,trace=true,timing=true,plot_field_maps=false,
        gmsh_verbosity=2,getdp_verbosity=4,verbosity=(default=1,))
    results = Any[]
    costs = Dict[]
    mesh_rows = NamedTuple[]
    for (label,layers) in cases
        directory = joinpath(root,label)
        mkpath(directory)
        println("START ",label," / 0.1 and 1e6 Hz"); flush(stdout)
        options = (;controls...,pml_layers=layers)
        write(joinpath(directory,"options.txt"),repr(options)*"\n")
        measured = @timed compute(problem,formulation;options)
        result = measured.value
        record = details(result).data.fem.run
        run_directory = record.run_directory
        cost_path = joinpath(directory,"cost.toml")
        previous = isfile(cost_path) ? TOML.parsefile(cost_path) : Dict()
        cost = if record.reused && get(previous,"run_directory",nothing)==run_directory
            previous
        else
            Dict("wall_seconds"=>measured.time,"compile_seconds"=>measured.compile_time,
                "recompile_seconds"=>measured.recompile_time,"gc_seconds"=>measured.gctime,
                "run_directory"=>run_directory,"reused"=>record.reused,
                "recovered_columns"=>record.recovered_columns)
        end
        open(io->TOML.print(io,cost),cost_path,"w")
        push!(costs,cost)
        push!(results,result)
        open(joinpath(directory,"matrices.csv"),"w") do io
            println(io,"quantity,frequency_hz,receiver,source,real,imaginary")
            for (name,matrix) in (("Z",Z(result)),("Y",Y(result)),
                    ("P",details(result).data.fem.primitive.P_primitive))
                for k in eachindex(frequencies), j in 1:2, i in 1:2
                    @printf(io,"%s,%.17g,%d,%d,%.17g,%.17g\n",name,frequencies[k],i,j,
                        real(matrix[i,j,k]),imag(matrix[i,j,k]))
                end
            end
        end
        # Inspect saved meshes after the independent public computation.
        gmsh = Gmsh.gmsh
        gmsh.initialize(String[],false,false)
        try
            gmsh.option.set_number("General.Verbosity",0)
            for (index,frequency) in enumerate(frequencies)
                stem = index==length(frequencies) ? "model" : @sprintf("frequency_%04d",index)
                path = joinpath(run_directory,"mesh",stem*".msh")
                metadata = JSON3.read(read(joinpath(run_directory,"mesh",stem*".json"),String))
                gmsh.open(path)
                nodes = length(first(gmsh.model.mesh.get_nodes()))
                triangles = sum(length,gmsh.model.mesh.get_elements(2)[2])
                pml_tag = only(g.tag for g in metadata.physical_groups if g.name=="LCM/domain/pml")
                pml,corners = 0,0
                for surface in gmsh.model.get_entities_for_physical_group(2,pml_tag)
                    count = sum(length,gmsh.model.mesh.get_elements(2,surface)[2])
                    pml += count
                    xmin,ymin,_,xmax,ymax,_ = gmsh.model.get_bounding_box(2,surface)
                    # Physical box is centred at x=0.5 for this fixed fixture.
                    if abs((xmin+xmax)/2-.5)>metadata.physical_domain_halfwidth_m &&
                            abs((ymin+ymax)/2)>metadata.physical_domain_halfwidth_m
                        corners += count
                    end
                end
                push!(mesh_rows,(;label,frequency,nodes,triangles,pml,corners))
                gmsh.clear()
            end
        finally
            gmsh.finalize()
        end
        println("DONE ",label," wall=",cost["wall_seconds"]," s, compile=",cost["compile_seconds"]," s; ",run_directory)
        flush(stdout)
    end
    open(joinpath(root,"meshes.csv"),"w") do io
        println(io,"case,frequency_hz,nodes,triangles,pml_triangles,corner_triangles")
        for row in mesh_rows
            println(io,join(values(row),','))
        end
    end
    # All comparisons follow the three independent computations. Values are raw;
    # P is the complete FEM inverse-admittance coefficient (ohm m), Y=inv(P).
    quantities(result) = (("R",real.(Z(result))),("X",imag.(Z(result))),
        ("G",real.(Y(result))),("B",imag.(Y(result))),
        ("P_real",real.(details(result).data.fem.primitive.P_primitive)),
        ("P_imag",imag.(details(result).data.fem.primitive.P_primitive)))
    reference = quantities(first(results))
    open(joinpath(root,"differences.csv"),"w") do io
        println(io,"case,quantity,frequency_hz,receiver,source,baseline,value,signed_difference,absolute_difference")
        for c in 2:3, (q,(name,matrix)) in enumerate(quantities(results[c])), k in 1:2, j in 1:2, i in 1:2
            base = reference[q][2][i,j,k]
            value = matrix[i,j,k]
            @printf(io,"%s,%s,%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g\n",
                cases[c][1],name,frequencies[k],i,j,base,value,value-base,abs(value-base))
        end
    end
    open(joinpath(root,"observations.md"),"w") do io
        println(io,"# Fixed PML controls: bounded observations\n")
        println(io,"Mixed bare copper wires, radius 0.0425 m, centres (0,+1) and (1,-1) m; earth 100 ohm m, relative permittivity and permeability one. Quasi-fw, 0.1 Hz and 1 MHz. Grading g=(192/191)log(1536) in every direction. Six frequency solves, twelve source columns; no additional scientific qualification.\n")
        println(io,"Physical and conductor controls match the active manual preset. Four frequency workers (only two available frequencies), one native thread. Counts below are measured, not projected.\n")
        println(io,"| Counts | Scan wall s | Compilation s (included) | Reused / recovered columns |\n|---|---:|---:|---|")
        for (c,(label,_)) in enumerate(cases)
            cost=costs[c]
            @printf(io,"| %s | %.3f | %.3f | %s / %d |\n",label,cost["wall_seconds"],cost["compile_seconds"],cost["reused"],cost["recovered_columns"])
        end
        println(io,"\nThe first fresh call includes compilation; later fresh calls use warmed Julia code. This is a single observation per configuration, not a repeated steady-state timing study. Native stage timings remain in each run's timing-summary.json.\n")
        println(io,"| Counts | Hz | Nodes | Triangles | PML triangles | Corner triangles |\n|---|---:|---:|---:|---:|---:|")
        for r in mesh_rows
            println(io,"| ",join(values(r)," | ")," |")
        end
        println(io,"\nRaw signed R/X/G/B and both parts of complete P appear in matrices.csv and differences.csv; no clipping or sign correction is applied. Maximum absolute differences below summarize the four entries, without mixing real and imaginary components.\n")
        println(io,"| Counts | Hz | Component | Maximum absolute difference |\n|---|---:|---|---:|")
        for c in 2:3, (q,(name,matrix)) in enumerate(quantities(results[c])), k in 1:2
            @printf(io,"| %s | %.6g | %s | %.9g |\n",cases[c][1],frequencies[k],name,maximum(abs.(matrix[:,:,k]-reference[q][2][:,:,k])))
        end
        println(io,"\nConductance matrices in S/m (rows are receivers; columns are sources):\n")
        for c in 1:3, k in 1:2
            matrix=real.(Y(results[c]))[:,:,k]
            @printf(io,"%s at %.6g Hz: `[%.17g %.17g; %.17g %.17g]`\n\n",cases[c][1],frequencies[k],matrix[1,1],matrix[1,2],matrix[2,1],matrix[2,2])
        end
        println(io,"These mesh-to-mesh differences are observations, not error bounds or scientific acceptance. No discretization bound on Z, P or Y is available here. The user chooses the operating controls; the active preset remains unchanged.")
    end
    println("COMPLETE six prescribed frequency solves / twelve source columns; ",joinpath(root,"observations.md")); flush(stdout)
end

if abspath(PROGRAM_FILE)==@__FILE__
    fixed_pml_observation(isempty(ARGS) ? joinpath(pwd(),".linecablemodels/fem/fixed-pml-controls/observation") : abspath(ARGS[1]))
end
