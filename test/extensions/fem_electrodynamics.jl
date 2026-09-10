@testitem "Gmsh FEM / electrodynamic earth matrices and independent voltage basis" tags=[:extension, :integration, :fem_numerical] begin
    using Gmsh
    using LineCableModels
    using LinearAlgebra
    using Printf
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)

    function bare_problem(rho)
        metal = Material(kind=:conductor, rho=rho)
        design = build(CableDesign, "electric-regression",
            terminal(:core, core(metal; r=0.0425)))
        system = build(LineCableSystem, [design, design], [(0.0,-1.0),(1.0,-1.0)];
            connections=[Dict(:core=>1),Dict(:core=>2)], line_length=1.0)
        LineParametersProblem(system; frequencies=[1e5,1e6],
            earth_props=homogeneous(rho=0.1, eps_r=1.0, mu_r=1.0))
    end
    problem = bare_problem(1/5.8e7)
    options = (reduce_bundle=false, kron_reduction=false, ideal_transposition=false)
    fem = Formulation(:LineCableModelsFEM; options,
        fem_options=(gmsh_verbosity=0, getdp_verbosity=0, keep_run_directory=true))
    actual = compute(problem, fem; options=(trace=true,))
    proposed = compute(problem, Formulation(; options))
    xue = compute(problem, Formulation(; earth_impedance=:Xue2018,
        earth_admittance=:Xue2018, options))
    for index in eachindex(problem.frequencies)
        @test Y(actual)[:,:,index] ≈ transpose(Y(actual)[:,:,index]) rtol=1e-11
        @test Z(actual)[:,:,index] ≈ transpose(Z(actual)[:,:,index]) rtol=1e-9
        # Check each complex entry: a matrix norm can hide mutual-Y errors.
        for i in 1:2, j in 1:2
            @test Z(actual)[i,j,index] ≈ Z(proposed)[i,j,index] rtol=0.025
            @test Y(actual)[i,j,index] ≈ Y(proposed)[i,j,index] rtol=0.025
        end
        @test Y(actual)[2,1,index] ≈ Y(xue)[2,1,index] rtol=0.025
        @test minimum(eigvals(Symmetric(real.(Y(actual)[:,:,index])))) > 0
    end
    @test actual.details.fem.run.getdp_invocations == 2
    @test actual.details.fem.run.completed_columns == 4
    @test maximum(actual.details.fem.inversion_residuals) < 1e-12

    # The electric domain excludes both metal interiors. Making the metal
    # nearly perfect must not change its equipotential terminal admittance.
    pec = compute(bare_problem(1e-12), fem)
    @test Y(pec) ≈ Y(actual) rtol=1e-10

    run_directory = actual.details.fem.run.run_directory
    for frequency in 1:2, basis in 1:2
        timing = joinpath(run_directory, "raw/jobs",
            @sprintf("getdp-f%04d-b%04d-timing.tsv",frequency,basis))
        @test parse(Int,last(split(read(timing,String)))) == (basis == 1 ? 1 : 0)
    end

    # Independently prescribe 1 V / 0 V and extract the associated terminal
    # currents. Their matrix must equal inv(P) from the production current basis.
    mktempdir() do directory
        for name in readdir(joinpath(run_directory,"input/getdp"))
            cp(joinpath(run_directory,"input/getdp",name),joinpath(directory,name))
        end
        path = joinpath(directory,"quasi_tem.pro")
        source = read(path,String)
        source = replace(source,
            "Value \$FEM_Q~{t};" => "Value -\$FEM_Q~{t};",
            "NameOfCoef Q; EntityType Auto; NameOfConstraint FEMTransverseCurrent;" =>
                "NameOfCoef V; EntityType Auto; NameOfConstraint FEMTransverseCurrent;",
            "Re[{V} / UnitTransverseSource]" => "-Re[{Q} / UnitTransverseSource]",
            "Im[{V} / UnitTransverseSource]" => "-Im[{Q} / UnitTransverseSource]")
        write(path,source)
        command = FEM._getdp_command(actual.details.fem.inputs.getdp_provenance.path,
            joinpath(directory,"model.pro"),joinpath(run_directory,"mesh/model.msh"),
            (path=run_directory,),fem,last(actual.details.fem.inputs.mesh_plans),
            [1,2],directory; reuse_factorization=false)
        open(joinpath(directory,"getdp.log"),"w") do log
            run(pipeline(command; stdout=log, stderr=log))
        end
        voltage_y = zeros(ComplexF64,2,2)
        for basis in 1:2
            path = joinpath(directory,"raw/jobs",
                @sprintf("getdp-f0002-b%04d-P.tsv",basis))
            @test FEM._valid_job_raw(path,2,2,1e6,basis)
            for line in eachline(path)
                row = split(line)
                voltage_y[parse(Int,row[3]),basis] = complex(
                    parse(Float64,row[5]),parse(Float64,row[6]))
            end
        end
        @test voltage_y ≈ Y(actual)[:,:,2] rtol=1e-10
    end
    rm(run_directory;recursive=true)
    rm(pec.details.fem.run.run_directory;recursive=true)
end

@testitem "Gmsh FEM / electrodynamic overhead and mixed terminal reciprocity" tags=[:extension, :integration, :fem_numerical] begin
    using Gmsh
    using LineCableModels
    copper = Material(kind=:conductor,rho=1/5.8e7)
    design = build(CableDesign,"electric-media",terminal(:core,core(copper;r=0.01)))
    formulation = Formulation(:LineCableModelsFEM;
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false),
        fem_options=(gmsh_verbosity=0,getdp_verbosity=0,
            plot_field_maps=true,keep_run_directory=true))
    for positions in ([(0.0,1.0),(1.0,1.0)],[(0.0,1.0),(1.0,-1.0)])
        system = build(LineCableSystem,[design,design],positions;
            connections=[Dict(:core=>1),Dict(:core=>2)])
        problem = LineParametersProblem(system;frequencies=[50.0,1e5],
            earth_props=homogeneous(rho=10.0,eps_r=10.0))
        result = compute(problem,formulation)
        for index in eachindex(problem.frequencies)
            @test all(isfinite,Y(result)[:,:,index])
            @test Y(result)[:,:,index] ≈ transpose(Y(result)[:,:,index]) rtol=1e-10
            @test all(real(Y(result)[i,i,index]) >= 0 for i in 1:2)
        end
        paths = result.details.fem.run.map_paths
        @test length(paths) == 9 * 2 * 2
        for (field,label) in (("e","E_t [V/m]; transverse drive 1 A/m"),
            ("ez","Ez [V/m]; axial drive 1 A"),
            ("jm","|J_t| [A/m2]; transverse drive 1 A/m"))
            path = only(filter(p -> basename(p) == "$(field)_f0001_b0001.pos",paths))
            @test occursin(label,read(path,String))
        end
        rm(result.details.fem.run.run_directory;recursive=true)
    end
end
