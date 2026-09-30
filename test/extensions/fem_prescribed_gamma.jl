@testitem "Gmsh FEM / prescribed Gamma options and frequency transport" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    @test Formulation(:LineCableModelsFEM).options.data.Γ == 0
    for Γ in (0, .2, .2+.3im, [0, .2+.3im])
        f = Formulation(:LineCableModelsFEM; options=(physics=:quasi_fw, Γ))
        @test f.options.data.Γ == Γ
        Γ isa AbstractVector && @test f.options.data.Γ !== Γ
    end
    for Γ in (true, NaN, Inf+im, [], [0, Inf], "zero", [true])
        @test_throws ArgumentError Formulation(:LineCableModelsFEM; options=(physics=:quasi_fw, Γ))
    end
    @test Formulation(:LineCableModelsFEM; options=(Γ=.1,)).options.data.Γ == .1
    wire = build(CableDesign,"gamma-input",terminal(:core,
        core(Material(kind=:conductor,rho=2e-8);r=.005)))
    system = build(LineCableSystem,[wire],[(0.,1.)]; connections=[Dict(:core=>1)])
    problem = LineParametersProblem(system;frequencies=[100.,10000.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    f = Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,Γ=[0.,.01+.02im]))
    model = FEM._resolved_fem_model(problem,f)
    @test getproperty.(model.mesh_plans,:Γ) == f.options.data.Γ
    @test_throws DimensionMismatch FEM._resolved_fem_model(problem,
        Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,Γ=[0.])))
    scalar = FEM._resolved_fem_model(problem,
        Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,Γ=.01+.02im)))
    @test all(p -> p.Γ == .01+.02im, scalar.mesh_plans)
    # The same complex ray must damp air and soil, including transverse
    # spectral components when Im(sqrt(gamma_medium^2-Gamma^2)) is negative.
    shifted = FEM._resolved_fem_model(problem,
        Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,Γ=.02+.01im)))
    for plan in shifted.mesh_plans, (sigma,epsr) in ((0.,1.),(.01,10.)), λ in (0.,.001,.1,1.)
        q = sqrt(complex(λ^2 + 2π*plan.frequency*im*4π*1e-7*
            (sigma+2π*plan.frequency*im*epsr*8.8541878128e-12)-plan.Γ^2))
        @test real(q*(1-im*plan.pml_slope)) > 0
    end
    physical = FEM._resolved_fem_model(problem,
        Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,Γ=.02+.01im)),
        computation_options(LineCableModelsFEM,ComputationOptions(
            pml_resolution=(interpolation_cells=72,coefficient_change=.12))))
    @test all(plan -> all(strips -> !isempty(strips),plan.pml_strips),physical.mesh_plans)
    mktempdir() do dir
        path = joinpath(dir,"model_data.pro")
        FEM._write_model_data(path,model)
        @test occursin("GammaReValues() = {0, 0.01};",read(path,String))
        @test occursin("GammaImValues() = {0, 0.02};",read(path,String))
    end
end

@testitem "Gmsh FEM / finite Gamma solve, resume, and detached export" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh, LinearAlgebra, DelimitedFiles
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"gamma-current",terminal(:core,
        core(Material(kind=:conductor,rho=2e-8);r=.01)))
    system = build(LineCableSystem,[wire,wire],[(0.,1.),(.5,-1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    frequency = 1e4
    Γ = .99sqrt(im*2π*frequency*4π*1e-7*(.01+im*2π*frequency*10*8.8541878128e-12))
    problem = LineParametersProblem(system;frequencies=[frequency,frequency],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    reductions = (reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    f = Formulation(:LineCableModelsFEM;options=(;reductions...,physics=:quasi_fw,Γ=[0.,Γ]))
    mesh_controls = (pml_layers=16,pml_grading=3.,pml_reflection=1e-6,mesh_size_factor=2.)
    controls = (;mesh_controls...,gmsh_verbosity=0,getdp_verbosity=2,
        keep_run_directory=true,trace=true,frequency_workers=1,plot_field_maps=true)
    result = compute(problem,f;options=controls)
    @test all(isfinite,Z(result)) && all(isfinite,Y(result))
    @test !isapprox(Z(result)[:,:,1],Z(result)[:,:,2];rtol=1e-4)
    root = result.details.data.fem.run.run_directory
    @test result.details.data.fem.inputs.options.Γ == [0.,Γ]
    # Integrate the actual finite-Gamma axial-current maps over each metal.
    # This checks that the shifted drive is reconstructed consistently in
    # physical outputs, independently of managed/detached matrix agreement.
    model = FEM._resolved_fem_model(problem,f)
    metal_tags = Set(m.physical_tag for m in model.material_plans if m.kind===:conductor)
    maps = result.details.data.fem.run.map_paths
    for basis in 1:2
        suffix = "f0002_b000$(basis).pos"
        current_path = only(filter(p->endswith(p,"/jz_"*suffix),maps))
        material_path = only(filter(p->endswith(p,"/material_region_"*suffix),maps))
        integrated = zeros(ComplexF64,2)
        for (current_line,material_line) in zip(eachline(current_path),eachline(material_path))
            current = match(r"^ST\(([^)]*)\)\{([^}]*)\};",current_line)
            current===nothing && continue
            material = match(r"^ST\(([^)]*)\)\{([^}]*)\};",material_line)
            @test material!==nothing && current[1]==material[1]
            tag = parse(Float64,first(split(material[2],',')))
            tag in metal_tags || continue
            xyz = reshape(parse.(Float64,split(current[1],',')),3,3)
            area = abs((xyz[1,2]-xyz[1,1])*(xyz[2,3]-xyz[2,1])-
                (xyz[1,3]-xyz[1,1])*(xyz[2,2]-xyz[2,1]))/2
            receiver = sum(xyz[2,:])>0 ? 1 : 2
            values = parse.(Float64,split(current[2],','))
            integrated[receiver] += area*complex(sum(values[1:3]),sum(values[4:6]))/3
        end
        @test integrated ≈ ComplexF64[basis==1,basis==2] atol=5e-7 rtol=0
    end
    repeated = compute(problem,f;options=(;controls...,resume_run_directory=root))
    @test Z(repeated) == Z(result) && Y(repeated) == Y(result)
    @test repeated.details.data.fem.run.reused
    other = Formulation(:LineCableModelsFEM;options=(;reductions...,physics=:quasi_fw,Γ=[0.,Γ/2]))
    @test_throws ArgumentError compute(problem,other;
        options=(;controls...,resume_run_directory=root))
    # Use the same supplied mesh in the detached native runtime. This isolates
    # Gamma transport and matrix extraction from independent mesher ordering.
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    mktempdir() do dir
        entry = export_data(:onelab,problem,f;file_name=joinpath(dir,"bundle","study.pro"),
            mesh_options=mesh_controls)
        bundle = dirname(entry)
        mesh = joinpath(root,"mesh","model.msh")
        cmd = `$getdp $entry -msh $mesh -solve LineCableModelsFEM -setnumber FrequencyIndex 2 -setnumber PlotFieldMaps 0 -v 2 -nt 1 -ksp_diagonal_scale -ksp_diagonal_scale_fix`
        ok = success(pipeline(cmd,stdout=joinpath(dir,"getdp.log"),stderr=joinpath(dir,"errors.log")))
        ok || error(read(joinpath(dir,"getdp.log"),String)*read(joinpath(dir,"errors.log"),String))
        matrices = joinpath(bundle,"results","f0002-quasi-fw-b0000","matrices")
        @test isfile(joinpath(dirname(matrices),"completed.txt"))
        for (quantity,expected) in (("Z",Z(result)[:,:,2]),("Y",Y(result)[:,:,2]))
            rows = [split(line,'\t') for line in readlines(joinpath(matrices,"$quantity.tsv"))[3:end]]
            actual = zeros(ComplexF64,2,2)
            for row in rows
                actual[parse(Int,row[1]),parse(Int,row[2])] =
                    complex(parse(Float64,row[5]),parse(Float64,row[6]))
            end
            @test actual ≈ expected rtol=1e-7
        end
    end
end

@testitem "Gmsh FEM / prescribed Gamma manufactured Maxwell fields" tags=[:extension,:fem_numerical] begin
    using Gmsh, DelimitedFiles
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    fixture = joinpath(pkgdir(LineCableModels),"test/fixtures/data/fem/prescribed_gamma")
    reference = read(joinpath(fixture,"model.pro"),String)
    production = read(FEM._getdp_assets().quasi_full,String)
    # Exercise the actual exterior equation block with analytic body sources.
    # Terminal constraints are tested separately on native conductor geometry.
    equations = split(split(production,"    Equation {";limit=2)[2],
        "\n    }\n  }\n}\n\nMacro";limit=2)[1]
    equations = replace(equations,r"GlobalTerm \{[^\n]+"=>"",
        "DomainMedia_Ele"=>"Domain","DomainFields"=>"Domain","gamma2[]"=>"g2[]")
    prefix = """
    PerfectConductors=1; FiniteGamma=GammaRe != 0 || GammaIm != 0;
    Function { nuPml[]=nuTensor[]; sePml[]=seTensor[]; seZ[]=se[]*detS[]; }
    """
    sources = """
      Galerkin { [-fa[]*Vector[0,0,1], {a}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-fb[], {bt}]; In Domain; Jacobian Vol; Integration I1; }
      Galerkin { [-fa[], {v}]; In Domain; Jacobian Vol; Integration I1; }
    """
    head,tail = split(reference,"Formulation {";limit=2)
    pre,tail = split(tail,"    Equation {";limit=2)
    _,post = split(tail,"\n    }\n  }\n}\nResolution";limit=2)
    actual = head*prefix*"Formulation {"*pre*"    Equation {"*equations*sources*
        "\n    }\n  }\n}\nResolution"*post
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    gammas = (0.,2+im,-2-im,1e-6*(2+im),2im,2.)
    mktempdir() do root
        errors = Dict()
        for n in (12,24)
            mesh = joinpath(root,"mesh$n.msh")
            geo = joinpath(root,"mesh$n.geo")
            write(geo,"N=$n;\n"*read(joinpath(fixture,"model.geo"),String))
            session = FEM._start_gmsh(0)
            try
                Gmsh.gmsh.open(geo)
                Gmsh.gmsh.model.mesh.generate(2)
                @test length(first(Gmsh.gmsh.model.mesh.get_nodes())) == (n+1)^2
                Gmsh.gmsh.write(mesh)
            finally
                FEM._finish_gmsh(session)
            end
            for (stretched,(sx,sy)) in enumerate(((0.,0.),(.3,-.5))),
                    (case,Γ) in enumerate(gammas), (kind,code) in enumerate((reference,actual))
                dir = mkpath(joinpath(root,"n$n-s$stretched-g$case-k$kind"))
                write(joinpath(dir,"model.pro"),code)
                cmd = `$getdp model.pro -msh $mesh -solve Check -setnumber GammaRe $(real(Γ)) -setnumber GammaIm $(imag(Γ)) -setnumber StretchX $sx -setnumber StretchY $sy -v 2 -nt 1 -ksp_diagonal_scale -ksp_diagonal_scale_fix`
                ok = success(pipeline(Cmd(cmd;dir),stdout=joinpath(dir,"getdp.log"),
                    stderr=joinpath(dir,"errors.log")))
                ok || error(read(joinpath(dir,"getdp.log"),String)*read(joinpath(dir,"errors.log"),String))
                e = [readdlm(joinpath(dir,"$field.tsv"))[2] for field in ("et","ez","bt","bz")]
                @test all(x -> isfinite(x) && x >= 0,e)
                errors[n,stretched,case,kind] = e
                kind==2 && @test e ≈ errors[n,stretched,case,1] rtol=1e-9
                if n==24
                    ratio=e./errors[12,stretched,case,kind]
                    @test all(<(.28),ratio[[1,3,4]]) # First-order E_t, B_t, B_z.
                    @test ratio[2] < .08             # Second-order E_z.
                end
            end
        end
        for stretch in 1:2, kind in 1:2
            @test errors[24,stretch,2,kind] == errors[24,stretch,3,kind]
            @test errors[24,stretch,1,kind] ≈ errors[24,stretch,4,kind] rtol=1e-8
        end
    end
end

@testitem "Gmsh FEM / prescribed Gamma cylindrical Maxwell reference" tags=[:extension,:fem_numerical] begin
    using Gmsh, SpecialFunctions, DelimitedFiles
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    fixture = joinpath(pkgdir(LineCableModels), "test/fixtures/data/fem/prescribed_gamma")
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM, ComputationOptions())).path
    model = FEM._getdp_assets().model
    data = joinpath(fixture, "coax-data.pro")
    q0 = sqrt(1.7im * 1.25 * (.7 + .68im))
    previous = Dict()
    mktempdir() do root
        for (level, h) in enumerate((.16, .08)), (case, Γ) in enumerate((0., .99q0, .6+.3im))
            dir = mkpath(joinpath(root, "h$level-g$case"))
            mesh = joinpath(dir, "coax.msh")
            geo = joinpath(dir, "coax.geo")
            write(geo, "H=$h;\n" * read(joinpath(fixture, "coax.geo"), String))
            session = FEM._start_gmsh(0)
            try
                Gmsh.gmsh.open(geo)
                Gmsh.gmsh.model.mesh.generate(2)
                Gmsh.gmsh.write(mesh)
            finally
                FEM._finish_gmsh(session)
            end
            cmd = `$getdp $model -msh $mesh -name $(joinpath(dir,"solver")) -solve LineCableModelsFEMScan -setnumber Physics 1 -setnumber GammaRe $(real(Γ)) -setnumber GammaIm $(imag(Γ)) -setstring ModelDataPath $data -setstring RunDirectory $dir -v 2 -nt 1 -ksp_diagonal_scale -ksp_diagonal_scale_fix`
            ok = success(pipeline(cmd, stdout=joinpath(dir,"getdp.log"), stderr=joinpath(dir,"errors.log")))
            ok || error(read(joinpath(dir,"getdp.log"),String) * read(joinpath(dir,"errors.log"),String))
            # Homogeneous annulus, current-carrying inner PEC and zero outer
            # potentials: the I0/K0 combination enforces the exact boundary.
            q = sqrt(1.7im*1.25*(.7+.68im) - Γ^2)
            r, R = .2, 2.
            ratio = besselk(0,q*R)/besseli(0,q*R)
            a = 1.25*(besselk(0,q*r)-ratio*besseli(0,q*r))/
                (2π*r*q*(besselk(1,q*r)+ratio*besseli(1,q*r)))
            exact = (1.7im*a, a/(1.25*(.7+.68im)))
            measured = map(("Z", "P")) do quantity
                row = readdlm(joinpath(dir,"raw/jobs/getdp-f0001-b0001-$quantity.tsv"))
                complex(row[5],row[6])
            end
            errors = abs.(measured ./ exact .- 1)
            @test all(<(.0025), errors)
            if level == 2
                @test all(<(.3), errors ./ previous[case])
                @test all(<(.0006), errors)
            end
            previous[case] = errors
        end
    end
end

@testitem "Gmsh FEM / prescribed Gamma outgoing Cartesian PML" tags=[:extension,:fem_numerical] begin
    using Gmsh, SpecialFunctions, DelimitedFiles
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    fixture = joinpath(pkgdir(LineCableModels), "test/fixtures/data/fem/prescribed_gamma")
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM, ComputationOptions())).path
    model = FEM._getdp_assets().model
    q0 = sqrt(1.7im * 1.25 * (.7 + .68im))
    previous = Dict()
    mktempdir() do root
        for (level, h) in enumerate((.16, .08)), (case, Γ) in enumerate((0., .99q0, 2+2im))
            dir = mkpath(joinpath(root, "h$level-g$case"))
            mesh = joinpath(dir, "outgoing.msh")
            geo = joinpath(dir, "outgoing.geo")
            write(geo, "H=$h;\n" * read(joinpath(fixture, "outgoing.geo"), String))
            session = FEM._start_gmsh(0)
            try
                Gmsh.gmsh.open(geo)
                Gmsh.gmsh.model.mesh.generate(2)
                Gmsh.gmsh.write(mesh)
            finally
                FEM._finish_gmsh(session)
            end
            q = sqrt(q0^2 - Γ^2)
            slope = imag(q) < 0 ? min(1., real(q)/(2abs(imag(q)))) : 1.
            strength = 4*(-log(1e-8)/2)/(2*(real(q)+slope*imag(q)))
            data = joinpath(dir, "data.pro")
            write(data, read(joinpath(fixture, "coax-data.pro"), String) * """
            DomainHalfwidth=1.; PmlSideThickness=2.; PmlTopThickness=2.; PmlBottomThickness=2.;
            PmlSlope=$slope; PmlSideStrength=$strength; PmlTopStrength=$strength; PmlBottomStrength=$strength;
            """)
            cmd = `$getdp $model -msh $mesh -name $(joinpath(dir,"solver")) -solve LineCableModelsFEMScan -setnumber Physics 1 -setnumber GammaRe $(real(Γ)) -setnumber GammaIm $(imag(Γ)) -setstring ModelDataPath $data -setstring RunDirectory $dir -v 2 -nt 1 -ksp_diagonal_scale -ksp_diagonal_scale_fix`
            ok = success(pipeline(cmd, stdout=joinpath(dir,"getdp.log"), stderr=joinpath(dir,"errors.log")))
            ok || error(read(joinpath(dir,"getdp.log"),String) * read(joinpath(dir,"errors.log"),String))
            # Infinite homogeneous exterior: only the decaying K0 solution is
            # admissible. This checks the variable PML tensors and the voltage
            # path through the PML against a reference without truncation.
            a = 1.25*besselk(0,.2q)/(2π*.2q*besselk(1,.2q))
            exact = (1.7im*a, a/(1.25*(.7+.68im)))
            measured = map(("Z", "P")) do quantity
                row = readdlm(joinpath(dir,"raw/jobs/getdp-f0001-b0001-$quantity.tsv"))
                complex(row[5],row[6])
            end
            errors = abs.(measured ./ exact .- 1)
            @test all(<(.03), errors)
            if level == 2
                @test all(<(.4), errors ./ previous[case])
                @test all(<(.01), errors)
            end
            previous[case] = errors
        end
    end
end
