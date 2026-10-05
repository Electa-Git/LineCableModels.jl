@testitem "Gmsh FEM / prescribed Gamma options and frequency transport" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    @test Formulation(:LineCableModelsFEM).options.data.Γ == 0
    for Γ in (0, .2, .2+.3im, [0, .2+.3im])
        f = Formulation(:LineCableModelsFEM; options=(physics=:helmholtz, Γ))
        @test f.options.data.Γ == Γ
        Γ isa AbstractVector && @test f.options.data.Γ !== Γ
    end
    for Γ in (true, NaN, Inf+im, [], [0, Inf], "zero", [true])
        @test_throws ArgumentError Formulation(:LineCableModelsFEM; options=(physics=:helmholtz, Γ))
    end
    @test Formulation(:LineCableModelsFEM; options=(Γ=.1,)).options.data.Γ == .1
    wire = build(CableDesign,"gamma-input",terminal(:core,
        core(Material(kind=:conductor,rho=2e-8);r=.005)))
    system = build(LineCableSystem,[wire],[(0.,1.)]; connections=[Dict(:core=>1)])
    problem = LineParametersProblem(system;frequencies=[100.,10000.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    f = Formulation(:LineCableModelsFEM;options=(physics=:helmholtz,Γ=[0.,.01+.02im]))
    model = FEM._resolved_fem_model(problem,f)
    @test model.prescribed_gamma == f.options.data.Γ
    @test_throws DimensionMismatch FEM._resolved_fem_model(problem,
        Formulation(:LineCableModelsFEM;options=(physics=:helmholtz,Γ=[0.])))
    scalar = FEM._resolved_fem_model(problem,
        Formulation(:LineCableModelsFEM;options=(physics=:helmholtz,Γ=.01+.02im)))
    @test all(==(.01+.02im), scalar.prescribed_gamma)
    mktempdir() do dir
        path = joinpath(dir,"model_data.pro")
        FEM._write_model_data(path,model,computation_options(LineCableModelsFEM,ComputationOptions()).data)
        @test occursin("GammaReValues() = {0, 0.01};",read(path,String))
        @test occursin("GammaImValues() = {0, 0.02};",read(path,String))
    end
end


@testitem "Gmsh FEM / prescribed Gamma manufactured Maxwell fields" tags=[:extension,:fem_numerical] begin
    using Gmsh, DelimitedFiles
    const LineCableModelsFEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt).LineCableModelsFEM
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    fixture = joinpath(pkgdir(LineCableModels),"test/fixtures/data/fem/prescribed_gamma")
    reference = read(joinpath(fixture,"model.pro"),String)
    production = read(FEM._getdp_assets().helmholtz,String)
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

