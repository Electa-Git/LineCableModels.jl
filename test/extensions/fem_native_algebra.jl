@testitem "Gmsh FEM / native matrix reductions and runtime complex inversion" tags=[:extension,:fem_numerical] begin
    using Gmsh, DelimitedFiles, LinearAlgebra
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    Engine=LineCableModels.Engine
    getdp=FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    fixture=joinpath(pkgdir(LineCableModels),"test","fixtures","data","fem","native_getdp")
    algebra=joinpath(pkgdir(LineCableModels),"ext","LineCableModelsGmshExt","getdp","line-parameters.pro")
    mktempdir() do dir
        gmsh=Gmsh.gmsh
        gmsh.initialize(String[],false,false)
        try
            gmsh.option.set_number("General.Terminal",0)
            gmsh.open(joinpath(fixture,"model.geo"))
            gmsh.model.mesh.generate(2)
            gmsh.write(joinpath(dir,"model.msh"))
        finally
            gmsh.finalize()
        end
        model=joinpath(dir,"matrices.pro")
        cp(joinpath(fixture,"matrices.pro"),model)
        mesh=joinpath(dir,"model.msh")
        for n in (1,2,4)
            maps=n==4 ? [[1,2,3,4],[2,1,2,0]] : [collect(1:n)]
            for connections in maps, bundle in (false,true), kron in (false,true), transposed in (false,true)
                case=connections==[2,1,2,0] ? 1 : 0
                cmd=`$getdp $model -msh $mesh -solve Matrices -name $(joinpath(dir,"algebra")) -setnumber LinearSolver 0 -setstring AlgebraPath $algebra -setnumber NumTerminals $n -setnumber ConnectionCase $case -setnumber ReduceBundle $(Int(bundle)) -setnumber KronReduction $(Int(kron)) -setnumber IdealTransposition $(Int(transposed)) -v 2`
                success(pipeline(Cmd(cmd;dir),stdout=joinpath(dir,"getdp.log"),stderr=joinpath(dir,"getdp-errors.log"))) ||
                    error(read(joinpath(dir,"getdp.log"),String)*read(joinpath(dir,"getdp-errors.log"),String))
                rows=readdlm(joinpath(dir,"matrices.tsv"))
                P=[complex((i==j)*5+.1*i+.2*j,.15*i-.07*j) for i in 1:n,j in 1:n]
                Zp=[complex((i==j)*7+.2*i-.1*j,.12*i+.06*j) for i in 1:n,j in 1:n]
                options=Formulation(:LineCableModelsFEM;options=(reduce_bundle=bundle,kron_reduction=kron,ideal_transposition=transposed)).options
                for pass in 1:2
                    expected=Engine.reduce_primitive_matrices(reshape(pass*Zp,n,n,1),reshape(pass*P,n,n,1),connections,options)
                    sizeout=size(expected.P,1)
                    observed=[zeros(ComplexF64,sizeout,sizeout) for _ in 1:3]
                    selected=rows[rows[:,1].==pass,:]
                    @test size(selected,1)==sizeout^2
                    for row in eachrow(selected), matrix in 1:3
                        observed[matrix][Int(row[2]),Int(row[3])]=complex(row[2matrix+2],row[2matrix+3])
                    end
                    @test observed[1] ≈ expected.Z[:,:,1] atol=1e-13 rtol=0
                    @test observed[2] ≈ expected.P[:,:,1] atol=1e-13 rtol=0
                    @test observed[3] ≈ inv(expected.P[:,:,1]) atol=1e-13 rtol=0
                    @test maximum(abs,observed[2]*observed[3]-I) <= 1e-13
                end
            end
        end
        singular=`$getdp $model -msh $mesh -solve Matrices -name $(joinpath(dir,"algebra")) -setnumber LinearSolver 0 -setstring AlgebraPath $algebra -setnumber Singular 1 -v 2`
        @test !success(pipeline(Cmd(singular;dir),stdout=joinpath(dir,"singular.log"),stderr=joinpath(dir,"singular-errors.log")))
    end
end
