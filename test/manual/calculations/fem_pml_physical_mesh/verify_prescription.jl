using LineCableModels, Gmsh, JSON3, TOML, Test
const FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
const ROOT=joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh")
const controls=(domain_skin_depths=24.,pml_resolution=(interpolation_cells=72,coefficient_change=.12),
    mesh_size_factor=3.,exterior_mesh_size_factor=8.,volume_quadrature=12,
    conductor_geometry_tolerance=1e-3,conductor_skin_depth_elements=6.,
    conductor_mesh_growth=sqrt(1.25),conductor_skin_depths=5.,conductor_thickness_elements=4)
const form=Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
@testset "Production matches qualified prescriptions" begin
    options=computation_options(LineCableModelsFEM,ComputationOptions(;controls...))
    @test options.data.pml_resolution==(interpolation_cells=72,coefficient_change=.12)
    for label in readdir(joinpath(ROOT,"qualification"))
        marker=TOML.parsefile(joinpath(ROOT,"qualification",label,"complete.toml"))
        problem=LineCableModels.ImportExport.deserialize_value(JSON3.read(read(joinpath(marker["run_directory"],"input/problem.json"),String),Dict{String,Any}))
        model=FEM._resolved_fem_model(problem,form,options)
        for plan in model.mesh_plans
            frozen=TOML.parsefile(joinpath(ROOT,"resolved",label,string("strict_guard-",plan.frequency,".toml")))
            for (direction,strips) in zip(("side","top","bottom"),plan.pml_strips)
                expected=frozen[direction]
                @test length(strips)==length(expected)
                for (strip,row) in zip(strips,expected)
                    @test strip.count==row["count"]
                    @test strip.start≈row["start"] rtol=2e-12 atol=1e-14
                    @test strip.stop≈row["stop"] rtol=2e-12 atol=1e-14
                    @test strip.ratio≈row["ratio"] rtol=2e-12 atol=1e-14
                end
            end
        end
        println("LAYOUT PRESERVATION ",label," ",length(model.mesh_plans)," frequencies");flush(stdout)
    end
end
