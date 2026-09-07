@testitem "Gmsh FEM / data-only GetDP input transport" tags=[:extension] begin
    using LineCableModels, Gmsh

    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    copper = Material(kind=:conductor, rho=1.72e-8, mu_r=0.999994)
    dielectric = Material(kind=:insulator, rho=1e8, eps_r=2.3, tan_delta=0.025)
    design = build(CableDesign, "transport", Stack(
        Group(:core, Region(:metal, Disk(0.005), copper)),
        Region(:insulation, Shell(0.005), dielectric)))
    system = build(LineCableSystem, [design, design], [(0.0, -0.1), (0.1, -0.1)];
        connections=[Dict(:core=>1), Dict(:core=>2)])
    for method in (:default, :Ametani2004)
        formulation = Formulation(:LineCableModelsFEM; insulation_admittance=method,
            options=(ideal_transposition=false,))
        line_counts = Int[]
        for frequencies in ([50.0, 1000.0], collect(10.0 .^ range(-1, 6; length=101)))
            problem = LineParametersProblem(system; frequencies,
                earth_props=LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0))
            model = extension._resolved_fem_model(problem, formulation)
            mktempdir() do directory
                path = joinpath(directory, "model_data.pro")
                @test extension._write_model_data(path, model) == path
                @test readdir(directory) == ["model_data.pro"]
                text = read(path, String)
                push!(line_counts, length(readlines(path)))
                # Input serialization must not grow a second solver/output program.
                @test !occursin(r"(?m)^\s*(Group|Function|Macro|PostOperation|For|If|Include)\b", text)
                @test !occursin("MaterialSigma()", text)
                @test !occursin("MaterialEpsilon()", text)
                @test !occursin("MaterialTanDelta", text)
                for (index, material) in enumerate(model.material_plans)
                    for (name, expected) in (
                        ("MaterialSigma_$index", real.(material.admittivity)),
                        ("MaterialEpsilon_$index", imag.(material.admittivity) ./ (2π .* frequencies)))
                        entry = match(Regex(name * raw"\(\) = \{([^}]*)\};"), text)
                        @test entry !== nothing
                        @test parse.(Float64, split(entry[1], ',')) == expected
                    end
                end
                @test extension._write_model_data(path, model) == path
                @test read(path, String) == text
            end
        end
        # More samples add numerical array entries, never per-job declarations.
        @test first(line_counts) == last(line_counts)
    end
end
