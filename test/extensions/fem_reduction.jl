@testitem "Gmsh FEM / P/(jω) through the FEM call site gives the coaxial admittance" tags=[:extension] begin
    using Gmsh
    include(joinpath(pkgdir(LineCableModels), "test", "support", "scenarios.jl"))
    using .CurrentScenarios
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    frequencies = [50.0, 1e3, 1e5]
    problem = three_bare_wires_problem(; heights = (10.0, 10.0, 10.0),
        horizontal = (0.0, 4.0, 8.0), rho = 100.0, frequencies)
    execution = computation_options(LineCableModelsFEM, ComputationOptions((;)))
    for ideal_transposition in (false, true)
        options = (; ideal_transposition)
        coaxial = compute(problem, Formulation(; options); options = (trace = true,))
        trace = details(coaxial).data.trace
        # FEM extracts the inverse-admittance coefficient P/(jω) in Ω·m from GetDP.
        extracted = similar(trace.P)
        for (k, f) in pairs(frequencies)
            extracted[:, :, k] = trace.P[:, :, k] / (im * 2π * f)
        end
        formulation = Formulation(:LineCableModelsFEM; options)
        model = FEM._resolved_fem_model(problem, formulation)
        run = FEM._create_run(mktempdir())
        # The input record only enters the result details.
        fem = FEM._line_parameters(run, model, formulation, execution,
            FEM.FEMScan(copy(trace.Z), extracted, String[]), (;))
        @test Z(fem) == Z(coaxial)
        @test Y(fem) ≈ Y(coaxial) rtol = 1e-12
        @test all(≤(1e-10), details(fem).data.fem.inversion_residuals)
    end
end
