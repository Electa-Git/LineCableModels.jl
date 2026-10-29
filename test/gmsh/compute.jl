@testset "Gmsh computation" begin
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    facts = Dict(
        "two_bare_wires" => (1e-6, (Z = 0.00568, Y = 0.00798)),
        "two_insulated_wires" => (1e-6, (Z = 0.00407, Y = 0.000447)),
        "coaxial_circular_strands" => (1e-2, (Z = 0.0492, Y = 0.105)),
        "coaxial_rectangular_strands" => (1e-2, (Z = 0.0522, Y = 0.104)),
        "overhead_bare_wires" => (1e-6, (Z = 0.00610, Y = 0.00817)),
        "mixed_bare_wires" => (1e-6, (Z = 0.00609, Y = 0.00617)))
    @testset "$case" for (case, reference) in REFERENCES
        tolerance, expected = facts[case]
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature, frequencies = reference.frequencies[1:10:end])
        formulation = FEM(; bundle = false, kron = false, transposition = false)
        n = length(problem.system.terminals)
        targets = Dict((source, quantity) => zeros(ComplexF64, n, n, 11)
        for source in ("fem", "analytical"), quantity in ("Z", "Y"))
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "matrices.csv"))
            fields = split(line, ',')
            haskey(targets, (fields[1], get(fields, 3, ""))) || continue
            f = parse(Float64, fields[2])
            k = findfirst(==(f), problem.frequencies)
            k === nothing && continue
            row = parse(Int, fields[4])
            column = parse(Int, fields[5])
            value = complex(parse(Float64, fields[6]), parse(Float64, fields[7]))
            targets[(fields[1], fields[3])][row, column, k] = value
        end
        mktempdir() do folder
            cd(folder) do
                path = extension.directory(problem)
                @test dirname(path) == joinpath(folder, ".linecablemodels", "fem", case)
                @test occursin(r"^\d{8}-\d{6}-[a-z0-9]{4}$", basename(path))
                @test !ispath(path)
                result = compute(
                    problem, formulation; trace = true, keep_run_directory = true,
                    getdp_executable = :download)
                root = dirname(path)
                retained = only(readdir(root; join = true))
                lines = readlines(joinpath(retained, "logs", "getdp-info.log"))
                @test lines[1] == realpath(extension.GetdpBinary{:download}())
                @test lines[2] == string(extension.GetdpBinary{:download}())
                @test result.frequencies == problem.frequencies
                @test result.terminals ==
                      ["cable:$(row.cable):$(problem.system.designs[row.cable].terminals[row.terminal].name)"
                       for row in problem.system.terminals]
                @test result.line_length == problem.system.line_length
                for (quantity, computed, status) in
                    (("Z", result.Z, expected.Z), ("Y", result.Y, expected.Y))
                    target = targets[("fem", quantity)]
                    for k in eachindex(problem.frequencies)
                        @test maximum(abs, computed[:, :, k] - target[:, :, k]) <=
                              tolerance * maximum(abs(target[i, i, k]) for i in 1:n)
                    end
                    analytical = targets[("analytical", quantity)]
                    errors = [maximum(abs, computed[:, :, k] - analytical[:, :, k]) /
                              maximum(abs(analytical[i, i, k]) for i in 1:n)
                              for k in eachindex(problem.frequencies)]
                    @test abs(maximum(errors) - status) <= 0.002
                end
                @test size(result.trace.Z) == (n, n, 11)
                @test size(result.trace.P) == (n, n, 11)
                model = extension.FieldModel(problem, formulation)
                for k in eachindex(problem.frequencies)
                    primitives = read(retained, model, k)
                    @test result.trace.Z[:, :, k] == primitives.Z
                    @test result.trace.P[:, :, k] == primitives.P
                end
                for name in
                    (joinpath("input", "getdp"), "mesh", joinpath("results", "jobs"),
                    "logs", "work")
                    @test isdir(joinpath(retained, name))
                end
                sources = readdir(joinpath(@__DIR__, "..", "..", "ext",
                    "LineCableModelsGmshExt", "getdp"))
                @test readdir(joinpath(retained, "input", "getdp")) == sources
                for name in ("model_data.pro", "physical.geo", "model.geo")
                    @test isfile(joinpath(retained, "input", name))
                end
                for name in ("Z.tsv", "P.tsv")
                    @test isfile(joinpath(retained, "results", name))
                end
                for k in eachindex(problem.frequencies)
                    index = lpad(k, 4, '0')
                    @test isfile(joinpath(retained, "mesh", "frequency-$index.msh"))
                    @test isfile(joinpath(retained, "logs", "mesh-$index.log"))
                    @test isfile(joinpath(retained, "logs", "mesh-$index-metadata.txt"))
                    @test isfile(joinpath(retained, "logs", "frequency-$index.log"))
                    for b in 1:n, suffix in (".done", "-Z.tsv", "-P.tsv")

                        basis = lpad(b, 4, '0')
                        @test isfile(joinpath(retained, "results", "jobs",
                            "getdp-f$index-b$basis$suffix"))
                    end
                end
                @test isempty(readdir(joinpath(retained, "work")))
            end
        end
        mktempdir() do folder
            cd(folder) do
                single = LineParametersProblem(reference.system, reference.earth;
                    temperature = reference.temperature,
                    frequencies = problem.frequencies[1:1])
                result = compute(
                    single, formulation; trace = false, keep_run_directory = false,
                    getdp_executable = :download)
                @test size(result.trace.Z) == (0, 0, 0)
                @test size(result.trace.P) == (0, 0, 0)
                @test isempty(readdir(joinpath(folder, ".linecablemodels", "fem")))
            end
        end
    end
    @testset "Caller-owned Gmsh session" begin
        reference = last(first(REFERENCES))
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature, frequencies = reference.frequencies[1:1])
        mktempdir() do folder
            cd(folder) do
                gmsh.initialize()
                try
                    gmsh.model.add("caller")
                    gmsh.option.setNumber("General.Terminal", 0)
                    gmsh.option.setNumber("General.Verbosity", 4)
                    models = gmsh.model.list()
                    current = gmsh.model.getCurrent()
                    terminal = gmsh.option.getNumber("General.Terminal")
                    verbosity = gmsh.option.getNumber("General.Verbosity")
                    compute(problem, FEM(); getdp_executable = :download)
                    @test Bool(gmsh.isInitialized())
                    @test gmsh.model.list() == models
                    @test gmsh.model.getCurrent() == current
                    @test gmsh.option.getNumber("General.Terminal") == terminal
                    @test gmsh.option.getNumber("General.Verbosity") == verbosity
                finally
                    gmsh.finalize()
                end
            end
        end
    end
    @testset "Float32 inputs" begin
        thin = build(CableDesign, "float32",
            Terminal(:core,
                Region(:core_metal, Circle(1.5f-3),
                    material(kind = :conductor, rho = 1.7241f-8))))
        system = build(CableSystem, [thin], [Placement(0.0f0, -1.0f0)];
            connections = [(core = 1,)], id = "float32", line_length = 1.0f0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        mktempdir() do folder
            cd(folder) do
                result = compute(problem, FEM(); getdp_executable = :download)
                @test result isa LineParameters{Float64}
                @test all(isfinite, result.Z)
                @test all(isfinite, result.Y)
            end
        end
    end
    @testset "Measurement inputs" begin
        reference = last(first(REFERENCES))
        copper = Material(:conductor, 1.7241e-8, 1, 0.999994, 20, 0.00393)
        wire = build(CableDesign, "bare_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425 ± 1e-4), copper)))
        system = build(CableSystem, [wire, wire],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
        frequencies = reference.frequencies[1:1]
        problem = LineParametersProblem(system, reference.earth; frequencies)
        deterministic = LineParametersProblem(reference.system, reference.earth; frequencies)
        mktempdir() do folder
            cd(folder) do
                result = compute(problem,
                    FEM(; bundle = false, kron = false, transposition = false);
                    getdp_executable = :download)
                expected = compute(deterministic,
                    FEM(; bundle = false, kron = false, transposition = false);
                    getdp_executable = :download)
                @test result isa LineParameters{Float64}
                for quantity in (:Z, :Y)
                    actual = getfield(result, quantity)
                    target = getfield(expected, quantity)
                    scaled = maximum(abs, actual - target) /
                             maximum(abs(target[i, i, 1]) for i in axes(target, 1))
                    @test scaled <= 1e-6
                end
            end
        end
    end
    @testset "Conductor reductions" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, 1.0e14, 2.3)
        design = build(CableDesign, "coaxial",
            Terminal(:core, Region(:core_metal, Circle(10e-3), copper)),
            Region(:core_insulation, Offset(10e-3), polymer),
            Terminal(:sheath, Region(:sheath_copper, Annulus(20e-3, 21e-3), copper)),
            Region(:sheath_jacket, Offset(3e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1, sheath = 0), (core = 1, sheath = 2)],
            id = "two cables/grounded", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        mktempdir() do folder
            cd(folder) do
                path = extension.directory(problem)
                @test basename(dirname(path)) == "two_cables_grounded"
                result = compute(problem, FEM(); trace = true, getdp_executable = :download)
                coaxial = compute(problem, CoaxialFormulation())
                @test result.terminals == coaxial.terminals
                @test result.terminals == ["cable:1:core+cable:2:core", "cable:2:sheath"]
                @test size(result.Z) == (2, 2, 1)
                @test size(result.Y) == (2, 2, 1)
                C = [1.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0]
                expected_Z = inv(transpose(C) * inv(result.trace.Z[:, :, 1]) * C)
                expected_Y = transpose(C) * inv(result.trace.P[:, :, 1]) * C
                @test norm(result.Z[:, :, 1] - expected_Z) <=
                      1e-12 * norm(result.Z[:, :, 1])
                @test norm(result.Y[:, :, 1] - expected_Y) <=
                      1e-12 * norm(result.Y[:, :, 1])
                @test isempty(readdir(joinpath(folder, ".linecablemodels", "fem")))
            end
        end
    end
    @testset "FEM transposition" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, 1.0e14, 2.3)
        design = build(CableDesign, "coaxial",
            Terminal(:core, Region(:core_metal, Circle(10e-3), copper)),
            Region(:core_insulation, Offset(10e-3), polymer),
            Terminal(:sheath, Region(:sheath_copper, Annulus(20e-3, 21e-3), copper)),
            Region(:sheath_jacket, Offset(3e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1, sheath = 0), (core = 1, sheath = 2)],
            id = "two cables/grounded", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        mktempdir() do folder
            cd(folder) do
                transposed = compute(problem,
                    FEM(; bundle = false, kron = false, transposition = true);
                    getdp_executable = :download)
                unreduced = compute(problem,
                    FEM(; bundle = false, kron = false, transposition = false);
                    getdp_executable = :download)
                for (computed, original) in ((transposed.Z, unreduced.Z),
                    (transposed.Y, unreduced.Y))
                    @test all(norm(computed[:, :, k] -
                                   sum(original[circshift(1:4, j), circshift(1:4, j), k]
                              for j in 0:3) / 4) <= 1e-14 * norm(original[:, :, k])
                    for k in eachindex(problem.frequencies))
                end
            end
        end
    end
end
