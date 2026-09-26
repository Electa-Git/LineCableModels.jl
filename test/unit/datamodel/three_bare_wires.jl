@testitem "DataModel / reusable three bare wire placements" tags=[:unit] setup=[TestFixtures] begin
    layouts = TestFixtures.three_bare_wires_layouts
    @test length(layouts) == 5
    for (name, heights) in pairs(layouts)
        problem = TestFixtures.three_bare_wires_problem(; heights, name=string(name))
        system = problem.system
        @test getproperty.(system.positions, :x) == [0.0, 1.0, 2.0]
        @test getproperty.(system.positions, :y) == collect(heights)
        @test system.positions == system.input_positions
        @test system.connection_order == [1, 2, 3]
        @test system.terminal_order == [(cable=i, terminal=:core) for i in 1:3]
        @test length(system.geometry) == 3
        @test all(design -> design.terminal_order == [:core], system.designs)
        @test system.line_length == 1.0
        @test problem.temperature == 20.0
        @test problem.frequencies == [0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]
        @test !hasproperty(problem, :Γ)
        @test problem.earth_props.layers[2].rho == 0.1
        @test problem.earth_props.layers[2].eps_r == 1.0
        @test problem.earth_props.layers[2].mu_r == 1.0
        mktempdir() do directory
            path = joinpath(directory, "problem.json")
            export_data(:json, problem; file_name=path)
            restored = import_data(:json, LineParametersProblem; file_name=path)
            @test restored.system.positions == system.positions
            @test restored.system.terminal_order == system.terminal_order
            @test restored.frequencies == problem.frequencies
        end
    end
    changed = TestFixtures.three_bare_wires_problem(; radius=0.01, rho=100.0,
        horizontal=(0.0, 2.0, 4.0), heights=(2.0, -3.0, -4.0), frequencies=[50.0])
    @test getproperty.(changed.system.positions, :x) == [0.0, 2.0, 4.0]
    @test getproperty.(changed.system.positions, :y) == [2.0, -3.0, -4.0]
    @test changed.earth_props.layers[2].rho == 100.0
    @test changed.frequencies == [50.0]
    @test_throws ArgumentError TestFixtures.three_bare_wires_problem(; heights=(1.0,))
end
