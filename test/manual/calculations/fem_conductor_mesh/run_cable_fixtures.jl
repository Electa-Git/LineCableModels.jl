# Serial public-API comparison of the documented screen, tube and sector fixtures.
# julia --project=. test/manual/calculations/fem_conductor_mesh/run_cable_fixtures.jl OUTPUT
# Include instead to choose frequencies, physics and prescribed mesh controls.
using LineCableModels, Gmsh, Printf, TOML
include(joinpath(@__DIR__, "fixtures.jl"))

function fixture_case(directory, problem, formulation, mesh_options)
    mkpath(directory)
    options = (; mesh_options..., mesh_policy=:reuse, resume_run_directory=:latest,
        keep_run_directory=true, trace=true, timing=true,
        frequency_workers=1, solver_threads=1, plot_field_maps=false,
        gmsh_verbosity=4, getdp_verbosity=4, verbosity=(default=1,))
    write(joinpath(directory, "options.txt"), repr(options) * "\n")
    println("START ", directory, " frequencies=", problem.frequencies); flush(stdout)
    measured = @timed compute(problem, formulation; options)
    result = measured.value
    record = details(result).data.fem.run
    frequencies = problem.frequencies
    open(joinpath(directory, "matrices.csv"), "w") do io
        println(io, "quantity,frequency_hz,receiver,source,real,imaginary")
        for (name, values) in (("Z", Z(result)), ("Y", Y(result)))
            for k in eachindex(frequencies), j in axes(values, 2), i in axes(values, 1)
                value = values[i, j, k]
                @printf(io, "%s,%.17g,%d,%d,%.17g,%.17g\n",
                    name, frequencies[k], i, j, real(value), imag(value))
            end
        end
    end
    for quantity in ("Z", "P")
        cp(joinpath(record.run_directory, "raw", "$quantity.tsv"),
            joinpath(directory, "$quantity-primitive.tsv"); force=true)
    end
    # A resumed read is not a fresh solver timing. Keep the original cost record
    # when this invocation reuses that same completed run.
    cost_path = joinpath(directory, "cost.toml")
    previous = isfile(cost_path) ? TOML.parsefile(cost_path) : Dict()
    if !record.reused || get(previous, "run_directory", nothing) != record.run_directory
        open(cost_path, "w") do io
            TOML.print(io, Dict("wall_seconds" => measured.time,
                "compile_seconds" => measured.compile_time,
                "recompile_seconds" => measured.recompile_time,
                "julia_gc_seconds" => measured.gctime,
                "run_directory" => record.run_directory, "reused" => record.reused))
        end
    end
    println("DONE ", directory, " in ", measured.time, " s; ", record.run_directory)
    flush(stdout)
    return nothing
end

function run_cable_fixtures(root; frequencies=[.1, 50., 1e4, 1e6],
        refined_frequencies=[1e6], physics=:quasi_fw,
        mesh_options=(domain_skin_depths=24., pml_layers=192,
            mesh_size_factor=3., exterior_mesh_size_factor=8.),
        refinements=(conductor_geometry_tolerance=2.5e-4,
            conductor_skin_depth_elements=12., conductor_mesh_growth=1.25^0.25))
    formulation = Formulation(:LineCableModelsFEM; options=(;
        physics, reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
    for (name, constructor, height) in (
            ("screen", ConductorMeshFixtures.screened_cable, -1.),
            ("tube", ConductorMeshFixtures.tubular_cable, -1.),
            ("sector", ConductorMeshFixtures.sector_cable, 1.))
        design = constructor()
        system = build(LineCableSystem, design, Pose2(0., height);
            connections=Dict(terminal => i for (i, terminal) in enumerate(design.terminal_order)),
            line_length=1., system_id="graded-$name")
        problem = LineParametersProblem(system; frequencies, temperature=20.,
            earth_props=homogeneous(rho=100., eps_r=1., mu_r=1.))
        directory = joinpath(root, name)
        # The exported bundle belongs to the caller. Restarting preserves edits;
        # choose a new output root when changing the experiment's inputs.
        entry = joinpath(directory, "detached", "study.pro")
        isfile(entry) || export_data(:onelab, problem, formulation;
            file_name=entry, mesh_options)
        fixture_case(joinpath(directory, "normal", String(physics)),
            problem, formulation, mesh_options)
        if !isempty(refined_frequencies)
            fine = LineParametersProblem(system; frequencies=refined_frequencies,
                temperature=20., earth_props=problem.earth_props)
            fixture_case(joinpath(directory, "refined", String(physics)),
                fine, formulation, (; mesh_options..., refinements...))
        end
    end
    println("COMPLETE prescribed cable fixture comparisons"); flush(stdout)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_cable_fixtures(abspath(only(ARGS)))
