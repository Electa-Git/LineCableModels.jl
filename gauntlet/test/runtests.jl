using Gauntlet
using LineCableModels
using LinearAlgebra
using Tables
using Test
using TOML
import JLD2

const IO = Gauntlet.PSCADIO

struct FixtureDatasource <: Datasource end
Gauntlet.datasource(::Val{:fixture}) = FixtureDatasource()
Gauntlet.decode(::FixtureDatasource, record) = (; decoded = record)
Gauntlet.load(source::FixtureDatasource, record) = decode(source, record)
Gauntlet.ingest(::FixtureDatasource, source, destination) = (; source, destination)

function temporary_file(test, content, extension)
    mktemp() do path, stream
        close(stream)
        target = path * extension
        open(target, "w") do output
            write(output, content)
        end
        try
            test(target)
        finally
            rm(target; force = true)
        end
    end
end

@testset "Datasource grammar" begin
    @test datasource(:pscad) isa PSCAD
    @test datasource(:fem) isa FEM
    @test datasource(:fixture) isa FixtureDatasource
    @test_throws ArgumentError datasource(:unknown)
    @test !isdefined(Gauntlet, :Source)
    @test !isdefined(Gauntlet, :Corpus)

    dataset = Dataset(:smoke)
    @test dataset.datasource isa PSCAD
    suite_fields = fieldnames(Base.unwrap_unionall(Suite))
    @test :dataset in suite_fields
    @test !(:corpus in suite_fields)

    id = first(keys(dataset.cases))
    path = joinpath(dataset.root, dataset.cases[id])
    @test load(:pscad, path).id == id
    record = Gauntlet.PSCADCase(JLD2.load(path, "record"))
    @test decode(:pscad, record) isa Reference
    @test load(:fixture, :record) == (; decoded = :record)
    @test decode(:fixture, :record) == (; decoded = :record)
    @test ingest(:fixture, :source, :destination) ==
          (; source = :source, destination = :destination)

    for action in (
        () -> ingest(:fem, "source", "destination"),
        () -> load(:fem, "record"),
        () -> decode(:fem, "record")
    )
        exception = try
            action()
            nothing
        catch error
            error
        end
        @test exception isa ErrorException
        @test occursin("FEM datasource is not implemented", sprint(showerror, exception))
    end

    mktempdir() do root
        open(joinpath(root, "index.toml"), "w") do io
            TOML.print(io, Dict("schema_version" => 1))
        end
        @test_throws ArgumentError Dataset(root)
        @test_throws ArgumentError ingest(:pscad, root, joinpath(root, "output"))
    end

    generic_dataset = read(joinpath(@__DIR__, "..", "src", "dataset.jl"), String)
    @test !occursin("PSCADIO", generic_dataset)
    @test !occursin(r"\.(cli|tli|clo|tlo)", lowercase(generic_dataset))
end

@testset "PSCAD input grammar" begin
    source = """
    Cable Summary:
    {
      Cable Name = sample
      Cable Length = 100.0
    }
    Frequency Dep. (Phase) Model Options:
    {
      Curve Fitting Start Frequency = 0.001
      Curve Fitting End Frequency = 1000.0
      Total Number of Frequency Increments = 6
    }
    Coax Cable:
    {
      P1 = 0.0 1.0
      Layers = 1
    }
    """
    temporary_file(source, ".cli") do path
        parsed = IO.read_input(path)
        @test parsed.kind === :coax
        @test only(IO.block(parsed.root, "Cable Summary")).fields["Cable Length"] == 100.0
    end
    temporary_file("{\n", ".cli") do path
        @test_throws ArgumentError IO.read_input(path)
    end
    temporary_file("immutable evidence\n", ".out") do path
        digest = Gauntlet._sha256(path)
        @test Gauntlet._verify(path, digest) == path
        @test_throws ArgumentError Gauntlet._verify(path, repeat("0", 64))
    end
end

@testset "Detailed polar matrices" begin
    magnitude = """
      LOG10(FN) FN Calc. Mag.: Z11, Z12 ...[ohm/m]
      0.0 1.0 1.0 2.0 3.0 4.0
      1.0 10.0 2.0 3.0 4.0 5.0
    """
    phase = """
      LOG10(FN) FN Calc. Phase: Z11, Z12 ...[degree]
      0.0 1.0 0.0 90.0 180.0 -90.0
      1.0 10.0 0.0 0.0 0.0 0.0
    """
    temporary_file(magnitude, ".out") do magnitude_path
        temporary_file(phase, ".out") do phase_path
            observed = IO.combine_polar(
                IO.read_detailed(magnitude_path), IO.read_detailed(phase_path)
            )
            @test size(observed) == (2, 2, 2)
            @test observed[1, 2, 1] ≈ 2im
            @test observed[2, 1, 1] ≈ -3.0
        end
    end

    overflow = """
      LOG10(FN) FN VAOPEN(M,PHI)
      0.0 1.0 0.15942974-142 0.0
    """
    temporary_file(overflow, ".out") do path
        parsed = IO.read_detailed(path)
        @test parsed.values[1, 1] == 0.15942974e-142
        @test IO.terminal_values(parsed)[1, 1] == 0.15942974e-142
    end
    temporary_file("LOG10(FN) FN empty\n", ".out") do path
        @test IO.read_detailed(path; allow_empty = true) === nothing
        @test_throws ArgumentError IO.read_detailed(path)
    end
end

@testset "Ordinary phase and sequence matrices" begin
    ordinary = """
      PHASE DOMAIN DATA @ 60.000 Hz:
      SERIES IMPEDANCE MATRIX (Z) [ohms/m]:
      1.0,2.0 3.0,4.0
      3.0,4.0 5.0,6.0

      SHUNT ADMITTANCE MATRIX (Y) [mhos/m]:
      1e-9,2e-9 0.0,0.0
      0.0,0.0 3e-9,4e-9

      SEQUENCE TRANSFORM MATRIX:
      1.0,0.0 0.0,0.0
      0.0,0.0 1.0,0.0

      SEQUENCE IMPEDANCE MATRIX (Zsq) [ohms/m]:
      1.0,2.0 3.0,4.0
      3.0,4.0 5.0,6.0

      SEQUENCE ADMITTANCE MATRIX (Ysq) [mhos/m]:
      1e-9,2e-9 0.0,0.0
      0.0,0.0 3e-9,4e-9
    """
    temporary_file(ordinary, ".out") do path
        parsed = IO.read_ordinary(path)
        @test parsed.frequency == 60.0
        @test size(parsed.Z) == (2, 2)
        @test parsed.sequence_transform == Matrix{ComplexF64}(I, 2, 2)
    end
end

@testset "CLO and TLO vector fits" begin
    fit_text = """
      Fitting parameters for the char. admittance:
      N.o. conductors:
      1
      Residues and poles:
      1
      -2.0
      0.0
      4.0
      0.0
      0.5
      Fitting parameters for the propagation function:
      1
      0.01
      1
      -3.0
      0.0
      6.0
      0.0
    """
    for extension in (".clo", ".tlo")
        temporary_file(fit_text, extension) do path
            fit = IO.read_fit(path; frequency_range = (1.0, 100.0))
            value = evaluate(fit, [1.0, 10.0])
            @test size(value.characteristic) == (1, 1, 2)
            @test size(value.propagation) == (1, 1, 2)
            @test all(isfinite, value.characteristic)
        end
    end
    temporary_file(first(fit_text, 120), ".clo") do path
        @test_throws EOFError IO.read_fit(path)
    end
end

@testset "Assignment and explicit port alignment" begin
    scores = Matrix{Float64}(I, 24, 24)
    @test Gauntlet._assignment(scores) == collect(1:24)

    frequency = [60.0]
    actual_values = reshape(ComplexF64[1, 2, 3, 4], 2, 2, 1)
    expected_values = actual_values[[2, 1], [2, 1], :]
    actual = LineCableModels.Engine.LineParameters(
        actual_values, actual_values, frequency; basis = :per_length
    )
    expected = LineCableModels.Engine.LineParameters(
        expected_values, expected_values, frequency; basis = :per_length
    )
    ports = [Port("second", 1, "conductor", 2), Port("first", 2, "conductor", 1)]
    comparison = compare(
        actual,
        expected,
        MatrixCheck{:Z}(),
        Tolerance(rtol = 1e-6, atol = 0.0);
        ports
    )
    @test comparison.verdict isa Gauntlet.Pass
    @test_throws ArgumentError compare(
        actual,
        expected,
        MatrixCheck{:Z}(),
        Tolerance(rtol = 1e-6, atol = 0.0);
        ports = [ports[1], ports[1]]
    )

    target = repeat(reshape(ComplexF64[1 0; 0 1], 2, 2, 1), 1, 1, 2)
    observed = similar(target)
    observed[:, :, 1] .= target[:, [2, 1], 1] * Diagonal(ComplexF64[im, -1])
    observed[:, :, 2] .= target[:, [2, 1], 2] * Diagonal(ComplexF64[-im, 1])
    modal_frequency = [1.0, 10.0]
    reference_modes = Modes(
        modal_frequency;
        transform = target,
        propagation = ComplexF64[1 1; 2 2]
    )
    actual_modes = Modes(
        modal_frequency;
        transform = observed,
        propagation = ComplexF64[2 2; 1 1]
    )
    alignment = Gauntlet._modal_alignment(
        actual_modes, reference_modes, Tolerance(rtol = 1e-6, atol = 1e-12)
    )
    @test alignment.aligned ≈ target
    @test all(==(Int[2, 1]), alignment.assignments)
    metrics = Gauntlet._modal_metrics(
        reference_modes.propagation,
        reference_modes.propagation,
        modal_frequency,
        Tolerance(rtol = 1e-6, atol = 0.0)
    )
    @test iszero(metrics.max_rel)
    @test iszero(metrics.symmetry)

    native = Gauntlet._native_modes(actual, 1.0)
    @test size(native.transform) == (2, 2, 1)
    @test size(native.propagation) == (2, 1)
    @test native.phase_propagation !== nothing

    sequence = ComplexF64[1 2 3; 4 5 6; 7 8 9]
    sequence_ports = [Port("phase:$index", index, "circuit:1:1", index) for index in 1:3]
    append!(sequence_ports, [Port("ground:$index", index, "ground:1:1:$index", index)
                             for index in 4:5])
    expanded = Gauntlet._sequence_matrix(sequence, sequence_ports, 5)
    @test expanded[1:3, 1:3] == sequence
    @test expanded[4:5, 4:5] == Matrix{ComplexF64}(I, 2, 2)
    @test iszero(expanded[1:3, 4:5])
end

@testset "Kron cohort diagnostics" begin
    frequency = [60.0]
    z = reshape(ComplexF64[2 0.5; 0.5 1], 2, 2, 1)
    y = reshape(ComplexF64[3 0.25; 0.25 2], 2, 2, 1)
    retained_parameters = LineCableModels.Engine.LineParameters(
        z, y, frequency; basis = :per_length
    )
    reduced_parameters = LineCableModels.Engine.LineParameters(
        Gauntlet._kron(z, [1], :Z),
        Gauntlet._kron(y, [1], :Y),
        frequency;
        basis = :per_length
    )
    identity = repeat("0", 64)
    source = Provenance(
        :pscad,
        "5.1.0",
        "synthetic",
        "cohort",
        identity,
        identity,
        identity,
        nothing,
        Dict{String, String}(),
        "shared-cohort"
    )
    retained_ports = [Port("a", 1, "a", 1), Port("b", 1, "b", 2)]
    reduced_ports = [Port("a", 1, "a", 1)]
    retained_case = Case(
        "retained",
        Coax(),
        Exact(),
        Retained(),
        nothing,
        nothing,
        Reference(phase = retained_parameters, ports = retained_ports),
        (),
        Assumption[],
        source
    )
    reduced_case = Case(
        "reduced",
        Coax(),
        Exact(),
        Reduced(),
        nothing,
        nothing,
        Reference(phase = reduced_parameters, ports = reduced_ports),
        (),
        Assumption[],
        source
    )
    trials = Trial[
        Trial(retained_case, retained_parameters, Comparison[], nothing),
        Trial(reduced_case, reduced_parameters, Comparison[], nothing)
    ]
    Gauntlet._append_kron!(trials)
    @test length(trials[2].comparisons) == 2
    @test all(comparison -> comparison.verdict isa Pass, trials[2].comparisons)
    @test Set(Gauntlet._check_name(comparison.check)
    for comparison in trials[2].comparisons) ==
          Set(("physical.kron_Z", "physical.kron_Y"))
end

@testset "Tracked smoke dataset" begin
    dataset = Dataset(:smoke)
    @test length(dataset) == 16
    @test length(keys(dataset)) == 16
    @test all(case -> case isa Case, dataset)
    @test count(case -> case.fidelity isa Rejected, dataset) == 1
    @test Set(nameof(typeof(case.family)) for case in dataset) ==
          Set((:Coax, :Overhead, :Mixed, :Pipe))
    rejected = only(case for case in dataset if case.fidelity isa Rejected)
    trial = gauntlet(rejected)
    @test only(trial.comparisons).verdict isa Gauntlet.ReferenceRejected
    suite = Suite(
        :typed;
        dataset,
        ids = [rejected.id],
        policy = AllCases(),
        performance = String[]
    )
    @test suite.policy isa AllCases
    @test isempty(suite.performance)
end

@testset "Comparisons, fits, reports, and archive identity" begin
    frequencies = [1.0, 10.0]
    values = reshape(ComplexF64[1 + im, 2 + 2im], 1, 1, 2)
    parameters = LineCableModels.Engine.LineParameters(
        values, values, frequencies; basis = :per_length
    )
    comparison = compare(
        parameters, parameters, MatrixCheck{:Z}(), Tolerance(rtol = 1e-6, atol = 0.0)
    )
    @test comparison.verdict isa Gauntlet.Pass
    @test iszero(comparison.metrics.max_abs)
    @test Gauntlet._check_name(comparison.check) == "matrix.Z"
    @test_throws DomainError Tolerance(rtol = 2e-3, atol = 0.0)
    @test Gauntlet._round125(1.01e-5) == 2e-5
    @test Gauntlet._round125(2.01e-5) == 5e-5
    @test 0 <= Gauntlet._calibration_bucket(repeat("a", 64)) <= 4
    source_frequency = [1.0, 10.0]
    source_values = reshape(ComplexF64[1, 10], 1, 1, 2)
    interpolated = Gauntlet._polar_log_interpolate(
        source_values, source_frequency, sqrt(10.0)
    )
    @test only(interpolated) ≈ sqrt(10.0)

    dataset = Dataset(:smoke)
    rejected = only(case for case in dataset if case.fidelity isa Rejected)
    suite = Suite(:parser; dataset, ids = [rejected.id])
    report = gauntlet(suite)
    @test length(collect(Tables.rows(report))) == 1
    mktempdir() do destination
        for extension in ("json", "csv", "md")
            path = write_report(joinpath(destination, "report.$extension"), report)
            @test filesize(path) > 0
        end
    end
    @test Gauntlet.tree_hash(joinpath(@__DIR__, "..", "fixtures", "smoke")) ==
          Gauntlet.tree_hash(joinpath(@__DIR__, "..", "fixtures", "smoke"))
end

@testset "Opt-in full artifact accounting" begin
    artifact = get(ENV, "GAUNTLET_ARTIFACT_DIR", "")
    if isempty(artifact)
        @test true
    else
        dataset = Dataset(artifact)
        index = TOML.parsefile(joinpath(artifact, "index.toml"))
        @test length(dataset.cases) == 869
        @test length(dataset.rejections) == 199
        @test length(dataset.aliases) == 36
        @test length(index["excluded"]) == 1
        @test count(case -> case.fidelity isa Rejected, dataset) == 199
    end
end
