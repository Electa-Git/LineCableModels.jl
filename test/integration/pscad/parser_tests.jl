@testitem "PSCAD / parser and fixed mappings" tags=[:integration, :pscad, :slow] setup=[PSCADStandIn] begin
    pscad = Formulation(:pscad)
    using Test
    import LineCableModels
    using LineCableModels: description
    using LineCableModels.Engine

    harness=LineCableModels.PSCAD
    overhead=Formulation(:pscad; earth_impedance = :gary1976)
    underground=Formulation(:pscad; earth_impedance = :wedepohl1973)
    @test overhead isa harness.PSCADFormulation
    @test underground isa harness.PSCADFormulation
    @test hasmethod(compute, Tuple{LineParametersProblem, harness.PSCADFormulation})
    @test harness.earth_impedance(Formulation(:pscad; earth_impedance=:gary1976).formulas.earth_impedance, Val(:mutual), Val(1), Val(1), pscad) ==
          (EarthForm2 = (value = 0, readback = "DERISEMLYEN"),)
    @test harness.earth_impedance(Formulation(:pscad; earth_impedance=:wedepohl1973).formulas.earth_impedance, Val(:mutual), Val(2), Val(2), pscad) ==
          (EarthForm = (value = 0, readback = "WEDEPOHL"),)
    @test EarthImpedance.formula_id(overhead.formulas.earth_impedance) === :gary1976
    @test EarthImpedance.formula_id(underground.formulas.earth_impedance) ===
          :wedepohl1973
    @test LineCableModels.formula_id(overhead.formulas.earth_admittance) === :ideal
    @test LineCableModels.formula_id(overhead.formulas.insulation_admittance) === :lossless
    @test occursin("lossless", lowercase(description(overhead.formulas.insulation_admittance)))
    # The selected owner supplies descriptions through current computation
    # details. A removed consumer-local label helper is not an asserted requirement.
    @test !isempty(string(computation_details(overhead)))
    @test overhead.options.data == (reduce_bundle = false, kron_reduction = false,
        ideal_transposition = false,base_frequency=50.0)
    @test_throws ArgumentError Formulation(:pscad; options = (output_stem = "invalid",))

    identifiers=(:gary1976, :carson1926, :pollaczek1926, :wedepohl1973,
        :saad1996, :ametani2009, :lucca1994)
    for identifier in identifiers, compact in (false, true)
        selected = Formulation(:pscad; earth_impedance=identifier).formulas.earth_impedance
        @test description(selected; compact) ==
            description(EarthImpedance.Formula{identifier}; compact)
    end
    @test all(
        id -> Formulation(:pscad; earth_impedance = LineCableModels.formula(id)) isa
              harness.PSCADFormulation,
        identifiers)
    @test haskey(harness.earth_impedance(Formulation(:pscad; earth_impedance=:saad1996).formulas.earth_impedance, Val(:mutual), Val(2), Val(2), pscad), :EarthForm)
    @test harness.earth_impedance(Formulation(:pscad; earth_impedance=:lucca1994).formulas.earth_impedance, Val(:mutual), Val(1), Val(2), pscad).EarthForm3.readback == "LUCCA"
    @test harness.earth_impedance(Formulation(:pscad; earth_impedance=:carson1926).formulas.earth_impedance, Val(:mutual), Val(1), Val(1), pscad).EarthForm2.value == 2
    @test harness.earth_impedance(Formulation(:pscad; earth_impedance=:pollaczek1926).formulas.earth_impedance, Val(:mutual), Val(2), Val(2), pscad).EarthForm.value == 2
    @test_throws ArgumentError harness.earth_impedance(Formulation(:pscad; earth_impedance=:pollaczek1926).formulas.earth_impedance, Val(:mutual), Val(1), Val(1), pscad)

    mktempdir() do directory
        frequency=[10.0, 100.0]
        magnitude=[.1*(1+2k+3i) for k in 1:2,i in 1:4]
        phase=[5.0*(2k+i) for k in 1:2,i in 1:4]
        ymag=[1e-6*(7+3k+i) for k in 1:2,i in 1:4]
        yphase=[-3.0*(k+2i) for k in 1:2,i in 1:4]
        function table(path, values)
            open(path, "w") do io
                println(io, " LOG10(FN) FN V11 V12 V21 V22")
                for index in eachindex(frequency)
                    println(
                        io,
                        join(
                            (log10(frequency[index]),
                                frequency[index], values[index, :]...), ' ')
                    )
                end
            end
        end
        table(joinpath(directory, "result_zm.out"), magnitude)
        table(joinpath(directory, "result_zp.out"), phase)
        table(joinpath(directory, "result_ym.out"), ymag)
        table(joinpath(directory, "result_yp.out"), yphase)
        result=harness.read_pscad_result(directory, frequency, (2, 2, 2), 1.0)
        @test result isa LineParameters
        @test size(Z(result)) == (2, 2, 2)
        @test Z(result)[:, :, 1] ≈ permutedims(reshape(magnitude[1,:].*cispi.(phase[1,:]./180),2,2))
        @test Y(result)[:, :, 2] ≈ permutedims(reshape(ymag[2,:].*cispi.(yphase[2,:]./180),2,2))
        rounded_frequency=[10.0, 100.0*(1+4.0e-8)]
        rounded_result=harness.read_pscad_result(
            directory, rounded_frequency, (2, 2, 2), 1.0
        )
        @test frequencies(rounded_result) == rounded_frequency
        @test_throws ArgumentError harness.read_pscad_result(
            directory, [1.0, 11.0], (2, 2, 2), 1.0
        )
        mv(joinpath(directory, "result_zm.out"), joinpath(directory, "raw_zm.out"))
        @test_throws ArgumentError harness.read_pscad_result(directory, frequency, (
            2, 2, 2), 1.0)
    end

    config=harness.RemoteConfig("host", 54330, raw"Z:\pscad\benchmarks\.work",
        raw"C:\PSCAD\Pscad.exe"; local_root=mktempdir(), timeout = 60)
    owner=harness.PSCADFormulation
    @test LineCableModels.formulation_options(owner, FormulationOptions()) == overhead.options
    @test_throws ArgumentError LineCableModels.formulation_options(
        owner,
        FormulationOptions(output_stem = "case")
    )
    @test_throws ArgumentError LineCableModels.computation_options(owner, ComputationOptions((;)))
    execution=LineCableModels.computation_options(owner, ComputationOptions((remote = config,)))
    # The verbosity and the callable of a computation are options of its run, not of PSCAD.
    for run_option in ((verbosity = (default = 0, PSCAD = 2),), (on_result = identity,))
        @test_throws ArgumentError LineCableModels.computation_options(owner,
            ComputationOptions(merge((; remote = config), run_option)))
    end
    @test_throws ArgumentError LineCableModels.computation_options(owner,
        ComputationOptions(; remote=config, output_basis=:total))
    for timing in (true, false, 1, nothing, :yes)
        @test_throws ArgumentError LineCableModels.computation_options(owner,
            ComputationOptions(; remote=config, timing))
    end
    @test execution.data == (
        work_root=config.local_root,
        remote = config,
        resume_run_directory = nothing,
        solver_identity = nothing
    )
    @test LineCableModels.computation_options(owner, execution) === execution
    # A run's directory gives the name of its line and of the files that PSCAD writes.
    @test_throws ArgumentError LineCableModels.computation_options(
        owner, ComputationOptions((output_stem = "case", remote = config)))
    root=config.local_root
    @test LineCableModels.computation_options(owner, ComputationOptions((remote=config,work_root=joinpath(root,"any","depth")))).data.work_root == joinpath(root,"any","depth")
    @test_throws ArgumentError LineCableModels.computation_options(owner, ComputationOptions((remote=config,work_root=dirname(root))))
    @test_throws ArgumentError harness.RemoteConfig(
        "host", 54330, "shared", "Pscad.exe"; local_root=mktempdir(), timeout = 0
    )
    @test_throws MethodError harness.RemoteConfig(
        "host", 54330, "shared", "Pscad.exe"; local_root=mktempdir(), verbosity = 2
    )
    @test_throws ArgumentError validate(
        collect(10.0 .^ range(0, stop = 6, length = 61)), harness.PSCADFormulation
    )
    frequency_probe=collect(10.0 .^ range(0, stop = 6, length = 101))
    @test validate(frequency_probe, harness.PSCADFormulation) === frequency_probe
    minimum_frequency_probe=collect(10.0 .^ range(-1, stop = 6, length = 101))
    @test validate(minimum_frequency_probe, harness.PSCADFormulation) ===
          minimum_frequency_probe
    @test_throws DomainError validate(
        collect(10.0 .^ range(-2, stop = 6, length = 101)), harness.PSCADFormulation
    )

    extension=Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    @test extension._remote_path(raw"Z:\pscad\exchange\\", "run-a", "case.pscx") ==
          raw"Z:\pscad\exchange\run-a\case.pscx"
    mktempdir() do directory
        partial=joinpath(directory, "fixture_zm.out")
        write(partial, "header\n")
        writer=@async begin
            sleep(0.05)
            open(partial, "a") do io
                println(io, "0.0 1.0 2.0")
                println(io, "1.0 10.0 3.0")
            end
        end
        completed=extension._wait_output(partial, 2; timeout = 1, poll_interval = 0.01)
        wait(writer)
        @test completed == partial
        @test extension._data_rows(completed) == 2
        open(partial, "a") do io
            println(io, "2.0 100.0 4.0")
        end
        @test_throws r"contains 3 rows; expected 2" extension._wait_output(partial, 2;
            timeout = 1, poll_interval = 0.01)
    end
    mktempdir() do directory
        partial=joinpath(directory, "fixture_zm.out")
        write(partial, "header\n")
        error=try
            extension._wait_output(partial, 2; timeout = 0.05, poll_interval = 0.01)
            nothing
        catch caught
            caught
        end
        @test error isa ArgumentError
        @test occursin("0 of 2 rows", sprint(showerror, error))
        @test_throws r"0 of 2 rows" extension._wait_output(joinpath(directory, "absent_zm.out"),
            2; timeout = 0.05, poll_interval = 0.01)
    end

    for count in (101,201,501,1001), T in (Float32,Float64)
        f=T.(10.0.^range(-1,6;length=count))
        @test validate(f, harness.PSCADFormulation) === f
    end
    @test_throws ArgumentError validate(collect(range(.1,1e6;length=201)), harness.PSCADFormulation)
end
