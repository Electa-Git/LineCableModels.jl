@testitem "Gauntlet / reviewed reference replay is read-only and detects matrix changes" tags=[:gauntlet_toolkit] setup=[TestFixtures] begin
    using .TestFixtures
    using JLD2, SHA, TOML
    using Pkg.Artifacts
    include(joinpath(pkgdir(LineCableModels), "test", "numerical", "references.jl"))
    owner = NumericalReferences
    problem = TestFixtures.line_parameters_problem(; frequencies=[0.1, 50.0, 1.0e6])
    formulation = Formulation(insulation_admittance=:Ametani2004,
        options=(reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
    parameters = compute(problem, formulation)
    # These computed arrays exercise the file/replay protocol, not scientific
    # accuracy. No generated fixture is an approved numerical baseline.
    document = (
        kind=:gauntlet_calculation, status=:complete, backend=:coaxial,
        numerical_reference_approval=:unreviewed,
        problem=LineCableModels.ImportExport.serialize_value(problem),
        formulation=(definitions=formulation.definitions, options=formulation.options),
        Z=parameters.Z.values, Y=parameters.Y.values, frequencies=problem.frequencies,
        port_order=string.(axes(parameters.Z, 1)), basis=:pul, domain=:PhaseDomain,
    )
    protected = [joinpath(pkgdir(LineCableModels), "test", "numerical", name)
        for name in ("approved.toml", "Artifacts.toml")]
    push!(protected, joinpath(pkgdir(LineCableModels), "test", "gauntlet", "Artifacts.toml"))
    original = [(read(path), stat(path).mtime) for path in protected]
    mktempdir() do root
        path = joinpath(root, "reference.jld2")
        JLD2.jldsave(path; document...)
        digest = bytes2hex(open(sha256, path))
        reference = owner.read_reference(path, digest)
        @test LineCableModels.ImportExport.serialize_value(reference.problem) == document.problem
        @test reference.formulation.definitions == formulation.definitions
        @test reference.formulation.options == formulation.options
        @test reference.parameters.Z.values == parameters.Z.values
        @test reference.parameters.Y.values == parameters.Y.values
        comparison = owner.compare_reference(reference)
        @test all(iszero, comparison.Z.absolute)
        @test all(iszero, comparison.Y.absolute)
        @test !isdefined(owner, :GauntletSupport)
        @test !isdefined(owner, :Gmsh)
        @test !isdefined(owner, :PSCADBenchmarks)

        hash = create_artifact() do directory
            cp(path, joinpath(directory, "reference.jld2"))
            altered = copy(document.Z)
            altered[1, 2, :] .+= 1
            JLD2.jldsave(joinpath(directory, "changed.jld2");
                merge(document, (Z=altered,))...)
        end
        bindings = joinpath(root, "Artifacts.toml")
        bind_artifact!(bindings, "protocol_fixture", hash)
        manifest = joinpath(root, "approved.toml")
        entry = Dict{String, Any}(
            "id"=>"protocol_fixture", "artifact"=>"protocol_fixture",
            "file"=>"reference.jld2", "sha256"=>digest,
            "review"=>"Temporary protocol fixture, not a scientific reference approval",
            "Z_atol"=>0.0, "Y_atol"=>0.0, "rtol"=>0.0)
        function write_manifest(entries)
            open(manifest, "w") do io
                TOML.print(io, Dict("schema_version"=>1, "references"=>entries))
            end
        end
        write_manifest([entry])
        pinned = joinpath(artifact_path(hash), "reference.jld2")
        retained = [(read(file), stat(file).mtime) for file in (path, pinned, manifest, bindings)]
        rows = owner.check(manifest)
        @test length(rows) == 2 * prod(size(parameters.Z)[1:2])
        @test all(row -> row.passed && iszero(row.absolute) && iszero(row.relative), rows)
        @test [(read(file), stat(file).mtime) for file in (path, pinned, manifest, bindings)] == retained
        @test JLD2.load(pinned, "numerical_reference_approval") === :unreviewed

        changed = joinpath(artifact_path(hash), "changed.jld2")
        write_manifest([merge(entry, Dict("file"=>"changed.jld2",
            "sha256"=>bytes2hex(open(sha256, changed))))])
        changed_bytes = read(changed)
        changed_mtime = stat(changed).mtime
        rows = owner.check(manifest)
        failure = only(filter(row -> !row.passed, rows))
        @test (failure.quantity, failure.row, failure.column) == (:Z, 1, 2)
        @test failure.absolute ≈ 1.0
        @test read(changed) == changed_bytes
        @test stat(changed).mtime == changed_mtime

        # Missing authority and malformed review settings fail before any replay.
        write_manifest([])
        @test_throws "no numerical references have been approved" owner.check(manifest)
        write_manifest([Dict("review"=>"incomplete")])
        @test_throws "entries require" owner.check(manifest)
        write_manifest([entry, entry])
        @test_throws "duplicate numerical-reference IDs" owner.check(manifest)
        for update in (Dict("review"=>" "), Dict("review"=>1), Dict("rtol"=>NaN),
                Dict("rtol"=>true), Dict("Z_atol"=>-1.0), Dict("Y_atol"=>Inf),
                Dict("artifact"=>"unbound"), Dict("file"=>"../reference.jld2"),
                Dict("file"=>path), Dict("file"=>"missing.jld2"))
            write_manifest([merge(entry, update)])
            @test_throws ArgumentError owner.check(manifest)
        end
        @test_throws "checksum mismatch" owner.read_reference(path, repeat("0", 64))
        @test_throws "reviewed SHA-256" owner.read_reference(path, "unrecorded")
        @test_throws "file is missing" owner.read_reference(path * ".missing", digest)
        for update in ((status=:failed,), (backend=:fem,), (basis=:pu,), (domain=:ModalDomain,),
                (formulation=(definitions=formulation.definitions,),),
                (frequencies=[0.1, 60.0, 1.0e6],), (port_order=fill("duplicate", size(parameters.Z, 1)),),
                (Z=fill(ComplexF64(NaN), size(parameters.Z)),))
            invalid = joinpath(root, "invalid.jld2")
            JLD2.jldsave(invalid; merge(document, update)...)
            @test_throws ArgumentError owner.read_reference(invalid, bytes2hex(open(sha256, invalid)))
        end
        legacy = Dict(pairs(document))
        delete!(legacy, :problem)
        path = joinpath(root, "legacy.jld2")
        JLD2.jldsave(path; legacy...)
        @test_throws "replay requires stored problem" owner.read_reference(path, bytes2hex(open(sha256, path)))
    end
    @test [(read(path), stat(path).mtime) for path in protected] == original
end
