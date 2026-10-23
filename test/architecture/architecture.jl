include("guards.jl")

# Loads a probe package from `probes/<directory>`: `src/<name>.jl` and, if given, the
# extension `ext/<extension>.jl`. Returns its tree and its extension modules.
function probe(directory, name; extension = nothing)
    root_directory = joinpath(@__DIR__, "probes", directory)
    root = Base.include(Main, joinpath(root_directory, "src", name * ".jl"))
    extensions = extension === nothing ? Module[] :
                 [Base.include(Main, joinpath(root_directory, "ext", extension * ".jl"))]
    tree = Base.invokelatest(Guards.PackageTree, root, extensions, root_directory)
    return tree, extensions
end

@testset "architecture" begin
    tree = Guards.PackageTree(LineCableModels, Module[], pkgdir(LineCableModels))
    files = Guards.source_files(pkgdir(LineCableModels))
    directory = pkgdir(LineCableModels)
    methods = Guards.package_methods(tree)
    homes = Guards.home_directories(tree)

    @testset "ownership" begin
        @test Guards.ownership(methods, tree) == Dict{String, Int}()
    end

    @testset "placement" begin
        @test Guards.placement(methods, tree, homes) == Dict{String, Int}()
    end

    @testset "names" begin
        @test Guards.shared_names(tree) == Dict{String, Int}()
    end

    @testset "shadowing" begin
        @test Guards.shadowing(tree, Guards.DEPENDENCIES) == Dict{String, Int}()
    end

    @testset "explicit imports" begin
        @test Guards.explicit_imports(tree, Module[]) == Dict{String, Int}()
    end

    @testset "direction" begin
        @test Set(nameof(M)
        for M in tree.core
        if M !== tree.root && parentmodule(M) === tree.root) ⊆ Set(Guards.ORDER)
        @test Guards.direction(methods, tree, Guards.ORDER) == Dict{String, Int}()
    end

    @testset "validate returns its subject" begin
        @test Guards.validate_returns(files, directory) == Dict{String, Int}()
    end

    @testset "reserved verbs" begin
        @test Guards.reserved_verbs(files, directory) == Dict{String, Int}()
    end

    @testset "storage vocabulary" begin
        @test Guards.storage_vocabulary(files, directory) == Dict{String, Int}()
    end

    @testset "repeated access names" begin
        @test Guards.repeated_access_names(files, directory) == Dict{String, Int}()
    end

    @testset "symbol switches and probes" begin
        @test Guards.symbol_switches(files, directory) == Dict{String, Int}()
    end

    @testset "decisions by absence" begin
        @test Guards.absences(files, directory) == Dict{String, Int}()
    end

    @testset "file hygiene" begin
        @test Guards.file_hygiene(files, directory) == Dict{String, Int}()
    end

    @testset "negative controls" begin
        planted, _ = probe("structure/planted", "ArchitectureProbePlanted")
        planted_methods = Base.invokelatest(Guards.package_methods, planted)
        planted_homes = Base.invokelatest(Guards.home_directories, planted)
        @test Base.invokelatest(Guards.ownership, planted_methods, planted) ==
              Dict{String, Int}(
            "Consumer | Early.extend | src/consumer/Consumer.jl" => 1,
            "Consumer | Early.Item | src/consumer/Consumer.jl" => 1)
        @test Base.invokelatest(Guards.placement, planted_methods, planted, planted_homes) ==
              Dict{String, Int}(
            "src/early/stray.jl | Consumer" => 1,
            "ext/core.jl | Consumer" => 1)
        @test Base.invokelatest(Guards.direction, planted_methods, planted,
            (:Early, :Late, :Consumer)) ==
              Dict{String, Int}("Early -> Late | src/early/Early.jl" => 1)
        @test Base.invokelatest(Guards.shared_names, planted) ==
              Dict{String, Int}("quantity" => 2)
        @test Base.invokelatest(Guards.shadowing, planted, (Base,)) ==
              Dict{String, Int}("Consumer.filter -> Base" => 1)

        clean,
        _ = probe("structure/clean", "ArchitectureProbeClean";
            extension = "ArchitectureProbeCleanExt")
        clean_methods = Base.invokelatest(Guards.package_methods, clean)
        clean_homes = Base.invokelatest(Guards.home_directories, clean)
        @test Base.invokelatest(Guards.ownership, clean_methods, clean) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.placement, clean_methods, clean, clean_homes) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.direction, clean_methods, clean, (:Early, :Late)) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.shared_names, clean) == Dict{String, Int}()
        @test Base.invokelatest(Guards.shadowing, clean, (Base,)) == Dict{String, Int}()

        planted,
        planted_extensions = probe("imports/planted", "ImportsProbePlanted";
            extension = "ImportsProbePlantedExt")
        @test Base.invokelatest(Guards.explicit_imports, planted, planted_extensions) ==
              Dict{String, Int}(
            "no_implicit_imports | Consumer | Exporter" => 1,
            "no_implicit_imports | Consumer | implicit_value" => 1,
            "no_stale_explicit_imports | Consumer | public_value" => 1,
            "all_explicit_imports_via_owners | Consumer | relayed_value" => 1,
            "all_explicit_imports_are_public | Consumer | relayed_value" => 1,
            "all_explicit_imports_are_public | Consumer | private_import" => 1,
            "all_qualified_accesses_via_owners | Consumer | relayed_value" => 1,
            "all_qualified_accesses_are_public | Consumer | relayed_value" => 1,
            "all_qualified_accesses_are_public | Consumer | private_value" => 1,
            "all_qualified_accesses_are_public | ImportsProbePlantedExt | private_value" => 1,
            "no_self_qualified_accesses | Consumer | implicit" => 1)

        clean,
        clean_extensions = probe("imports/clean", "ImportsProbeClean";
            extension = "ImportsProbeCleanExt")
        @test Base.invokelatest(Guards.explicit_imports, clean, clean_extensions) ==
              Dict{String, Int}()

        directory = joinpath(@__DIR__, "probes", "sources", "planted")
        files = Guards.source_files(directory)
        @test Guards.validate_returns(files, directory) ==
              Dict{String, Int}("src/validate.jl" => 6)
        @test Guards.reserved_verbs(files, directory) == Dict{String, Int}(
            "src/validate.jl | __validate_shape" => 1,
            "src/validate.jl | _check_input" => 1,
            "src/validate.jl | ensure_ready" => 1,
            "src/validate.jl | require_kind" => 1)
        @test Guards.storage_vocabulary(files, directory) == Dict{String, Int}(
            "src/storage.jl | AbstractStorage" => 1,
            "src/storage.jl | FormulaBuffers" => 1,
            "src/storage.jl | KernelCache" => 1,
            "src/storage.jl | ProbeWorkspace" => 2,
            "src/storage.jl | cache" => 2,
            "src/storage.jl | kernel_workspace" => 1,
            "src/storage.jl | sample_storage" => 1,
            "src/storage.jl | scratch" => 1,
            "src/storage.jl | storage" => 1,
            "src/storage.jl | trace_buffers" => 1,
            "src/storage.jl | work" => 2)
        @test Guards.repeated_access_names(files, directory) ==
              Dict{String, Int}("src/repeats.jl" => 3)
        @test Guards.symbol_switches(files, directory) == Dict{String, Int}(
            "src/switches.jl | applicable" => 2,
            "src/switches.jl | eval" => 2,
            "src/switches.jl | symbol" => 9,
            "src/switches.jl | type test" => 9)
        @test Guards.absences(files, directory) == Dict{String, Int}(
            "src/absence.jl | argument" => 9,
            "src/absence.jl | comparison" => 9,
            "src/absence.jl | ismissing" => 1,
            "src/absence.jl | isnothing" => 4,
            "src/absence.jl | parameter" => 7,
            "src/absence.jl | something" => 1)
        @test Guards.file_hygiene(files, directory) == Dict{String, Int}(
            "ext/hygiene.jl | homedir" => 1,
            "src/hygiene.jl | homedir" => 1,
            "src/hygiene.jl | mktempdir" => 1,
            "src/hygiene.jl | tempdir" => 1,
            "src/hygiene.jl | tempname" => 1)

        directory = joinpath(@__DIR__, "probes", "sources", "clean")
        files = Guards.source_files(directory)
        @test Guards.validate_returns(files, directory) == Dict{String, Int}()
        @test Guards.reserved_verbs(files, directory) == Dict{String, Int}()
        @test Guards.storage_vocabulary(files, directory) == Dict{String, Int}()
        @test Guards.repeated_access_names(files, directory) == Dict{String, Int}()
        @test Guards.symbol_switches(files, directory) == Dict{String, Int}()
        @test Guards.absences(files, directory) == Dict{String, Int}()
        @test Guards.file_hygiene(files, directory) == Dict{String, Int}()
    end
end
