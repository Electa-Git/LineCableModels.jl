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
    called = Guards.called_names(keys(Guards.test_sources(
        joinpath(directory, "test"); excluded = Guards.NOT_TESTS)))
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

    @testset "Commons admission" begin
        @test Guards.commons_admission(methods, tree, homes) == Dict{String, Int}()
    end

    @testset "reserved vocabulary" begin
        @test Guards.reserved_vocabulary(files, directory, Guards.VOCABULARY) ==
              Dict{String, Int}()
    end

    @testset "literal fingerprints" begin
        @test Guards.fingerprints(files, directory) == Dict{String, Int}()
    end

    @testset "clones of Commons bodies" begin
        @test Guards.clones(files, methods, tree) == Dict{String, Int}()
    end

    @testset "tiny private helpers" begin
        @test Guards.tiny_helpers(methods, tree, called) == Dict{String, Int}()
    end

    @testset "root module" begin
        @test Guards.root_definitions(tree) == Dict{String, Int}()
    end

    @testset "boxed captures" begin
        @test Guards.boxed_captures(methods, tree) == Dict{String, Int}()
    end

    @testset "typed values" begin
        @test Guards.typed_values(tree) == Dict{String, Int}()
    end

    @testset "naming rule" begin
        names = Guards.TOML.parsefile(joinpath(@__DIR__, "names.toml"))
        @test Guards.naming_rule(methods, tree, names) == Dict{String, Int}()
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

        vocabulary = Dict{Symbol, Pair{String, Regex}}(
            :ideal_transposition! => "src/commons" => r"transpos",
            :vacuum_permittivity => "src/commons" => r"(?i)^_*(\w+_)?(vacuum_permittivity|eps(ilon)?_?[0₀]|[εϵ]_?[0₀])$",
            :vacuum_permeability => "src/commons" => r"(?i)^_*(\w+_)?(vacuum_permeability|mu_?[0₀]|[μµ]_?[0₀])$")
        planted,
        _ = probe("commons/planted", "CommonsProbePlanted";
            extension = "CommonsProbePlantedExt")
        planted_methods = Base.invokelatest(Guards.package_methods, planted)
        planted_homes = Base.invokelatest(Guards.home_directories, planted)
        planted_files = Guards.source_files(planted.directory)
        planted_called = Guards.called_names(keys(Guards.test_sources(
            joinpath(planted.directory, "test"); excluded = Guards.NOT_TESTS)))
        @test Base.invokelatest(Guards.commons_admission, planted_methods, planted, planted_homes) ==
              Dict{String, Int}(
            "Lone | owners" => 1,
            "extension_shared | owners" => 1,
            "hidden | docstring" => 1,
            "hidden | owners" => 1,
            "hidden | public" => 1,
            "hidden | tests" => 1,
            "hidden_wrapper | docstring" => 1,
            "hidden_wrapper | public" => 1,
            "hidden_wrapper | tests" => 1,
            "inner_only | owners" => 1,
            "lonely | owners" => 1,
            "undocumented | docstring" => 1,
            "untested | tests" => 1,
            "welcome | tests" => 1)
        @test Base.invokelatest(Guards.reserved_vocabulary, planted_files, planted.directory, vocabulary) ==
              Dict{String, Int}(
            "src/alpha/Alpha.jl | EPS0" => 1,
            "src/alpha/Alpha.jl | transposed!" => 1,
            "src/alpha/Alpha.jl | μ0" => 1)
        @test Base.invokelatest(Guards.fingerprints, planted_files, planted.directory) ==
              Dict{String, Int}("src/alpha/Alpha.jl" => 2, "src/beta/Beta.jl" => 2)
        @test Base.invokelatest(Guards.clones, planted_files, planted_methods, planted) ==
              Dict{String, Int}("src/beta/Beta.jl | average_cyclic! ~ ideal_transposition!" => 1)
        @test Base.invokelatest(Guards.tiny_helpers, planted_methods, planted, planted_called) ==
              Dict{String, Int}("src/alpha/Alpha.jl" => 2, "src/beta/Beta.jl" => 1)
        @test Base.invokelatest(Guards.root_definitions, planted) ==
              Dict{String, Int}("CommonsProbePlanted" => 2)

        clean,
        _ = probe("commons/clean", "CommonsProbeClean";
            extension = "CommonsProbeCleanExt")
        clean_methods = Base.invokelatest(Guards.package_methods, clean)
        clean_homes = Base.invokelatest(Guards.home_directories, clean)
        clean_files = Guards.source_files(clean.directory)
        clean_called = Guards.called_names(keys(Guards.test_sources(
            joinpath(clean.directory, "test"); excluded = Guards.NOT_TESTS)))
        @test Base.invokelatest(Guards.commons_admission, clean_methods, clean, clean_homes) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.reserved_vocabulary, clean_files, clean.directory, vocabulary) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.fingerprints, clean_files, clean.directory) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.clones, clean_files, clean_methods, clean) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.tiny_helpers, clean_methods, clean, clean_called) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.root_definitions, clean) == Dict{String, Int}()

        planted, _ = probe("values/planted", "ValuesProbePlanted")
        planted_methods = Base.invokelatest(Guards.package_methods, planted)
        @test Base.invokelatest(Guards.boxed_captures, planted_methods, planted) ==
              Dict{String, Int}("ValuesProbePlanted | boxed | src/ValuesProbePlanted.jl" => 1)
        @test Base.invokelatest(Guards.typed_values, planted) == Dict{String, Int}(
            "ValuesProbePlanted | Loose | constructor | src/ValuesProbePlanted.jl" => 1,
            "ValuesProbePlanted | Loose | x" => 1,
            "ValuesProbePlanted | Open | constructor | src/ValuesProbePlanted.jl" => 1,
            "ValuesProbePlanted | Open | x" => 1,
            "ValuesProbePlanted | Plain | constructor | src/ValuesProbePlanted.jl" => 3,
            "ValuesProbePlanted | Wide | constructor | src/ValuesProbePlanted.jl" => 1,
            "ValuesProbePlanted | Wide | x" => 1)

        clean, _ = probe("values/clean", "ValuesProbeClean")
        clean_methods = Base.invokelatest(Guards.package_methods, clean)
        @test Base.invokelatest(Guards.boxed_captures, clean_methods, clean) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.typed_values, clean) == Dict{String, Int}()

        planted, _ = probe("names/planted", "NamesProbePlanted")
        planted_methods = Base.invokelatest(Guards.package_methods, planted)
        planted_names = Guards.TOML.parsefile(joinpath(planted.directory, "names.toml"))
        @test Base.invokelatest(Guards.naming_rule, planted_methods, planted, planted_names) ==
              Dict{String, Int}(
            "base | function | Earth | length" => 1,
            "context | function | Earth | earth_layer | earth" => 1,
            "file | Earth | formula | sanctioned twice" => 1,
            "file | global | check | rejected" => 1,
            "file | global | check | unused" => 1,
            "file | global | of | preposition" => 1,
            "file | global | of | unused" => 1,
            "file | global | resistivity | one module" => 1,
            "file | global | unused | unused" => 1,
            "shape | field | Earth | mean_free_path | length" => 1,
            "shape | function | Earth | mesh_build | form" => 1,
            "shape | function | Earth | radial_half_extent | length" => 1,
            "shape | keyword | Earth | outer_layer_count | length" => 1,
            "shape | type | Earth | AbstractDisk | form" => 1,
            "shape | type | Earth | BuildPlan | form" => 1,
            "subject | function | Earth | ring_poses | ring" => 1,
            "terms | function | Earth | abc | ambiguous" => 1,
            "terms | function | Earth | frobnicate | unknown" => 1,
            "terms | function | Earth | get_radius | rejected get -> the name of the quantity" => 1,
            "terms | function | Earth | get_radius | unknown" => 1,
            "verbs | function | Earth | compute | compute" => 1)
        @test Base.invokelatest(Guards.shared_names, planted) ==
              Dict{String, Int}("Plan" => 2, "compute" => 2)

        clean, _ = probe("names/clean", "NamesProbeClean")
        clean_methods = Base.invokelatest(Guards.package_methods, clean)
        clean_names = Guards.TOML.parsefile(joinpath(clean.directory, "names.toml"))
        @test Base.invokelatest(Guards.naming_rule, clean_methods, clean, clean_names) ==
              Dict{String, Int}()
        @test Base.invokelatest(Guards.shared_names, clean) == Dict{String, Int}()
    end
end
