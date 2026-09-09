function write_source_toml(path, data)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
end

function source_fixture(directory)
    root = joinpath(directory, "profile")
    dep = joinpath(directory, "dependency")
    root_id = "b7f013d4-6e73-4b8e-ab87-30a5f4a17f8b"
    dep_id = "21c0ee66-688f-4126-a3d8-a40b79bcedc5"
    for (path, name, id) in ((root, "FixtureProfile", root_id), (dep, "FixtureDependency", dep_id))
        mkpath(joinpath(path, "src"))
        write_source_toml(joinpath(path, "Project.toml"), Dict("name"=>name, "uuid"=>id, "version"=>"0.1.0"))
        # Loading either package would fail. Source inspection must not evaluate it.
        write(joinpath(path, "src", name * ".jl"), "error(\"source inspection evaluated code\")\n")
    end
    project = TOML.parsefile(joinpath(root, "Project.toml"))
    project["deps"] = Dict("FixtureDependency"=>dep_id)
    write_source_toml(joinpath(root, "Project.toml"), project)
    manifest = Dict("manifest_format"=>"2.0", "julia_version"=>string(VERSION),
        "deps"=>Dict("FixtureDependency"=>[Dict("uuid"=>dep_id, "path"=>"../dependency", "version"=>"0.1.0")]))
    write_source_toml(joinpath(root, "Manifest.toml"), manifest)
    return root, dep
end

@testset "native source fingerprint is passive and complete for declared inputs" begin
    mktempdir() do directory
        root, dep = source_fixture(directory)
        first = native_environment_fingerprint(root)
        @test first.package == "FixtureProfile"
        @test first.uuid == UUID("b7f013d4-6e73-4b8e-ab87-30a5f4a17f8b")
        @test first.files == 5
        @test first.bytes > 0
        @test length(first.digest) == 64
        @test native_environment_fingerprint(root).digest == first.digest
        @test !any(id.name in ("FixtureProfile", "FixtureDependency") for id in keys(Base.loaded_modules))
        mktempdir() do copy_root
            copied, _ = source_fixture(copy_root)
            @test native_environment_fingerprint(copied).digest == first.digest
        end
        write(joinpath(root, "README.md"), "Non-executed documentation")
        @test native_environment_fingerprint(root).digest == first.digest
        source = joinpath(dep, "src", "FixtureDependency.jl")
        write(source, read(source, String) * "# changed local dependency\n")
        changed = native_environment_fingerprint(root)
        @test changed.digest != first.digest
        profile = ProfileDefinition("fixture", root, changed.digest; operations=("fixture.operation",))
        @test verify_native_environment(profile).digest == changed.digest
        @test_throws ArgumentError verify_native_environment(ProfileDefinition("fixture", root, first.digest; operations=("fixture.operation",)))
        @test_throws ArgumentError verify_native_environment(ProfileDefinition("terminal",
            "registry.test/julia@sha256:" * changed.digest, changed.digest; kind=:terminal, isolation=:container))
        previous = changed.digest
        for (name, bytes) in (("LocalPreferences.toml", "[FixtureDependency]\nvalue = 2\n"),
                ("Artifacts.toml", "[fixture]\ngit-tree-sha1 = '" * repeat("a", 40) * "'\n"),
                ("ext/FixtureExtension.jl", "# optional extension\n"),
                ("deps/build.jl", "# build input\n"))
            path = joinpath(dep, name)
            mkpath(dirname(path))
            write(path, bytes)
            next = native_environment_fingerprint(root).digest
            @test next != previous
            previous = next
        end
        shared = joinpath(directory, "shared.jl")
        write(shared, "# shared executor implementation\n")
        write_source_toml(joinpath(root, "RuntimeSources.toml"), Dict("schema_version"=>1, "files"=>["../shared.jl"]))
        next = native_environment_fingerprint(root).digest
        @test next != previous
        write(shared, "# changed shared executor implementation\n")
        @test native_environment_fingerprint(root).digest != next
        @test_throws ArgumentError native_environment_fingerprint(root; max_files=true)
        @test_throws ArgumentError native_environment_fingerprint(root; max_files=2)
        @test_throws ArgumentError native_environment_fingerprint(root; max_file_bytes=16)
        @test_throws ArgumentError native_environment_fingerprint(root; max_bytes=100)
        @test_throws ArgumentError native_environment_fingerprint(root; max_bytes=typemax(Int))
        output = joinpath(directory, "fingerprint.txt")
        open(output, "w") do io
            redirect_stdout(io) do
                @test runtime_cli(["runtime", "fingerprint", "--project", root]) === nothing
            end
        end
        @test strip(read(output, String)) == native_environment_fingerprint(root).digest
        @test_throws ArgumentError runtime_cli(["runtime", "fingerprint", "--config", "unused"])
        @test_throws ArgumentError runtime_cli(["runtime", "fingerprint", "--project", root, "--xray"])
    end
end

@testset "native source inspection fails closed" begin
    mktempdir() do directory
        root, dep = source_fixture(directory)
        path = joinpath(root, "Manifest.toml")
        good = read(path, String)
        manifest = TOML.parse(good)
        manifest["julia_version"] = "0.0.0"
        write_source_toml(path, manifest)
        @test_throws ArgumentError native_environment_fingerprint(root)
        write(path, good)
        versioned = joinpath(root, "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
        write(versioned, good)
        chosen = native_environment_fingerprint(root).digest
        write(path, "unselected manifest cannot parse")
        @test native_environment_fingerprint(root).digest == chosen
        rm(versioned)
        @test_throws ArgumentError native_environment_fingerprint(root)
        write(path, good)
        manifest = TOML.parse(good)
        manifest["deps"]["Unpinned"] = [Dict("uuid"=>string(uuid4()), "version"=>"1.0.0")]
        write_source_toml(path, manifest)
        @test_throws ArgumentError native_environment_fingerprint(root)
        manifest["deps"]["Unpinned"][1]["git-tree-sha1"] = repeat("a", 40)
        write_source_toml(path, manifest)
        @test native_environment_fingerprint(root).files == 5
        write(path, good)
        dependency_project = joinpath(dep, "Project.toml")
        dep_good = read(dependency_project, String)
        bad_version = TOML.parse(dep_good)
        bad_version["version"] = "9.0.0"
        write_source_toml(dependency_project, bad_version)
        @test_throws ArgumentError native_environment_fingerprint(root)
        write(dependency_project, dep_good)
        spoofed_stdlib = TOML.parse(good)
        spoofed_stdlib["deps"]["ForeignDates"] = [Dict("uuid"=>"ade2ca70-3891-5945-98fb-dc099432e06a")]
        write_source_toml(path, spoofed_stdlib)
        @test_throws ArgumentError native_environment_fingerprint(root)
        write(path, good)
        write_source_toml(dependency_project, Dict("name"=>"Foreign", "uuid"=>string(uuid4())))
        @test_throws ArgumentError native_environment_fingerprint(root)
        write(dependency_project, dep_good)
        extra = joinpath(root, "RuntimeSources.toml")
        for invalid in (Dict("schema_version"=>true, "files"=>String[]),
                Dict("schema_version"=>1, "files"=>["../dependency"]),
                Dict("schema_version"=>1, "files"=>["/etc/passwd"]),
                Dict("schema_version"=>1, "files"=>["../*.jl"]),
                Dict("schema_version"=>1, "files"=>["Project.toml", "./Project.toml"]),
                Dict("schema_version"=>1, "files"=>String[], "typo"=>1))
            write_source_toml(extra, invalid)
            @test_throws ArgumentError native_environment_fingerprint(root)
        end
        rm(extra)
        link = joinpath(root, "src", "linked.jl")
        symlink(joinpath(dep, "src", "FixtureDependency.jl"), link)
        @test_throws ArgumentError native_environment_fingerprint(root)
        rm(link)
        symlink(dep, joinpath(directory, "linked-package"))
        @test_throws ArgumentError native_environment_fingerprint(joinpath(directory, "linked-package"))
        write(joinpath(root, "JuliaProject.toml"), dep_good)
        @test_throws ArgumentError native_environment_fingerprint(root)
        rm(joinpath(root, "JuliaProject.toml"))
        write(joinpath(root, "JuliaManifest.toml"), good)
        @test_throws ArgumentError native_environment_fingerprint(root)
        rm(joinpath(root, "JuliaManifest.toml"))
        inventory = RT.SourceInventory(Dict{String,String}(), Dict{String,Tuple}(), 0, 0, 100, 10000, 100000)
        RT.fingerprint_file!(inventory, path, "manifest")
        write(path, good * "# concurrent operator edit\n")
        @test_throws ArgumentError RT.verify_source_stamps(inventory)
        @test_throws ArgumentError RT.remember_source!(inventory, path)
    end
end

@testset "shipped numerical source closures" begin
    worker = normpath(joinpath(@__DIR__, "..", "..", "worker"))
    line = native_environment_fingerprint(joinpath(worker, "profiles", "line-parameters"))
    flow = native_environment_fingerprint(joinpath(worker, "profiles", "power-flow"))
    @test line.package == "LineCableModelsLineParameters"
    @test flow.package == "LineCableModelsPowerFlow"
    @test line.digest != flow.digest
    @test line.files > flow.files > 10
    @test !any(id.name in ("LineCableModels", "PowerImpedance", "Bonito") for id in keys(Base.loaded_modules))
    for (package, files) in (("core", ["../src/OperationRegistry.jl", "../src/Cache.jl", "../src/Executor.jl"]),
            ("profiles/line-parameters", ["../../src/operations/linecablemodels.jl"]),
            ("profiles/power-flow", ["../../src/operations/powerimpedance.jl"]))
        @test TOML.parsefile(joinpath(worker, package, "RuntimeSources.toml"))["files"] == files
        package_root = joinpath(worker, package)
        declared = Set(abspath(joinpath(package_root, file)) for file in files)
        # Derive literal cross-package includes from owned code. An added shared
        # include must be declared too, not silently omitted from its fingerprint.
        for (directory, _, names) in walkdir(joinpath(package_root, "src")), name in names
            endswith(name, ".jl") || continue
            for match in eachmatch(r"include\(\"([^\"]+)\"\)", read(joinpath(directory, name), String))
                target = abspath(joinpath(directory, match.captures[1]))
                startswith(target, abspath(joinpath(package_root, "src")) * "/") && continue
                @test target in declared
            end
        end
    end
end
