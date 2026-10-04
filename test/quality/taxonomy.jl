# Taxonomy guards. Each guard takes the item inventory as an argument, and the controls
# apply it to probe items. `test/support/taxonomy.jl` defines the taxonomy.
@testmodule TaxonomyGuards begin
    using TestItemRunner
    include(joinpath(@__DIR__, "..", "support", "runner.jl"))
    const T = ValidationTestRunner
    const REPOSITORY = dirname(dirname(@__DIR__))
    const KNOWN = union(Set(T.OWNERS), Set(T.KINDS), T.ORDINARY_EXCLUDED_TAGS, Set((:slow,)))

    label(item) = string(item.file, " | ", item.name)

    # Owner and kind tags. An item outside `quality` has exactly one owner tag, and a
    # `quality` item has none. Each item has a kind tag and only known tags.
    function owner_and_kind_tags(items)
        found = String[]
        for item in items
            owners = count(in(T.OWNERS), item.tags)
            expected = :quality in item.tags ? 0 : 1
            owners == expected ||
                push!(found, "$(label(item)): $owners owner tags, expected $expected")
            any(in(T.KINDS), item.tags) || push!(found, "$(label(item)): no kind tag")
            unknown = setdiff(item.tags, KNOWN)
            isempty(unknown) ||
                push!(found, "$(label(item)): unknown tags $(join(unknown, ", "))")
        end
        return found
    end

    # Unit directories bound the owner. An item under `test/unit/<dir>/` has the owner of
    # `src/<dir>` or a later one.
    function unit_directory_owners(items, owners)
        found = String[]
        for item in items
            parts = split(item.file, '/')
            length(parts) > 3 && parts[1] == "test" && parts[2] == "unit" || continue
            floor = T.directory_owner(owners, "src/" * parts[3])
            owner = T.owner_tag(item.tags)
            (floor === nothing || owner === nothing) && continue
            T.rank(owner) < T.rank(floor) && push!(found,
                "$(label(item)): $owner is earlier than $floor, the owner of src/$(parts[3])")
        end
        return found
    end

    report(found) = (foreach(println, found); found)
end

@testitem "Quality / taxonomy / single definition" tags=[:quality] setup=[TaxonomyGuards] begin
    T = TaxonomyGuards.T
    @test allunique(T.OWNERS) && allunique(T.KINDS)
    # Module owners follow load order, and each core owner covers at least one module.
    ranks = map(T.rank, values(T.MODULE_OWNERS))
    @test all(>(0), ranks) && issorted(ranks)
    @test Set(values(T.MODULE_OWNERS)) == Set(T.OWNERS[1:T.rank(:pscad)])
    @test Set(values(T.EXTENSION_OWNERS)) == Set(T.OWNERS[T.rank(:pscad)+1:end])
    @test isempty(intersect(Set(T.OWNERS), union(Set(T.KINDS), T.ORDINARY_EXCLUDED_TAGS)))
    # Every Julia file under `src/` and `ext/` has an owner.
    owners = T.loaded_owners(TaxonomyGuards.REPOSITORY)
    for base in ("src", "ext"), file in T.repository_files(joinpath(TaxonomyGuards.REPOSITORY, base))
        path = relpath(file, TaxonomyGuards.REPOSITORY)
        @test T.path_owner(owners, path) in T.OWNERS
    end
end

@testitem "Quality / taxonomy / owner and kind tags" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    @test G.report(G.owner_and_kind_tags(G.T.inventory(G.REPOSITORY).items)) == String[]
end

@testitem "Quality / taxonomy / unit directories bound the owner" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    owners = G.T.loaded_owners(G.REPOSITORY)
    @test G.report(G.unit_directory_owners(G.T.inventory(G.REPOSITORY).items, owners)) == String[]
end

@testitem "Quality / taxonomy / negative controls" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    T = G.T
    item(file, name, tags...) = (; file, name, tags = collect(Symbol, tags), setups = Symbol[])

    # The owner and kind check reports each planted defect once and accepts the valid items.
    valid = [item("test/unit/engine/a.jl", "valid", :unit, :engine),
        item("test/quality/q.jl", "quality", :quality),
        item("test/extensions/x.jl", "slow extension", :extension, :fem_numerical, :fem, :slow)]
    @test G.owner_and_kind_tags(valid) == String[]
    planted = [item("test/unit/engine/a.jl", "no owner", :unit),
        item("test/unit/engine/a.jl", "two owners", :unit, :engine, :uq),
        item("test/quality/q.jl", "quality owner", :quality, :engine),
        item("test/unit/engine/a.jl", "no kind", :engine),
        item("test/unit/engine/a.jl", "unknown", :unit, :engine, :slwo)]
    found = G.owner_and_kind_tags(planted)
    @test length(found) == 5
    for (probe, text) in zip(planted, ("0 owner tags", "2 owner tags", "1 owner tags",
        "no kind tag", "unknown tags slwo"))
        @test count(f -> startswith(f, G.label(probe)) && occursin(text, f), found) == 1
    end

    # Owners under `test/unit/<dir>/` are bounded by the earliest owner loaded from `src/<dir>`.
    owners = Dict("src/units/Units.jl" => :units, "src/engine/Engine.jl" => :engine,
        "src/engine/late.jl" => :uq)
    accepted = [item("test/unit/engine/a.jl", "same", :unit, :engine),
        item("test/unit/engine/a.jl", "later", :unit, :uq),
        item("test/unit/core/a.jl", "no source directory", :unit, :units),
        item("test/integration/a.jl", "not a unit directory", :integration, :units),
        item("test/unit/a.jl", "top level", :unit, :units)]
    @test G.unit_directory_owners(accepted, owners) == String[]
    early = item("test/unit/engine/a.jl", "earlier", :unit, :commons)
    @test G.unit_directory_owners([early], owners) ==
        ["$(G.label(early)): commons is earlier than engine, the owner of src/engine"]

    # Files are owned at their load position. A root file takes the owner of the last
    # module included before it. A module entry behind a docstring counts as a module. A
    # file that the package does not load takes the earliest owner of its directory.
    mktempdir() do repository
        files = Dict(
            "src/LineCableModels.jl" => """
                module LineCableModels
                include("first.jl")
                include("units/Units.jl")
                include("between.jl")
                include("engine/Engine.jl")
                include("modalanalysis/late.jl")
                end
                """,
            "src/first.jl" => "", "src/between.jl" => "", "src/modalanalysis/late.jl" => "",
            "src/units/Units.jl" => "module Units\ninclude(\"inner.jl\")\nend\n",
            "src/units/inner.jl" => "",
            "src/engine/Engine.jl" => "\"Docstring.\"\nmodule Engine\nend\n",
            "ext/LineCableModelsXLSXExt.jl" => "module LineCableModelsXLSXExt\nend\n",
            "ext/LineCableModelsGmshExt/LineCableModelsGmshExt.jl" =>
                "module LineCableModelsGmshExt\ninclude(\"mesh.jl\")\nend\n",
            "ext/LineCableModelsGmshExt/mesh.jl" => "")
        for (path, text) in files
            mkpath(dirname(joinpath(repository, path)))
            write(joinpath(repository, path), text)
        end
        loaded = T.loaded_owners(repository)
        @test loaded == Dict("src/LineCableModels.jl" => :units, "src/first.jl" => :units,
            "src/units/Units.jl" => :units, "src/units/inner.jl" => :units,
            "src/between.jl" => :units, "src/engine/Engine.jl" => :engine,
            "src/modalanalysis/late.jl" => :engine, "ext/LineCableModelsXLSXExt.jl" => :xlsx,
            "ext/LineCableModelsGmshExt/LineCableModelsGmshExt.jl" => :fem,
            "ext/LineCableModelsGmshExt/mesh.jl" => :fem)
        @test T.path_owner(loaded, "src/engine/data.toml") === :engine
        @test T.path_owner(loaded, "ext/LineCableModelsGmshExt/remote/run.py") === :fem
        @test T.path_owner(loaded, "src/new/file.jl") === :units
        @test T.path_owner(loaded, "test/support/runner.jl") === nothing
        @test T.directory_owner(loaded, "src/modalanalysis") === :engine
        write(joinpath(repository, "src/twice.jl"), "")
        write(joinpath(repository, "src/LineCableModels.jl"),
            "module LineCableModels\ninclude(\"twice.jl\")\ninclude(\"twice.jl\")\nend\n")
        @test_throws ErrorException T.loaded_owners(repository)
    end

    # The guards read the tracked files and the new files that git would add, never the
    # ignored ones. Outside a git work tree, they read every file.
    mktempdir() do repository
        git(arguments...) = run(pipeline(Cmd(["git", "-C", repository, arguments...]);
            stdout = devnull, stderr = devnull))
        for (path, text) in ("test/tracked.jl" => "", "test/new.jl" => "",
            "test/captured/old.jl" => "", ".gitignore" => "captured/\n")
            mkpath(dirname(joinpath(repository, path)))
            write(joinpath(repository, path), text)
        end
        directory = joinpath(repository, "test")
        @test T.repository_files(directory) == joinpath.(directory,
            ["captured/old.jl", "new.jl", "tracked.jl"])
        git("init", "-q")
        git("add", "test/tracked.jl", ".gitignore")
        @test T.repository_files(directory) == joinpath.(directory, ["new.jl", "tracked.jl"])
        rm(joinpath(directory, "tracked.jl"))
        @test T.repository_files(directory) == [joinpath(directory, "new.jl")]
    end
end

@testitem "Quality / taxonomy / changed: selection controls" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    T = G.T
    (; items, setups) = T.inventory(G.REPOSITORY)
    owners = T.loaded_owners(G.REPOSITORY)
    select(paths...) = T.changed_items(collect(String, paths), items, setups, owners)
    key(item) = (item.file, item.name)
    environment(item) = any(in(T.ENVIRONMENT_TAGS), item.tags)
    quality = Set(key(i) for i in items if :quality in i.tags)
    owner(item) = T.owner_tag(item.tags)

    # A clean tree, or a change outside `src/`, `ext/` and `test/`, selects only quality.
    @test select() == quality
    @test select("docs/src/developers.md", "Project.toml") == quality

    # A change in `src/engine/` selects the non-slow items owned by engine or later.
    engine = select("src/engine/lineparameters.jl")
    expected = Set(key(i) for i in items if !environment(i) && :slow ∉ i.tags &&
        owner(i) !== nothing && T.rank(owner(i)) >= T.rank(:engine))
    @test engine == union(quality, expected)
    @test any(i -> key(i) in engine && owner(i) === :engine, items)
    @test any(i -> key(i) in engine && owner(i) === :makie, items)
    @test !any(i -> key(i) in engine && owner(i) === :datamodel, items)
    @test !any(i -> key(i) in engine && :slow in i.tags, items)
    # Root files and extensions count at their load position.
    @test select("src/performance.jl") == select("src/uq/UQ.jl")
    @test select("src/modalanalysis/delegation.jl") == select("src/parametricbuilder/ParametricBuilder.jl")
    @test select("src/engine/lineparameters.jl", "src/uq/UQ.jl") == engine
    fem = select("ext/LineCableModelsGmshExt/LineCableModelsGmshExt.jl")
    @test all(k -> k in quality || owner(only(i for i in items if key(i) == k)) in (:fem, :makie), fem)

    # A changed test file selects all its items, slow ones too, and the files that use
    # a setup it defines. Its environment items remain excluded.
    file = "test/unit/engine/unified_earth_return.jl"
    @test any(i -> i.file == file && :slow in i.tags, items)
    @test select(file) == union(quality, Set(key(i) for i in items if i.file == file))
    tracking = select("test/unit/modalanalysis/vieira2026.jl")
    @test any(k -> first(k) == "test/unit/modalanalysis/wedepohl1996.jl", tracking)
    @test select("test/extensions/fem_artifact.jl") == quality

    # A change under `test/support/` selects every non-slow item and every item under
    # `test/unit/core/`, so a planted change to the runner selects the harness test.
    support = select("test/support/runner.jl")
    harness = [i for i in items if i.file == "test/unit/core/test_harness.jl"]
    @test length(harness) == 2 && all(i -> :slow in i.tags && key(i) in support, harness)
    @test support == union(quality, Set(key(i) for i in items if !environment(i) &&
        (:slow ∉ i.tags || startswith(i.file, "test/unit/core/"))))

    # The changed paths come from git: tracked changes, deletions and untracked files.
    mktempdir() do repository
        git(arguments...) = run(pipeline(Cmd(["git", "-C", repository, arguments...]);
            stdout = devnull, stderr = devnull))
        git("init", "--quiet")
        git("config", "user.email", "probe@example.invalid")
        git("config", "user.name", "probe")
        mkpath(joinpath(repository, "src", "engine"))
        foreach(f -> write(joinpath(repository, f), "x"), ("src/engine/a.jl", "src/b.jl", "c.md"))
        git("add", "--all")
        git("commit", "--quiet", "-m", "probe")
        @test T.changed_paths("HEAD", repository) == String[]
        write(joinpath(repository, "src", "engine", "a.jl"), "y")
        rm(joinpath(repository, "src", "b.jl"))
        write(joinpath(repository, "d.jl"), "z")
        git("mv", "c.md", "e.md")
        @test T.changed_paths("HEAD", repository) == ["c.md", "d.jl", "e.md", "src/b.jl", "src/engine/a.jl"]
        @test_throws ErrorException T.changed_paths("no-such-revision", repository)
    end
    # `changed:REF` takes no other selector and runs only through `execute`.
    @test T.changed_reference(["changed:HEAD", "--list"]) == "HEAD"
    @test T.changed_reference(["tag:unit", "--list"]) === nothing
    @test_throws ErrorException T.changed_reference(["changed:HEAD", "tag:unit"])
    @test_throws ErrorException T.changed_reference(["changed:"])
    @test_throws ErrorException T.selection(["changed:HEAD"], joinpath(G.REPOSITORY, "test"))
end
