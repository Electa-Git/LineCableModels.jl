# Taxonomy guards. Each guard takes the item inventory as an argument; the controls
# apply it to probe items. The taxonomy is defined in `test/support/taxonomy.jl`.
@testmodule TaxonomyGuards begin
    using TestItemRunner
    include(joinpath(@__DIR__, "..", "support", "runner.jl"))
    const T = ValidationTestRunner
    const REPOSITORY = dirname(dirname(@__DIR__))
    const KNOWN = union(Set(T.OWNERS), Set(T.KINDS), T.ORDINARY_EXCLUDED_TAGS, Set((:slow,)))

    label(item) = string(item.file, " | ", item.name)

    # T1. An item outside `quality` carries exactly one owner tag and a `quality` item
    # none. Every item carries a kind tag, and every tag is known.
    function t1(items)
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

    # T2. An item under `test/unit/<dir>/` carries the owner of `src/<dir>` or a later one.
    function t2(items, owners)
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
    # Module owners follow load order and name every core owner.
    ranks = map(T.rank, values(T.MODULE_OWNERS))
    @test all(>(0), ranks) && issorted(ranks)
    @test Set(values(T.MODULE_OWNERS)) == Set(T.OWNERS[1:T.rank(:pscad)])
    @test Set(values(T.EXTENSION_OWNERS)) == Set(T.OWNERS[T.rank(:pscad)+1:end])
    @test isempty(intersect(Set(T.OWNERS), union(Set(T.KINDS), T.ORDINARY_EXCLUDED_TAGS)))
    # Every Julia file under `src/` and `ext/` has an owner.
    owners = T.loaded_owners(TaxonomyGuards.REPOSITORY)
    for base in ("src", "ext"), (directory, _, names) in walkdir(joinpath(TaxonomyGuards.REPOSITORY, base))
        for name in names
            path = relpath(joinpath(directory, name), TaxonomyGuards.REPOSITORY)
            @test T.path_owner(owners, path) in T.OWNERS
        end
    end
end

@testitem "Quality / taxonomy / T1 owner and kind tags" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    @test G.report(G.t1(G.T.inventory(G.REPOSITORY).items)) == String[]
end

@testitem "Quality / taxonomy / T2 unit directories bound the owner" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    owners = G.T.loaded_owners(G.REPOSITORY)
    @test G.report(G.t2(G.T.inventory(G.REPOSITORY).items, owners)) == String[]
end

@testitem "Quality / taxonomy / negative controls" tags=[:quality] setup=[TaxonomyGuards] begin
    G = TaxonomyGuards
    T = G.T
    item(file, name, tags...) = (; file, name, tags = collect(Symbol, tags), setups = Symbol[])

    # T1 reports each planted defect once and accepts the valid items.
    valid = [item("test/unit/engine/a.jl", "valid", :unit, :engine),
        item("test/quality/q.jl", "quality", :quality),
        item("test/extensions/x.jl", "slow extension", :extension, :fem_numerical, :fem, :slow)]
    @test G.t1(valid) == String[]
    planted = [item("test/unit/engine/a.jl", "no owner", :unit),
        item("test/unit/engine/a.jl", "two owners", :unit, :engine, :uq),
        item("test/quality/q.jl", "quality owner", :quality, :engine),
        item("test/unit/engine/a.jl", "no kind", :engine),
        item("test/unit/engine/a.jl", "unknown", :unit, :engine, :slwo)]
    found = G.t1(planted)
    @test length(found) == 5
    for (probe, text) in zip(planted, ("0 owner tags", "2 owner tags", "1 owner tags",
        "no kind tag", "unknown tags slwo"))
        @test count(f -> startswith(f, G.label(probe)) && occursin(text, f), found) == 1
    end

    # T2 bounds owners under `test/unit/<dir>/` by the earliest owner loaded from `src/<dir>`.
    owners = Dict("src/units/Units.jl" => :units, "src/engine/Engine.jl" => :engine,
        "src/engine/late.jl" => :uq)
    accepted = [item("test/unit/engine/a.jl", "same", :unit, :engine),
        item("test/unit/engine/a.jl", "later", :unit, :uq),
        item("test/unit/core/a.jl", "no source directory", :unit, :units),
        item("test/integration/a.jl", "not a unit directory", :integration, :units),
        item("test/unit/a.jl", "top level", :unit, :units)]
    @test G.t2(accepted, owners) == String[]
    early = item("test/unit/engine/a.jl", "earlier", :unit, :commons)
    @test G.t2([early], owners) ==
        ["$(G.label(early)): commons is earlier than engine, the owner of src/engine"]

    # Files are owned at their load position: a root file takes the owner of the last
    # module included before it, a module entry behind a docstring counts, and a file
    # the package does not load takes the earliest owner of its directory.
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
end
