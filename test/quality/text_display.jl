@testitem "Quality / TextDisplay / ownership and side-effect boundaries" tags=[:quality] begin
    using DataFrames

    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    display_paths=String[
        joinpath(source_root, "textdisplay", "TextDisplay.jl"),
        joinpath(source_root, "materials", "base.jl"),
        joinpath(source_root, "earthprops", "base.jl")
    ]
    append!(display_paths,
        [joinpath(source_root, owner, "textdisplay.jl")
         for owner in (
            "datamodel",
            "engine",
            "parametricbuilder",
            "reportbuilder",
            "uq",
            "validation"
        )])
    display_source=join(read(path, String) for path in display_paths)

    for forbidden in (
        r"\bbuild\s*\(",
        r"\bresolve\s*\(",
        r"\bcompute\s*\(",
        r"\bobservables\s*\(",
        r"\bDataFrame\s*\(",
        r"\bmaterialize\s*\(",
        r"\btessellate\s*\(",
        r"\b(?:CairoMakie|GLMakie|Makie)\."
    )
        @test !occursin(forbidden, display_source)
    end
    @test !occursin("fieldnames", display_source)
    @test !occursin("PrettyTables", display_source)
    @test length(collect(eachmatch(r"\bmacro\s+showfields\b", display_source))) == 1
    @test !Base.ispublic(LineCableModels.TextDisplay, Symbol("@showfields"))

    owner_modules=(
        LineCableModels.Units,
        LineCableModels.Validation,
        LineCableModels.Materials,
        LineCableModels.EarthProps,
        LineCableModels.DataModel,
        LineCableModels.Engine,
        LineCableModels.ParametricBuilder,
        LineCableModels.UQ,
        LineCableModels.ReportBuilder
    )
    for binding in names(LineCableModels; all = false, imported = true)
        isdefined(LineCableModels, binding) || continue
        owned_type=getfield(LineCableModels, binding)
        owned_type isa Union{DataType, UnionAll} || continue
        Base.isabstracttype(owned_type) && continue
        parentmodule(Base.unwrap_unionall(owned_type)) in owner_modules || continue
        @test which(summary, (IO, owned_type)).module !== Base
        @test which(show, (IO, owned_type)).module !== Base
        @test which(show, (IO, MIME"text/plain", owned_type)).module !== Base
    end
    for owner in owner_modules
        for binding in names(owner; all = false, imported = false)
            isdefined(owner, binding) || continue
            owned_type=getfield(owner, binding)
            owned_type isa Union{DataType, UnionAll} || continue
            Base.isabstracttype(owned_type) && continue
            parentmodule(Base.unwrap_unionall(owned_type)) in owner_modules || continue
            @test which(summary, (IO, owned_type)).module !== Base
            @test which(show, (IO, owned_type)).module !== Base
            @test which(show, (IO, MIME"text/plain", owned_type)).module !== Base
        end
    end

    report_methods=filter(method -> method.module === LineCableModels.ReportBuilder,
        collect(methods(DataFrame)))
    @test length(report_methods) == 1
    signature=Base.unwrap_unionall(only(report_methods).sig)
    @test signature.parameters[2] <: LineCableModels.Grammar.ObservationPublication

    tables_owners=String[]
    for (directory, _, files) in walkdir(source_root)
        for file in files
            endswith(file, ".jl") || continue
            path=joinpath(directory, file)
            occursin(r"Tables\.(?:istable|rowaccess|columnaccess|columns|schema)\s*\(",
                read(path, String)) && push!(tables_owners, relpath(path, root))
        end
    end
    @test tables_owners == [joinpath("src", "grammar", "observables.jl")]

    for owner in ("datamodel", "engine")
        owner_source=join(
            read(joinpath(directory, file), String)
        for (directory, _, files) in walkdir(joinpath(source_root, owner))
        for file in files if endswith(file, ".jl")
        )
        @test !occursin("DataFrames", owner_source)
    end

    tutorial_source=join(
        read(joinpath(root, "examples", tutorial), String)
    for tutorial in ("tutorial2.jl", "tutorial3.jl")
    )
    for object in (
        "materials",
        "constants",
        "equivalent_design",
        "cable_design",
        "library",
        "earth_params",
        "cable_system",
        "line_parameters",
        "sequence_parameters"
    )
        @test !occursin(Regex("DataFrame\\s*\\(\\s*" * object * "\\s*[),]"),
            tutorial_source)
    end
end
