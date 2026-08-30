@testitem "Core / architecture / removed infrastructure stays removed" tags=[:unit] begin
    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    source_files=filter(endswith(".jl"),
        [joinpath(directory, file)
         for (directory, _, files) in walkdir(source_root)
         for file in files])
    source=Dict(relpath(path, root)=>read(path, String) for path in source_files)

    @test !isdefined(LineCableModels, :Commons)
    @test !isdefined(LineCableModels, :Utils)
    @test !ispath(joinpath(source_root, "commons"))
    @test !ispath(joinpath(source_root, "utils"))
    @test isempty(filter(path -> basename(path) == "helpers.jl", source_files))

    forbidden=(
        "resolve_T", "coerce_to_T", "BASE_FLOAT", "REALSCALAR", "COMPLEXSCALAR",
        "@construct", "@parameterize", "@measurify", "Meta.parse", "global_logger(",
        "line_component_quantity", "line_component_unit", "parse_kwargs",
        "resolve_input", "axis_quantity", "axis_unit", "series_data",
        "ComponentMetadata", "UnitSpec", "DEFAULT_QUANTITY_UNITS",
        "_mc_selector", "AbstractObservable", "ObservableDescriptor",
        "ObservationRequest", "AbstractDetails", "DetailsRegistry",
        "_valdispatch", "symbol_dispatch", "@symbol_facade"
    )
    for token in forbidden
        @test all(!occursin(token, contents) for contents in values(source))
    end
    @test all(!occursin("@assert", contents) for contents in values(source))

    macro_path=joinpath("src", "parametricbuilder", "macros.jl")
    macro_facades=Set((
        joinpath("src", "LineCableModels.jl"),
        joinpath("src", "parametricbuilder", "ParametricBuilder.jl")
    ))
    for (path, contents) in source
        path == macro_path && continue
        if path in macro_facades
            @test occursin("@relax", contents)
            @test !occursin(r"\brecast\b", contents)
            continue
        end
        @test !occursin("@relax", contents)
        @test !occursin(r"\brecast\b", contents)
    end

    @test !isfile(joinpath(source_root, "retired.jl"))
    @test all(!occursin("calc_", contents) for contents in values(source))

    construction_files=(
        joinpath("src", "datamodel", "design", "region.jl"),
        joinpath("src", "datamodel", "design", "stack.jl"),
        joinpath("src", "datamodel", "design", "group.jl"),
        joinpath("src", "datamodel", "design", "assembly.jl"),
        joinpath("src", "datamodel", "design", "enclosure.jl"),
        joinpath("src", "datamodel", "design", "cabledesign.jl"),
        joinpath("src", "datamodel", "linecablesystem", "linecablesystem.jl")
    )
    for path in construction_files
        contents=source[path]
        @test !occursin(r"\btemperature\s*::", contents)
        @test !occursin(r";\s*temperature\s*=", contents)
    end
end

@testitem "Core / selectors / owner-local symbol facades" tags=[:unit] begin
    const Engine=LineCableModels.Engine
    const ImportExport=LineCableModels.ImportExport
    const PlotBuilder=LineCableModels.PlotBuilder
    const WirePatterns=LineCableModels.ParametricBuilder.WirePatterns

    @test which(Engine.Formulation, (Symbol,)).module === Engine
    @test which(ImportExport.export_data, (Symbol, Nothing)).module === ImportExport
    @test which(ImportExport.import_data, (Symbol, Nothing)).module === ImportExport
    @test which(Engine.LineParameters, (Symbol, String)).module === ImportExport
    @test !isdefined(PlotBuilder, :layout_preset)
    @test which(
        getindex,
        (WirePatterns.WireEstimate{Float64, WirePatterns.HexaPattern{Float64}}, Symbol)
    ).module === WirePatterns

    @test Engine.Formulation() isa Engine.LineParametersFormulation
    @test_throws MethodError Engine.Formulation(:analytical)
    @test_throws MethodError Engine.Formulation(:line_cable_models)
    @test_throws MethodError ImportExport.export_data(:unregistered, nothing)
    @test_throws MethodError ImportExport.import_data(:unregistered, nothing)
end

@testitem "Core / architecture / ownership-centered recursive layout" tags=[:unit] begin
    using TOML

    root=pkgdir(LineCableModels)
    source_root=joinpath(root, "src")
    extension_root=joinpath(root, "ext")
    source_files=filter(endswith(".jl"),
        [joinpath(directory, file)
         for (directory, _, files) in walkdir(source_root)
         for file in files])
    source=Dict(relpath(path, root)=>read(path, String) for path in source_files)

    expected_entries=(
        joinpath("src", "grammar", "Grammar.jl"),
        joinpath("src", "units", "Units.jl"),
        joinpath("src", "validation", "Validation.jl"),
        joinpath("src", "plotbuilder", "PlotBuilder.jl"),
        joinpath("src", "datamodel", "DataModel.jl"),
        joinpath("src", "engine", "Engine.jl"),
        joinpath("src", "transforms", "Transforms.jl"),
        joinpath("src", "uq", "UQ.jl"),
        joinpath("src", "reportbuilder", "ReportBuilder.jl"),
        joinpath("src", "importexport", "ImportExport.jl"),
        joinpath("ext", "LineCableModelsMakieExt", "LineCableModelsMakieExt.jl"),
        joinpath("ext", "LineCableModelsXLSXExt.jl")
    )
    @test all(isfile(joinpath(root, path)) for path in expected_entries)

    expected_owned_files=(
        joinpath("src", "datamodel", "geometry", "primitives.jl"),
        joinpath("src", "datamodel", "geometry", "shell.jl"),
        joinpath("src", "datamodel", "geometry", "roundedsector.jl"),
        joinpath("src", "datamodel", "geometry", "resolve.jl"),
        joinpath("src", "datamodel", "design", "region.jl"),
        joinpath("src", "datamodel", "design", "cabledesign.jl"),
        joinpath("src", "datamodel", "placement", "patterns.jl"),
        joinpath("src", "datamodel", "placement", "paths.jl"),
        joinpath("src", "parametricbuilder", "physicaltree.jl"),
        joinpath("src", "datamodel", "preview", "cable.jl"),
        joinpath("src", "datamodel", "preview", "cables.jl"),
        joinpath("src", "engine", "lineparameters", "lineparameters.jl"),
        joinpath("src", "engine", "lineparameters", "quantities.jl"),
        joinpath("src", "engine", "earthreturn.jl"),
        joinpath("src", "engine", "insulationimpedance", "interface.jl"),
        joinpath("src", "engine", "insulationimpedance", "formulas", "lossless.jl"),
        joinpath("src", "engine", "insulationadmittance", "interface.jl"),
        joinpath("src", "engine", "insulationadmittance", "formulas", "gustavsen2013.jl"),
        joinpath("src", "engine", "insulationadmittance", "formulas", "marti2001.jl"),
        joinpath("src", "engine", "semiconadmittance", "interface.jl"),
        joinpath("src", "engine", "semiconadmittance", "formulas", "ametani2004.jl"),
        joinpath("src", "transforms", "compute.jl"),
        joinpath("src", "transforms", "eigensystems.jl"),
        joinpath("src", "transforms", "quantities.jl"),
        joinpath("src", "transforms", "formulas", "chrysochos2014.jl"),
        joinpath("src", "transforms", "formulas", "fan2009.jl"),
        joinpath("src", "transforms", "formulas", "fortescue.jl"),
        joinpath("src", "transforms", "formulas", "wedepohl1996.jl"),
        joinpath("src", "uq", "montecarlo", "compute.jl"),
        joinpath("src", "reportbuilder", "grammar.jl"),
        joinpath("src", "reportbuilder", "tables.jl"),
        joinpath("src", "reportbuilder", "xlsx.jl"),
        joinpath("src", "importexport", "pscad", "pscad.jl"),
        joinpath("ext", "LineCableModelsMakieExt", "UIComponents.jl"),
        joinpath("dev", "plotting", "Project.toml"),
        joinpath("schemas", "cable-design.schema.json")
    )
    @test all(isfile(joinpath(root, path)) for path in expected_owned_files)

    obsolete_paths=(
        joinpath("src", "Grammar.jl"),
        joinpath("src", "unithandler"),
        joinpath("src", "plotbuilder", "backendhandler"),
        joinpath("src", "plotbuilder", "uicomponents"),
        joinpath("src", "datamodel", "strands_handler.jl"),
        joinpath("src", "engine", "solver.jl"),
        joinpath("src", "engine", "insulationimpedance", "lossless.jl"),
        joinpath("src", "engine", "insulationadmittance", "lossless.jl"),
        joinpath("src", "engine", "insulationadmittance", "parallelrc.jl"),
        joinpath("src", "uq", "plotspecs.jl"),
        joinpath("src", "engine", "lineparameters", "dataframe.jl"),
        joinpath("src", "uq", "dataframe.jl"),
        joinpath("src", "importexport", "xlsx.jl"),
        joinpath("integration", "plotting")
    )
    @test all(!ispath(joinpath(root, path)) for path in obsolete_paths)

    @test Set(readdir(joinpath(source_root, "plotbuilder"))) == Set((
        "PlotBuilder.jl",
        "backends.jl",
        "base.jl",
        "interfaces.jl",
        "render.jl",
        "types.jl"
    ))
    @test Set(readdir(joinpath(extension_root, "LineCableModelsMakieExt"))) == Set((
        "LineCableModelsMakieExt.jl",
        "UIComponents.jl",
        "context.jl",
        "lineparameters.jl",
        "montecarlo.jl",
        "native.jl",
        "previews.jl",
        "shell.jl",
        "uicomponents"
    ))
    @test Set(readdir(joinpath(
        extension_root,
        "LineCableModelsMakieExt",
        "uicomponents"
    ))) == Set((
        "axes.jl",
        "docks.jl",
        "drawing.jl",
        "export.jl",
        "theme.jl"
    ))

    @test !isdefined(LineCableModels, :UnitHandler)
    @test !isdefined(LineCableModels.PlotBuilder, :BackendHandler)
    @test !isdefined(LineCableModels.DataModel, :PhysicalFillLimit)
    @test !isdefined(LineCableModels.Validation, :PhysicalFillLimit)

    @test parentmodule(LineCableModels.description) === LineCableModels
    @test parentmodule(LineCableModels.formula_id) === LineCableModels
    @test parentmodule(LineCableModels.constitutive) === LineCableModels
    @test parentmodule(LineCableModels.Z) === LineCableModels.Engine
    @test parentmodule(LineCableModels.frequencies) === LineCableModels.Engine
    @test parentmodule(LineCableModels.ncables) === LineCableModels.DataModel
    @test parentmodule(LineCableModels.nphases) === LineCableModels.DataModel
    @test parentmodule(LineCableModels.DataModel.preview_shapes) ===
          LineCableModels.DataModel
    @test parentmodule(LineCableModels.DataModel.preview_materials) ===
          LineCableModels.DataModel
    @test Base.ispublic(LineCableModels.DataModel, :preview_shapes)
    @test Base.ispublic(LineCableModels.DataModel, :preview_materials)
    @test parentmodule(LineCableModels.compute) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.observe) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.observables) === LineCableModels.Grammar
    @test parentmodule(LineCableModels.Units.quantity) === LineCableModels.Units
    @test parentmodule(LineCableModels.report) === LineCableModels.ReportBuilder
    @test parentmodule(LineCableModels.ReportArtifact) === LineCableModels.ReportBuilder
    @test parentmodule(LineCableModels.XLSXReportDefinition) ===
          LineCableModels.ReportBuilder
    @test parentmodule(LineCableModels.validate) === LineCableModels.Validation
    @test parentmodule(LineCableModels.build) === LineCableModels
    @test parentmodule(LineCableModels.homogenize) === LineCableModels
    @test :homogenize in names(LineCableModels.DataModel)
    @test :flatten ∉ names(LineCableModels.Engine)
    @test Base.ispublic(LineCableModels.DataModel, :flatten)
    @test LineCableModels.build === LineCableModels.DataModel.build ===
          LineCableModels.ParametricBuilder.build
    @test fieldtype(LineCableModels.Material{Float64}, :kind) === Symbol
    @test !hasfield(LineCableModels.CableDesign, :effective)
    @test !hasfield(LineCableModels.CableDesign, :reference_frequency)
    @test !isdefined(LineCableModels.DataModel, :ResolvedRegion)
    @test !isdefined(LineCableModels.DataModel, :ResolvedPart)
    @test !isdefined(LineCableModels.DataModel, :AbstractPrimitiveDefinition)
    for name in (
        :DiskDefinition, :RectangleDefinition, :EllipseDefinition,
        :SectorDefinition, :AnnulusDefinition, :ShellDefinition, :PolygonDefinition,
        :RoundedSectorDefinition
    )
        @test !isdefined(LineCableModels.DataModel, name)
        @test !isdefined(LineCableModels, name)
    end
    @test !hasmethod(
        LineCableModels.CableDesign,
        Tuple{LineCableModels.AbstractCablePart}
    )

    rejected_cable_model_tokens=(
        "Material{Kind",
        "DesignMaterializer",
        "SystemMaterializer",
        "PartBuilder",
        "PositionDefinition",
        "BuildContext",
        "BuildManager",
        "BuildPlan",
        "BuildSpec",
        "BuildResult",
        "_gridspace_axis",
        "_current_radial_equivalence",
        "compute_cable_constants",
        "CableConstantsProblem",
        "ResolvedRegion",
        "ResolvedPart",
        "AbstractPrimitiveDefinition",
        "DiskDefinition",
        "RectangleDefinition",
        "EllipseDefinition",
        "SectorDefinition",
        "AnnulusDefinition",
        "ShellDefinition",
        "PolygonDefinition",
        "RoundedSectorDefinition",
        "SectorParams",
        "SectorInsulator",
        "_squashed_sector_poly",
        ".effective"
    )
    for token in rejected_cable_model_tokens
        @test all(!occursin(token, contents) for contents in values(source))
    end

    dataframe_sources=filter(keys(source)) do path
        occursin("DataFrame", source[path])
    end
    @test all(startswith(path, joinpath("src", "reportbuilder"))
    for path in dataframe_sources)

    rounded_sector_source=source[joinpath(
        "src", "datamodel", "geometry", "roundedsector.jl"
    )]
    flatten_source=source[joinpath("src", "datamodel", "flatten.jl")]
    for token in ("Point2f", "GeometryBasics", "Makie")
        @test !occursin(token, rounded_sector_source)
        @test !occursin(token, flatten_source)
    end
    @test !occursin("tessellate", flatten_source)
    for token in (
        "LineParametersProblem",
        "LineParametersFormulation",
        "LineCableModelsCoaxial",
        "Matrix{"
    )
        @test !occursin(token, flatten_source)
    end

    ci_workflow=read(joinpath(root, ".github", "workflows", "CI.yml"), String)
    @test occursin(r"(?m)^  quality:\s*$", ci_workflow)
    @test occursin("push!(ARGS, \"tag:quality\")", ci_workflow)
    orchestrator_quality=read(
        joinpath(root, "test", "quality", "orchestrators.jl"),
        String
    )
    @test occursin("@test isempty(findings)", orchestrator_quality)
    @test !occursin("@test_skip", orchestrator_quality)
    @test !occursin("Advisory orchestrator", orchestrator_quality)

    report_grammar=source[joinpath("src", "reportbuilder", "grammar.jl")]
    @test occursin("@required AbstractReportDefinition begin", report_grammar)
    @test occursin(
        "encode(::AbstractReportDefinition, source, published, table, illustration)",
        report_grammar
    )
    @test occursin(
        "write(::AbstractReportDefinition, source, published, table, illustration, encoded)",
        report_grammar
    )
    report_monte_carlo=source[joinpath("src", "reportbuilder", "montecarlo.jl")]
    @test !occursin(
        r"\bsource\.(?:formulation|values|stats|sample_values|histogram_values|root_seed|point_seeds|trial_counts|details)\b",
        report_monte_carlo
    )
    @test !occursin(r"\bUQ\.(?:statistics|samples|histograms)\(source\)",
        report_monte_carlo)
    @test length(findall("observables(", report_monte_carlo)) == 1
    report_grammar=source[joinpath("src", "reportbuilder", "grammar.jl")]
    @test occursin(
        "PlotBuilder.make_render(\n        definition.illustration,\n        published;",
        report_grammar
    )
    report_xlsx=source[joinpath("src", "reportbuilder", "xlsx.jl")]
    @test !occursin(r"\b(?:Z|Y|frequencies)\(source", report_xlsx)
    @test !occursin("_xlsx_table_definition", report_xlsx)
    @test !occursin("_select_line_table", report_xlsx)
    @test !occursin("_tabulate_line_table", report_xlsx)
    @test occursin("observation_columns(table)", report_xlsx)
    report_tables=source[joinpath("src", "reportbuilder", "tables.jl")]
    @test !occursin("zero_atol", report_tables)
    @test !occursin("relative_status", report_tables)
    @test !occursin(r"\bfamily\s*=", report_tables)
    @test !occursin(r"\bquantity\s*=", report_tables)
    @test !occursin(r"\bvalue\s*=", report_tables)
    importexport_index=source[joinpath("src", "importexport", "ImportExport.jl")]
    reportbuilder_index=source[joinpath("src", "reportbuilder", "ReportBuilder.jl")]
    @test !occursin("using XLSX", importexport_index)
    @test !occursin("include(\"xlsx.jl\")", importexport_index)
    @test !occursin("using XLSX", reportbuilder_index)
    @test occursin("include(\"xlsx.jl\")", reportbuilder_index)
    @test all(!occursin(r"^\s*(using|import)\s+XLSX"m, contents)
    for contents in values(source))
    xlsx_extension=read(joinpath(extension_root, "LineCableModelsXLSXExt.jl"), String)
    @test occursin("using XLSX", xlsx_extension)
    @test occursin("function ReportBuilder.write", xlsx_extension)

    plotbuilder_source=join(
        (source[path]
        for path in keys(source)
        if startswith(path, joinpath("src", "plotbuilder"))),
        "\n")
    makie_source=join(
        (read(path, String)
        for (directory, _, files) in walkdir(
            joinpath(extension_root, "LineCableModelsMakieExt"))
        for file in files if endswith(file, ".jl")
        for path in (joinpath(directory, file),)),
        "\n")
    maintained_implementation=join((values(source)..., makie_source), "\n")
    for token in (
        "input_kwargs",
        "renderer_kwargs",
        "_symbol_tuple",
        "_validate_defaults",
        "_BACKEND_EXTENSIONS",
        "_backend_extension",
        "_backend_package"
    )
        @test !occursin(token, plotbuilder_source)
    end
    for token in (
        "hasproperty(page.payload",
        "_page_legend",
        "_page_colorbars",
        "_page_export",
        "dispatch_on",
        "_shell_kind",
        "format_axes!",
        "MCDistributionPlotDefinition",
        "_mc_layers",
        "line_requests",
        "line_parent"
    )
        @test !occursin(token, maintained_implementation)
    end
    @test !occursin(r"\bquantities\s*=", maintained_implementation)
    @test !occursin(r"\bcon\s*=", maintained_implementation)

    preview_materials=source[joinpath("src", "datamodel", "preview", "materials.jl")]
    @test !occursin(r"\blayer\s+isa\b", preview_materials)
    @test !occursin("hasproperty(layer", preview_materials)
    @test occursin("preview_shapes(region::PlacedRegion", preview_materials)
    @test occursin("preview_materials(region::PlacedRegion", preview_materials)

    shell_source=read(
        joinpath(extension_root, "LineCableModelsMakieExt", "shell.jl"),
        String
    )
    @test occursin("function build_widget! end", shell_source)
    @test occursin("function toolbar_button!", shell_source)
    @test occursin("function bind_widget_callback!", shell_source)
    @test occursin("build_shell(context, recipe.definition, page)", shell_source)

    monte_carlo_renderer=read(
        joinpath(extension_root, "LineCableModelsMakieExt", "montecarlo.jl"),
        String
    )
    for verb in ("hist", "stairs", "ecdfplot", "lines", "qqplot")
        @test occursin("function Makie.$verb", monte_carlo_renderer)
    end
    @test !occursin(r"\b(?:mode|ijk)\s*=", monte_carlo_renderer)
    @test !occursin(r"Val\(:(?:histogram|stairs|line|scatter)", monte_carlo_renderer)
    @test !occursin("Figure(", monte_carlo_renderer)
    @test occursin("PlotBuilder.plotwindow", monte_carlo_renderer)
    @test !occursin("LineCableModels.result(", monte_carlo_renderer)
    @test !occursin("LineCableModels.histograms(", monte_carlo_renderer)
    @test !occursin("LineCableModels.UQ._histogram", monte_carlo_renderer)
    for token in (
        "_monte_carlo_request",
        "_distribution_observation",
        "_register_distribution!",
        "_distribution_window"
    )
        @test !occursin(token, monte_carlo_renderer)
    end

    grammar_observables=source[joinpath("src", "grammar", "observables.jl")]
    @test !occursin("applicable(", grammar_observables)
    positions=source[joinpath("src", "parametricbuilder", "positions.jl")]
    system=source[joinpath("src", "parametricbuilder", "system.jl")]
    @test !occursin("_position_kind", positions)
    @test !occursin("_position_kind", system)
    @test !occursin("shell.kind", makie_source)

    backends=source[joinpath("src", "plotbuilder", "backends.jl")]
    @test !occursin("_parent_package", backends)
    @test !hasmethod(LineCableModels.PlotBuilder.ensure_backend!, Tuple{Symbol})
    @test !hasmethod(LineCableModels.PlotBuilder.ensure_backend!, Tuple{Val{:cairo}})

    @test !isdefined(LineCableModels.Engine, :FEM)
    @test !isdefined(LineCableModels.DataModel, :SectorParams)
    @test !isdefined(LineCableModels.DataModel, :SectorInsulator)

    root_project=TOML.parsefile(joinpath(root, "Project.toml"))
    @test root_project["workspace"]["projects"] == ["test"]
    test_project=TOML.parsefile(joinpath(root, "test", "Project.toml"))
    @test test_project["sources"]["LineCableModels"]["path"] == ".."
    @test !haskey(root_project, "extras")
    @test !haskey(root_project, "targets")

    for pattern in (
        r"details\s*::\s*Dict",
        r"details\s*::\s*Any",
        r"function\s+computation_details\s*\([^,)]*\)"
    )
        @test all(!occursin(pattern, contents) for contents in values(source))
    end

    presentation_paths=filter(keys(source)) do path
        startswith(path, joinpath("src", "plotbuilder")) ||
            startswith(path, joinpath("src", "reportbuilder")) ||
            endswith(path, "plotdefinition.jl") ||
            endswith(path, "comparisonplot.jl") ||
            endswith(path, joinpath("montecarlo", "plot.jl"))
    end
    for path in presentation_paths
        contents=source[path]
        @test !occursin(r"\.(Z|Y|f)\b", contents)
        @test !occursin(r"\b(Ω/m|Ω/km|H/m|mH/km|S/m|S/km|F/m|μF/km)\b", contents)
    end

    @test all(!occursin("@eval", contents) for contents in values(source))
    @test all(!occursin(r"^\s*(using|import)\s+.*Makie"m, contents)
    for contents in values(source))
    @test all(!occursin("@eval", read(path, String))
    for (directory, _, files) in walkdir(extension_root)
    for file in files if endswith(file, ".jl")
    for path in (joinpath(directory, file),))
end

@testitem "ParametricBuilder / relax / isolated future provision" tags=[:unit] begin
    using LineCableModels.ParametricBuilder: @relax

    @relax struct FutureValue{T <: Real}
        scalar::T
        values::Vector{T}
        label::Symbol
    end

    value=FutureValue(1.0f0, [2.0, 3.0], :future)
    @test value isa FutureValue{Float64}
    @test value.scalar == 1.0
    @test value.values == [2.0, 3.0]
    converted=convert(FutureValue{Float32}, value)
    @test converted isa FutureValue{Float32}
    @test converted.values == Float32[2, 3]
    @test converted.label === :future
end
