using Documenter
using DocumenterCitations
using LineCableModels
using Literate
using TOML

include("type_trees.jl")

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DOCS_SRC_DIR = joinpath(@__DIR__, "src")
const REPOSITORY = "Electa-Git/LineCableModels.jl"
const REPOSITORY_URL = "https://github.com/$(REPOSITORY)"
const SITE_URL = "https://electa-git.github.io/LineCableModels.jl"
const TUTORIAL_SOURCE = joinpath(ROOT_DIR, "examples")
const TUTORIAL_OUTPUT = joinpath(DOCS_SRC_DIR, "tutorials")
const PLOTBUILDER_SOURCE = joinpath(@__DIR__, "literate", "plotbuilder.jl")
const GAUNTLET_SOURCE = joinpath(@__DIR__, "literate", "gauntlet.jl")

const CONVENIENCE_API_OBJECTS = ()

const BENCHMARK_API_OBJECTS = (
    LineCableModels.Engine.RMSError,
    LineCableModels.Engine.LineParametersBenchmark,
    LineCableModels.Engine.compare,
    LineCableModels.Engine.absolute_error,
    LineCableModels.Engine.relative_error,
)

const EXTENSION_API_OBJECTS = (
    LineCableModels.Validation,
    LineCableModels.Validation.Rule,
    LineCableModels.Validation.OwnerRule,
    LineCableModels.Validation.rules,
    LineCableModels.Validation.check,
    getfield(LineCableModels.Grammar, Symbol("@orchestrator")),
    LineCableModels.Grammar.check_core_result,
    LineCableModels.Grammar.computation_owner,
    LineCableModels.Grammar.validate_observables,
    LineCableModels.Grammar.unit_targets,
    LineCableModels.Grammar.detach,
    LineCableModels.Grammar.observation_request,
    LineCableModels.Grammar.observation_indices,
    LineCableModels.Grammar.materialize_observation,
    LineCableModels.Grammar.request_identity,
    LineCableModels.Grammar.request_quantity,
    LineCableModels.Grammar.request_indices,
    LineCableModels.Grammar.ObservationPublication,
    LineCableModels.Grammar.publication_table,
    LineCableModels.Grammar.orchestrator_root,
    LineCableModels.Grammar.orchestrator_method,
    LineCableModels.Units.family,
    LineCableModels.DataModel.preview_shapes,
    LineCableModels.DataModel.preview_materials,
    LineCableModels.DataModel.PreviewPolygon,
    LineCableModels.DataModel.PreviewReferenceLine,
    LineCableModels.DataModel.PreviewPayload,
    LineCableModels.Engine.has_uncertainty_type,
    LineCableModels.ParametricBuilder.AbstractProjectionDefinition,
    LineCableModels.ParametricBuilder.entitle,
    LineCableModels.ParametricBuilder.select,
    LineCableModels.ParametricBuilder.derive,
    LineCableModels.ParametricBuilder.materialize,
    LineCableModels.ParametricBuilder.traverse,
    LineCableModels.ParametricBuilder.sample_uncertainty,
    LineCableModels.PlotBuilder,
    LineCableModels.PlotBuilder.AbstractPlotDefinition,
    LineCableModels.PlotBuilder.PlotPage,
    LineCableModels.PlotBuilder.PlotRecipe,
    LineCableModels.PlotBuilder.LegendDefinition,
    LineCableModels.PlotBuilder.ColorbarDefinition,
    LineCableModels.PlotBuilder.ExportDefinition,
    LineCableModels.PlotBuilder.AbstractWidgetDefinition,
    LineCableModels.PlotBuilder.make_render,
    LineCableModels.PlotBuilder.plotwindow,
    LineCableModels.PlotBuilder.axis!,
    LineCableModels.PlotBuilder.register!,
    LineCableModels.PlotBuilder.backend_available,
    LineCableModels.PlotBuilder.current_backend_symbol,
    LineCableModels.PlotBuilder.ensure_backend!,
    LineCableModels.PlotBuilder.make_screen,
    LineCableModels.PlotBuilder.next_fignum,
    LineCableModels.PlotBuilder.renderfig,
    LineCableModels.PlotBuilder.with_backend,
    LineCableModels.PlotBuilder.input_defaults,
    LineCableModels.PlotBuilder.renderer_defaults,
    LineCableModels.PlotBuilder.entitle,
    LineCableModels.PlotBuilder.parse,
    LineCableModels.PlotBuilder.resolve,
    LineCableModels.PlotBuilder.fetch,
    LineCableModels.PlotBuilder.finish,
    LineCableModels.PlotBuilder.validate_export_theme,
    LineCableModels.ReportBuilder,
    LineCableModels.ReportBuilder.AbstractReportDefinition,
    LineCableModels.ReportBuilder.CableConstantsTableDefinition,
    LineCableModels.ReportBuilder.LineParametersTableDefinition,
    LineCableModels.ReportBuilder.BenchmarkTableDefinition,
    LineCableModels.ReportBuilder.MonteCarloTableDefinition,
    LineCableModels.ReportBuilder.entitle,
    LineCableModels.ReportBuilder.select,
    LineCableModels.ReportBuilder.tabulate,
    LineCableModels.ReportBuilder.illustrate,
    LineCableModels.ReportBuilder.encode,
    LineCableModels.ReportBuilder.write,
    LineCableModels.ReportBuilder.finish,
    LineCableModels.ReportBuilder.observation_columns,
    LineCableModels.ReportBuilder.encode_cell,
    LineCableModels.ReportBuilder.XLSXSheet,
    LineCableModels.ReportBuilder.XLSXWorkbook,
    LineCableModels.ImportExport.serialize_value,
    LineCableModels.ImportExport.deserialize_value,
    LineCableModels.ImportExport.deserialize_extension,
)

_contains_identity(collection, object) = any(candidate -> candidate === object, collection)
api_reference_entry(object) =
    !_contains_identity(CONVENIENCE_API_OBJECTS, object) &&
    !_contains_identity(BENCHMARK_API_OBJECTS, object) &&
    !_contains_identity(EXTENSION_API_OBJECTS, object)
developer_reference_entry(object) = _contains_identity(EXTENSION_API_OBJECTS, object)

function project_metadata()
    project = TOML.parsefile(joinpath(ROOT_DIR, "Project.toml"))
    authors = get(project, "authors", String[])
    return (
        name = get(project, "name", "LineCableModels"),
        version = get(project, "version", "dev"),
        authors = isempty(authors) ? "LineCableModels contributors" : join(authors, ", ")
    )
end

function strip_literate_footer(content::AbstractString)
    return replace(
        content,
        r"(?ms)^---\s*\n\*This page was generated using \[Literate\.jl\]\(.*?\)\.\*\s*$" => "Back to [Tutorials](../tutorials.md)\n"
    )
end

normalize_literate_page(content::AbstractString) = rstrip(content) * "\n"

function tutorial_title(path::AbstractString)
    content = read(path, String)
    matchobj = match(r"(?m)^#\s+(.+)$", content)
    isnothing(matchobj) ||
        return replace(String(matchobj.captures[1]), r"^#+\s*" => "")
    stem = splitext(basename(path))[1]
    return titlecase(replace(stem, "_" => " ", "-" => " "))
end

function build_tutorials!()
    rm(TUTORIAL_OUTPUT; recursive = true, force = true)
    mkpath(TUTORIAL_OUTPUT)

    for file in sort(readdir(TUTORIAL_SOURCE))
        endswith(file, ".jl") || continue
        Literate.markdown(
            joinpath(TUTORIAL_SOURCE, file),
            TUTORIAL_OUTPUT;
            documenter = true,
            postprocess = strip_literate_footer
        )
    end

    files = sort(filter(file -> endswith(file, ".md"), readdir(TUTORIAL_OUTPUT)))
    return [tutorial_title(joinpath(TUTORIAL_OUTPUT, file)) => joinpath("tutorials", file)
            for
            file in files]
end

function generate_maintained_pages!()
    cp(joinpath(ROOT_DIR, "TODO.md"), joinpath(DOCS_SRC_DIR, "TODO.md"); force = true)
    Literate.markdown(
        PLOTBUILDER_SOURCE,
        DOCS_SRC_DIR;
        documenter = true,
        credit = false,
        postprocess = normalize_literate_page
    )
    Literate.markdown(
        GAUNTLET_SOURCE,
        DOCS_SRC_DIR;
        documenter = true,
        credit = false,
        postprocess = normalize_literate_page
    )
    return nothing
end

metadata = project_metadata()
tutorials = build_tutorials!()
tutorial_pages = last.(tutorials)
generate_maintained_pages!()

DocMeta.setdocmeta!(
    LineCableModels,
    :DocTestSetup,
    quote
        using LineCableModels
        using LineCableModels.DataModel.BaseParams
    end;
    recursive = true
)

bibliography = CitationBibliography(joinpath(DOCS_SRC_DIR, "refs.bib"); style = :numeric)

makedocs(;
    modules = [LineCableModels],
    authors = metadata.authors,
    sitename = "$(metadata.name).jl",
    format = Documenter.HTML(;
        canonical = SITE_URL,
        edit_link = "main",
        assets = [
            "assets/citations.css",
            "assets/favicon.ico",
            "assets/custom.css",
            "assets/custom.js"
        ],
        mathengine = MathJax3(
            Dict(
            :loader => Dict("load" => ["[tex]/physics"]),
            :tex => Dict(
                "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                "tags" => "ams",
                "packages" => ["base", "ams", "autoload", "physics"]
            ),
            :chtml => Dict(:scale => 1.1)
        ),
        ),
        prettyurls = get(ENV, "CI", "false") == "true",
        footer = "[$(metadata.name).jl]($(REPOSITORY_URL)) v$(metadata.version) supported by the Etch Competence Hub of EnergyVille, financed by the Flemish Government.",
        size_threshold_warn = 700 * 1024,
        size_threshold = 1024 * 1024
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any["Contents" => "tutorials.md", tutorials...],
        "User guide" => Any[
            "Cable data model" => "data-model.md",
            "Modelling and results" => "usage.md",
            "Gmsh/GetDP FEM backend" => "fem.md",
            "Gridspace and uncertainty" => "gridspace.md",
            "Global sensitivity" => "sensitivity.md"
        ],
        "API reference" => "reference.md",
        "Conveniences" => Any[
            "Overview" => "conveniences.md",
            "Data entry validation" => "validation.md"
        ],
        "Benchmarks" => "gauntlet.md",
        "Developers" => Any[
            "Grammar invariants" => "developers.md",
            "Extension API" => "extensions.md",
            "Conventions" => "conventions.md",
            "Computational engine" => "engine.md",
            "PlotBuilder guide" => "plotbuilder.md",
            "Contributing" => "contributing.md",
            "TODO" => "TODO.md"
        ],
        "Bibliography" => "bib.md"
    ],
    clean = true,
    plugins = [bibliography],
    checkdocs = :exports,
    pagesonly = true
)

@info "Finished documentation build."
