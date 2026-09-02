using Documenter
using DocumenterCitations
using CairoMakie
using JLD2
using LineCableModels
using Literate
using Markdown
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
const CASE_OUTPUT = joinpath(DOCS_SRC_DIR, "cases")
const CASE_ASSETS = joinpath(DOCS_SRC_DIR, "assets", "cases")
const CASE_MANIFEST = joinpath(@__DIR__, ".generated", "case_catalogue.jld2")
const LCM_CLI = joinpath(ROOT_DIR, "cli", "lcm")

const CONVENIENCE_API_OBJECTS = ()

const BENCHMARK_API_OBJECTS = (
    LineCableModels.Engine.RMSError,
    LineCableModels.Engine.LineParametersBenchmark,
    LineCableModels.Engine.compare,
    LineCableModels.Engine.absolute_error,
    LineCableModels.Engine.relative_error
)

const EXTENSION_API_OBJECTS = (
    LineCableModels.InputValidation,
    LineCableModels.InputValidation.validate,
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
    LineCableModels.ImportExport.deserialize_extension
)

_contains_identity(collection, object) = any(candidate -> candidate === object, collection)
function api_reference_entry(object)
    !_contains_identity(CONVENIENCE_API_OBJECTS, object) &&
        !_contains_identity(BENCHMARK_API_OBJECTS, object) &&
        !_contains_identity(EXTENSION_API_OBJECTS, object)
end
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

function formulation_catalogue(module_owner::Module, path::AbstractString...)
    directory = normpath(joinpath(ROOT_DIR, "src", path...))
    entries = NamedTuple[]

    for (binding, multidoc) in Base.Docs.meta(module_owner)
        binding.var === :description || continue
        for (typesig, docstring) in multidoc.docs
            source = normpath(String(docstring.data[:path]))
            dirname(source) == directory || continue

            formula_type = only(typesig.parameters)
            identifier = first(Base.unwrap_unionall(formula_type).parameters)
            body = join(filter(
                value -> value isa AbstractString,
                collect(docstring.text)
            ))
            push!(entries, (; source, identifier, body = strip(body)))
        end
    end

    sort!(entries; by = entry -> basename(entry.source))
    isempty(entries) && error("no formulation docstrings found in $directory")
    return Markdown.parse(join(
        ("### `:$(entry.identifier)`\n\n$(entry.body)" for entry in entries),
        "\n\n"
    ))
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

function markdown_cell(value)
    return replace(string(value), '|' => "\\|", '\n' => "<br>")
end

function markdown_table(rows)
    isempty(rows) && return "_None._\n"
    names = collect(keys(first(rows)))
    lines = String[
        "| " * join(string.(names), " | ") * " |",
        "| " * join(fill("---", length(names)), " | ") * " |"
    ]
    for row in rows
        push!(lines, "| " *
                     join(
                         (markdown_cell(getproperty(row, name)) for name in names),
                         " | "
                     ) * " |")
    end
    return join(lines, '\n') * "\n"
end

function case_page(record)
    id = record.id
    problem = LineCableModels.ImportExport.deserialize_value(record.problem)
    system = problem.system
    design_image = "$(id)-designs.png"
    system_image = "$(id)-system.png"

    design_plot = preview(
        system.designs;
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    system_plot = preview(
        system;
        backend = :cairo,
        display_plot = false,
        controls = false
    )
    CairoMakie.save(joinpath(CASE_ASSETS, design_image), design_plot.figure)
    CairoMakie.save(joinpath(CASE_ASSETS, system_image), system_plot.figure)

    parameters = if isempty(record.parameters)
        "_This is a fixed materialised case; it declares no variable case parameters._\n"
    else
        markdown_table([(
                            id = parameter.id,
                            nominal = parameter.nominal,
                            tags = join(string.(parameter.tags), ", ")
                        ) for parameter in record.parameters])
    end
    problem_text = sprint(show, MIME("text/plain"), problem)
    return """
# `$(id)`

$(record.description)

Input fingerprint: `$(record.input_sha256)`.

## Summary

$(markdown_table(record.summary))

## Cable designs

$(markdown_table(record.designs))

![Cable-design cross-sections](../assets/cases/$(design_image))

## System cross-section

The horizontal reference line is the air–earth interface declared by the line-parameter geometry.

![Cable-system cross-section](../assets/cases/$(system_image))

## Case parameters

$(parameters)

## Materialised problem

```text
$(problem_text)
```

[Back to the case catalogue](index.md).
"""
end

function build_case_catalogue!()
    mkpath(dirname(CASE_MANIFEST))
    run(`$LCM_CLI gauntlet case catalogue --output $CASE_MANIFEST`)
    document = JLD2.load(CASE_MANIFEST)
    document["schema_version"] == 1 || error("unsupported Gauntlet case catalogue")
    records = document["cases"]

    rm(CASE_OUTPUT; recursive = true, force = true)
    rm(CASE_ASSETS; recursive = true, force = true)
    mkpath(CASE_OUTPUT)
    mkpath(CASE_ASSETS)

    index_rows = NamedTuple[]
    pages = Pair{String, String}[]
    for record in records
        id = record.id
        path = joinpath(CASE_OUTPUT, "$(id).md")
        write(path, case_page(record))
        push!(pages, id => joinpath("cases", "$(id).md"))
        push!(index_rows,
            (
                case = "[`$(id)`]($(id).md)",
                description = record.description,
                terminals = length(record.port_order),
                frequencies = length(
                    LineCableModels.ImportExport.deserialize_value(record.problem).frequencies
                )
            ))
    end
    write(
        joinpath(CASE_OUTPUT, "index.md"),
        """# Gauntlet case catalogue

This catalogue is generated from the indexed, materialised Gauntlet cases. It does not run a numerical formulation.

$(markdown_table(index_rows))

The input fingerprint covers the complete materialised `LineParametersProblem`; result artefacts separately record the selected numerical implementation blobs.
"""
    )
    return ["Contents" => "cases/index.md"; pages]
end

metadata = project_metadata()
tutorials = build_tutorials!()
tutorial_pages = last.(tutorials)
generate_maintained_pages!()
case_pages = build_case_catalogue!()

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
            "Transmission line parameters" => "transmission-line-parameters.md",
            "Gmsh/GetDP FEM backend" => "fem.md",
            "Gridspace and uncertainty" => "gridspace.md"
        ],
        "API reference" => "reference.md",
        "Conveniences" => Any[
            "Overview" => "conveniences.md",
            "Data entry validation" => "validation.md"
        ],
        "Benchmarks" => Any[
            "Gauntlet" => "gauntlet.md",
            "Case catalogue" => case_pages
        ],
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
