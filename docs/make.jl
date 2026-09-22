using Documenter
using DocumenterCitations
using CairoMakie
using LineCableModels
using Literate
using TOML

include("type_trees.jl")
include("owned_doc_links.jl")

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DOCS_SRC_DIR = joinpath(@__DIR__, "src")
const REPOSITORY = "Electa-Git/LineCableModels.jl"
const REPOSITORY_URL = "https://github.com/$(REPOSITORY)"
const SITE_URL = "https://electa-git.github.io/LineCableModels.jl"
const TUTORIAL_SOURCE = joinpath(ROOT_DIR, "examples")
const TUTORIAL_OUTPUT = joinpath(DOCS_SRC_DIR, "tutorials")
const PLOTTING_SOURCE = joinpath(@__DIR__, "literate", "plotting.jl")

const CONVENIENCE_API_OBJECTS = ()

const EXTENSION_API_OBJECTS = (
    LineCableModels.InputValidation,
    LineCableModels.InputValidation.validate,
    LineCableModels.Grammar.check_core_result,
    LineCableModels.Grammar.validate_observables,
    LineCableModels.Grammar.unit_targets,
    LineCableModels.Grammar.detach,
    LineCableModels.Grammar.observation_request,
    LineCableModels.Grammar.observation_indices,
    LineCableModels.Grammar.materialize_observation,
    LineCableModels.Grammar.request_identity,
    LineCableModels.Grammar.request_quantity,
    LineCableModels.Grammar.request_indices,
    LineCableModels.ObservedResult,
    LineCableModels.Grammar.observation_quantity,
    LineCableModels.Grammar.observation_gridpoint,
    LineCableModels.Grammar.observation_groups,
    LineCableModels.Units.family,
    LineCableModels.DataModel.preview_shapes,
    LineCableModels.DataModel.preview_materials,
    LineCableModels.DataModel.PreviewShape,
    LineCableModels.DataModel.material_property_ranges,
    LineCableModels.materialcolors,
    LineCableModels.materialscale!,
    LineCableModels.Engine.has_uncertainty_type,
    LineCableModels.materialize,
    LineCableModels.ParametricBuilder.traverse,
    LineCableModels.sample_uncertainty,
    LineCableModels.UIPlot,
    LineCableModels.plotwindow,
    LineCableModels.export_svg,
    LineCableModels.figurelegend!,
    LineCableModels.panellegend!,
    LineCableModels.figuretitle!,
    LineCableModels.paneltitle!,
    LineCableModels.figurecolorbars!,
    LineCableModels.axisscale!,
    LineCableModels.resetview!,
    LineCableModels.addwidget!,
    LineCableModels.removewidget!,
    LineCableModels.ReportBuilder,
    LineCableModels.ReportBuilder.AbstractReportDefinition,
    LineCableModels.ReportBuilder.CableConstantsTableDefinition,
    LineCableModels.ReportBuilder.LineParametersTableDefinition,
    LineCableModels.ReportBuilder.BenchmarkTableDefinition,
    LineCableModels.ReportBuilder.MonteCarloTableDefinition,
    LineCableModels.ReportBuilder.select,
    LineCableModels.ReportBuilder.tabulate,
    LineCableModels.ReportBuilder.illustrate,
    LineCableModels.ReportBuilder.encode,
    LineCableModels.ReportBuilder.write,
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
    Literate.markdown(
        PLOTTING_SOURCE,
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

bibliography = CitationBibliography(joinpath(DOCS_SRC_DIR, "bibliography.bib"); style = :numeric)
owned_doc_links = OwnedDocLinks.OwnedDocLinker(LineCableModels)

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
        "Theory" => Any[
            "Contents" => "theory/contents.md",
            "Matrix formulation" => "theory/matrix_formulation.md",
            "Modal decomposition" => Any[
                "Overview" => "theory/modal_decomposition.md",
                "Default modal decomposition and eigenvalue tracking" =>
                    "theory/modal-decomposition/default.md"
            ],
            "Earth properties" => Any[
                "Overview" => "theory/earth_properties.md",
                "Default frequency-dependent earth material" =>
                    "theory/earth-properties/frequency-dependent/default.md",
                "Registered frequency-dependent soil relations" =>
                    "theory/earth-properties/frequency-dependent/formulas.md",
                "Default equivalent homogeneous-earth rule" =>
                    "theory/earth-properties/equivalent-homogeneous/default.md"
            ],
            "Earth return admittance" => Any[
                "Overview" => "theory/earth_return_admittance.md",
                "Default two-half-space earth-return admittance" =>
                    "theory/external-admittance/default.md",
                "Pollaczek underground earth-return admittance" =>
                    "theory/external-admittance/1926/homogeneous-earth-generalized-induction-green-function/Pollaczek1926.md",
                "Wise homogeneous-earth overhead potential coefficient" =>
                    "theory/external-admittance/1948/homogeneous-earth-overhead-potential-coefficient/Wise1948.md",
                "Xue underground earth-return admittance" =>
                    "theory/external-admittance/2018/complete-field-and-quasi-tem-underground/Xue2018.md"
            ],
            "Earth return impedance" => Any[
                "Overview" => "theory/earth_return_impedance.md",
                "Default two-half-space earth-return impedance" =>
                    "theory/external-impedance/default.md",
                "Carson homogeneous-earth overhead correction integral" =>
                    "theory/external-impedance/1926/homogeneous-earth-overhead-integral/Carson1926.md",
                "Gary complex-depth overhead approximation" =>
                    "theory/external-impedance/1976/complex-depth-overhead/Gary1976.md",
                "Lucca mixed-pair homogeneous-earth impedance" =>
                    "theory/external-impedance/1994/mixed-pair/Lucca1994.md",
                "Pollaczek generalized induction coefficients" =>
                    "theory/external-impedance/1926/homogeneous-earth-generalized-induction-green-function/Pollaczek1926.md",
                "Saad homogeneous-earth underground closed form" =>
                    "theory/external-impedance/1996/homogeneous-earth-underground-closed-form/Saad1996.md",
                "Wise high-frequency overhead displacement-current integral" =>
                    "theory/external-impedance/1934/homogeneous-earth-overhead-displacement-current-integral/Wise1934.md",
                "Wedepohl–Wilcox underground low-order impedance" =>
                    "theory/external-impedance/1973/wedepohl-wilcox-low-order/WedepohlWilcox1973.md",
                "Xue underground earth-return impedance" =>
                    "theory/external-impedance/2018/complete-field-and-quasi-tem-underground/Xue2018.md",
                "Ametani mixed-pair exponential-image approximation" =>
                    "theory/external-impedance/2009/homogeneous-earth-mixed-exponential-image/Ametani2009.md"
            ],
            "Insulation parameters" => Any[
                "Overview" => "theory/insulation_parameters.md",
                "Default cable-insulation admittivity" =>
                    "theory/insulation-admittance/default.md",
                "Lossless cable-layer admittivity" =>
                    "theory/insulation-admittance/lossless.md",
                "Lossy cable-layer admittivity (Ametani 2004 application)" =>
                    "theory/insulation-admittance/2004/semiconducting-screen-complex-permittivity/Ametani2004.md",
                "Default coaxial-insulation magnetic series impedance" =>
                    "theory/insulation-impedance/default.md"
            ],
            "Internal impedance" => Any[
                "Overview" => "theory/internal_impedance.md",
                "Default cylindrical-conductor surface impedances" =>
                    "theory/internal-impedance/default.md",
                "Default analytical pipe-type treatment" =>
                    "theory/internal-impedance/pipe-default.md"
            ]
        ],
        "Tutorials" => Any["Contents" => "tutorials.md", tutorials...],
        "User guide" => Any[
            "Cable data model" => "data-model.md",
            "Modeling and results" => "usage.md",
            "Gmsh/GetDP FEM backend" => "fem.md",
            "Gridspace and uncertainty" => "gridspace.md"
        ],
        "API reference" => "reference.md",
        "Conveniences" => Any[
            "Overview" => "conveniences.md",
            "Data entry validation" => "validation.md"
        ],
        "Developers" => Any[
            "Grammar invariants" => "developers.md",
            "Extension API" => "extensions.md",
            "Conventions" => "conventions.md",
            "Computational engine" => "engine.md",
            "Makie plotting" => "plotting.md",
            "Contributing" => "contributing.md"
        ],
        "Bibliography" => "bibliography.md"
    ],
    clean = true,
    plugins = [bibliography, owned_doc_links],
    checkdocs = :exports,
    pagesonly = true
)

owned_doc_links.linked > 0 || error("owned documentation linker did not process any names")

@info "Finished documentation build."
