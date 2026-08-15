using Changelog
using Documenter
using DocumenterCitations
using LineCableModels
using Literate
using TOML

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DOCS_SRC_DIR = joinpath(@__DIR__, "src")
const REPOSITORY = "Electa-Git/LineCableModels.jl"
const REPOSITORY_URL = "https://github.com/$(REPOSITORY)"
const CANONICAL_URL = "https://electa-git.github.io/LineCableModels.jl"
const TUTORIAL_SOURCE = joinpath(ROOT_DIR, "examples")
const TUTORIAL_OUTPUT = joinpath(DOCS_SRC_DIR, "tutorials")
const PLOTBUILDER_SOURCE = joinpath(@__DIR__, "literate", "plotbuilder.jl")

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
        return String(matchobj.captures[1])
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
            documenter = false,
            postprocess = strip_literate_footer
        )
    end

    files = sort(filter(file -> endswith(file, ".md"), readdir(TUTORIAL_OUTPUT)))
    return [tutorial_title(joinpath(TUTORIAL_OUTPUT, file)) => joinpath("tutorials", file)
            for
            file in files]
end

function generate_maintained_pages!()
    Changelog.generate(
        Changelog.Documenter(),
        joinpath(ROOT_DIR, "CHANGELOG.md"),
        joinpath(DOCS_SRC_DIR, "CHANGELOG.md");
        repo = REPOSITORY
    )
    cp(joinpath(ROOT_DIR, "TODO.md"), joinpath(DOCS_SRC_DIR, "TODO.md"); force = true)
    Literate.markdown(
        PLOTBUILDER_SOURCE,
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
        canonical = CANONICAL_URL,
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
        "API reference" => "reference.md",
        "Development" => Any[
            "Conventions" => "conventions.md",
            "PlotBuilder guide" => "plotbuilder.md",
            "Validation module" => "validation.md",
            "Docstrings" => "docstrings.md",
            "Contributing" => "contributing.md",
            "TODO" => "TODO.md",
            "Changelog" => "CHANGELOG.md"
        ],
        "Bibliography" => "bib.md"
    ],
    clean = true,
    plugins = [bibliography],
    checkdocs = :exports,
    pagesonly = true
)

@info "Finished documentation build."
