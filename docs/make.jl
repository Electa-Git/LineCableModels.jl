using Documenter, DocumenterCitations, Literate, LineCableModels, Measurements, Gmsh
using TOML

const GmshExt = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
const MeasurementsExt = Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt)
const REPOSITORY_URL = "https://github.com/Electa-Git/LineCableModels.jl"
const SITE_URL = "https://electa-git.github.io/LineCableModels.jl"

project = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))

# The examples, in the order of the contents, rendered from their scripts at every build.
example_files = ["solid_1000mm2_single.jl", "cable_380kv_no_armour_ac_flat.jl",
    "cable_320kv_armoured_dc_bipole.jl", "cable_525kv_subsea_armoured_ac_flat.jl"]
rm(joinpath(@__DIR__, "src", "examples"); recursive = true, force = true)
for file in example_files
    Literate.markdown(joinpath(@__DIR__, "..", "examples", file),
        joinpath(@__DIR__, "src", "examples"); documenter = true, credit = false)
end
examples = [String(chopprefix(readline(joinpath(@__DIR__, "..", "examples", file)),
                "# # ")) => joinpath("examples", replace(file, ".jl" => ".md"))
            for file in example_files]
example_pages = last.(examples)

DocMeta.setdocmeta!(LineCableModels, :DocTestSetup, :(using LineCableModels);
    recursive = true)

bibliography = CitationBibliography(joinpath(@__DIR__, "src", "bibliography.bib");
    style = :numeric)

makedocs(;
    modules = [LineCableModels, GmshExt, MeasurementsExt],
    authors = join(project["authors"], ", "),
    sitename = "LineCableModels.jl",
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
        footer = "[LineCableModels.jl]($(REPOSITORY_URL)) v$(project["version"]) supported by the Etch Competence Hub of EnergyVille, financed by the Flemish Government.",
        size_threshold_warn = 700 * 1024,
        size_threshold = 1024 * 1024
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => Any["Contents" => "examples.md", examples...],
        "User guide" => Any[
            "Cable data model" => "data-model.md",
            "Line parameters" => "computing.md",
            "Finite element solver" => "fem.md",
            "Uncertainty quantification" => "uncertainty.md"
        ],
        "Theory" => Any[
            "Contents" => "theory/contents.md",
            "Matrix formulation" => "theory/matrix.md",
            "Conductors" => "theory/conductors.md",
            "Dielectrics" => "theory/dielectrics.md",
            "Earth return" => "theory/earth.md"
        ],
        "API reference" => "reference.md",
        "Developers" => Any[
            "Architecture" => "developers/architecture.md",
            "Rules and patterns" => "developers/rules.md",
            "Conventions" => "developers/conventions.md"
        ],
        "Bibliography" => "bibliography.md"
    ],
    plugins = [bibliography],
    checkdocs = :public,
    doctest = true,
    pagesonly = true,
    clean = true
)
