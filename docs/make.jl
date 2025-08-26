using Documenter
using DocumenterCitations
using Literate
using Pkg
using Changelog

function get_project_toml()
    # Get the current active environment (docs)
    docs_env = Pkg.project().path

    # Path to the main project (one level up from docs)
    main_project_path = joinpath(dirname(docs_env), "..")

    # Parse the main project's TOML
    project_toml = Pkg.TOML.parsefile(joinpath(main_project_path, "Project.toml"))

    return project_toml
end

function open_in_default_browser(url::AbstractString)::Bool
    try
        if Sys.isapple()
            Base.run(`open $url`)
            true
        elseif Sys.iswindows()
            Base.run(`powershell.exe Start "'$url'"`)
            true
        elseif Sys.islinux()
            Base.run(`xdg-open $url`, devnull, devnull, devnull)
            true
        else
            false
        end
    catch ex
        false
    end
end



# Get project data
PROJECT_TOML = get_project_toml()
PROJECT_VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors."
GITHUB = PROJECT_TOML["git_url"]

@eval using $(Symbol(NAME))
main_module = @eval $(Symbol(NAME))

function customize_literate_footer(content, custom_footer="")
    if isempty(custom_footer)
        return replace(
            content,
            r"---\s*\n\*This page was generated using \[Literate\.jl\]\(.*?\)\.\*\s*$" => "",
        )
    else
        return replace(content,
            r"\*This page was generated using \[Literate\.jl\]\(.*?\)\.\*\s*$" =>
                custom_footer)
    end
end

function post_process_literate(content)
    content = customize_literate_footer(content, "🏠 Back to [Tutorials](@ref)\n")
    return content
end

tutorial_source = joinpath(@__DIR__, "..", "examples")
tutorial_output = joinpath(@__DIR__, "src", "tutorials")
# Remove the directory if it exists and then create it fresh
if isdir(tutorial_output)
    rm(tutorial_output, recursive=true)
end
mkpath(tutorial_output)

for file in readdir(tutorial_source)
    if endswith(file, ".jl")
        Literate.markdown(
            joinpath(tutorial_source, file),
            tutorial_output,
            documenter=true,
            postprocess=content ->
                post_process_literate(
                    content,
                ),
        )
    end
end

# Get all .md files in tutorial_output
tutorial_files = filter(
    file -> endswith(file, ".md") && file != "index.md",
    readdir(tutorial_output),
)

# Build menu from existing files only
tutorial_menu = ["Contents" => "tutorials.md"]
for file in tutorial_files
    relative_path = String(joinpath("tutorials", file))  # Convert to full String
    # Get title from file content
    content = read(joinpath(tutorial_output, file), String)
    m = match(r"#\s+(.*)", content)
    # Make sure title is a full String too, not SubString
    if m !== nothing
        title = String(m.captures[1])
    else
        title = String(titlecase(replace(basename(file)[1:end-3], "_" => " ")))
    end
    push!(tutorial_menu, title => relative_path)
end

tutorial_pages = [String(joinpath("tutorials", file)) for file in tutorial_files]

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib"),
    style=:numeric,  # default
)

DocMeta.setdocmeta!(
    main_module,
    :DocTestSetup,
    :(using $(Symbol(NAME)));
    recursive=true,
)

mathengine = MathJax3(
    Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
            "tags" => "ams",
            "packages" => ["base", "ams", "autoload", "physics"],
        ),
        :chtml => Dict(
            :scale => 1.1,
        ),
    ),
)

Changelog.generate(
    Changelog.Documenter(),                 # output type
    joinpath(@__DIR__, "..", "CHANGELOG.md"),  # input file
    joinpath(@__DIR__, "src", "CHANGELOG.md"); # output file
    repo="Electa-Git/LineCableModels.jl",        # default repository for links
)

todo_src = joinpath(@__DIR__, "..", "TODO.md")
todo_dest = joinpath(@__DIR__, "src", "TODO.md")
cp(todo_src, todo_dest, force=true)

makedocs(;
    modules=[main_module],
    authors="Amauri Martins",
    sitename="$NAME.jl",
    format=Documenter.HTML(;
        mathengine=mathengine,
        edit_link="main",
        assets=[
            "assets/citations.css",
            "assets/favicon.ico",
            "assets/custom.css",
            "assets/custom.js",
        ],
        prettyurls=get(ENV, "CI", "false") == "true",
        ansicolor=true,
        collapselevel=1,
        footer="[$NAME.jl]($GITHUB) v$PROJECT_VERSION supported by the Etch Competence Hub of EnergyVille, financed by the Flemish Government.",
        size_threshold=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => tutorial_menu,
        "API reference" => "reference.md",
        "Development" => Any[
            "Validation module"=>"validation.md",
            "Docstrings"=>"docstrings.md",
            "TODO"=>"TODO.md",
            "Changelog"=>"CHANGELOG.md",
        ],
        "Bibliography" => "bib.md",
    ],
    clean=true,
    plugins=[bib],
    checkdocs=:exports,
    pagesonly=true,
    warnonly=true,
)

if haskey(ENV, "CI")
    deploydocs(
        repo="github.com/Electa-Git/LineCableModels.jl.git",
        devbranch="main",
        versions=["stable" => "v^", "dev" => "main"],
        branch="gh-pages",
    )
else
    open_in_default_browser(
        "file://$(abspath(joinpath(@__DIR__, "build", "index.html")))",
    ) ||
        println("Failed to open the documentation in the browser.")
end
@info "Finished docs build." # Good to know the script completed

