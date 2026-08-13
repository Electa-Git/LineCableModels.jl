using Documenter

deploydocs(;
    repo = "github.com/Electa-Git/LineCableModels.jl.git",
    devbranch = "main",
    versions = ["stable" => "v^", "dev" => "main"],
    branch = "gh-pages"
)
