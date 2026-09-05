include("artifacts.jl")
using .GauntletArtifacts

length(ARGS) == 3 || throw(ArgumentError(
    "usage: bind.jl COLLECTION VERSION PUBLISHED_URL",
))

binding = bind_published_artifact(
    Symbol(ARGS[1]),
    VersionNumber(ARGS[2]),
    ARGS[3]
)

println("artifact=$(binding.artifact)")
println("version=$(binding.version)")
println("tag=$(binding.tag)")
println("tree_hash=$(binding.tree_hash)")
println("archive_sha256=$(binding.archive_sha256)")
println("url=$(binding.url)")
