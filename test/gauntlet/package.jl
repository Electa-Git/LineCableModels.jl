include("artifacts.jl")
using .GauntletArtifacts

length(ARGS) in (4, 5) || throw(ArgumentError(
    "usage: package.jl COLLECTION VERSION REASON GIT_COMMIT [FORCE]",
))

force = length(ARGS) == 5 ? parse(Bool, ARGS[5]) : false
package = package_collection(
    Symbol(ARGS[1]),
    VersionNumber(ARGS[2]);
    reason = ARGS[3],
    git_commit = ARGS[4],
    force
)

println("collection=$(package.collection)")
println("version=$(package.version)")
println("tag=$(package.tag)")
println("archive=$(abspath(package.archive))")
println("archive_sha256=$(package.archive_sha256)")
println("tree_hash=$(package.tree_hash)")
println("package=$(abspath(package.package_path))")
