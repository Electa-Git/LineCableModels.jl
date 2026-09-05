include("artifacts.jl")
using .GauntletArtifacts

length(ARGS) in (4, 5) || throw(ArgumentError(
    "usage: package.jl COLLECTION VERSION REASON GIT_COMMIT [FORCE]",
))

repository = normpath(joinpath(@__DIR__, "..", ".."))
head = lowercase(readchomp(`git -C $repository rev-parse HEAD`))
requested_commit = lowercase(strip(ARGS[4]))
requested_commit == head || throw(ArgumentError(
    "official Gauntlet packages must record the current full Git commit $head",
))
isempty(readchomp(`git -C $repository status --porcelain`)) || throw(ArgumentError(
    "official Gauntlet packages require a clean repository worktree",
))

force = length(ARGS) == 5 ? parse(Bool, ARGS[5]) : false
package = package_collection(
    Symbol(ARGS[1]),
    VersionNumber(ARGS[2]);
    reason = ARGS[3],
    git_commit = requested_commit,
    force
)

println("collection=$(package.collection)")
println("version=$(package.version)")
println("tag=$(package.tag)")
println("archive=$(abspath(package.archive))")
println("archive_sha256=$(package.archive_sha256)")
println("tree_hash=$(package.tree_hash)")
println("package=$(abspath(package.package_path))")
