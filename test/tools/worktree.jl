# A git revision checked out in a temporary worktree, for the tools that compare it with
# the working tree (`performance.jl`, `equivalence.jl`). Include this file after
# defining `REPOSITORY`.

git(arguments...) = Cmd(["git", "-C", REPOSITORY, arguments...])
julia(project, program, arguments...) =
    `$(Base.julia_cmd()) --startup-file=no --project=$project -e $program $arguments`

# The standard output of `command`, which must succeed.
function output(command)
    out, err = IOBuffer(), IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = out, stderr = err))
    success(process) || error("Command failed: $(command)\n$(String(take!(err)))")
    return String(take!(out))
end

const VERSIONS = raw"""
using Pkg
for (_, p) in Pkg.dependencies()
    p.version === nothing || println(p.name, "\t", p.version)
end
"""
versions(project) = Dict(split(line, '\t') for line in eachline(IOBuffer(output(julia(project, VERSIONS)))))

# Runs `f(project)` with the test project of `reference` in a temporary worktree. Its
# environment starts from this repository's `Manifest.toml`, resolved again for
# `reference`. The tool prints the dependency versions that then differ from the working
# tree. The worktree is removed afterwards.
function at_revision(f, reference)
    success(pipeline(git("rev-parse", "--verify", "--quiet", reference * "^{commit}");
        stdout = devnull)) || error("Unknown git revision: $reference")
    directory = mktempdir()
    run(pipeline(git("worktree", "add", "--detach", "--quiet", directory, reference);
        stderr = devnull))
    try
        manifest = joinpath(REPOSITORY, "Manifest.toml")
        isfile(manifest) && cp(manifest, joinpath(directory, "Manifest.toml"))
        project = joinpath(directory, "test")
        println("Preparing $reference in a temporary worktree...")
        output(julia(project, "using Pkg; Pkg.resolve(); Pkg.instantiate()"))
        before, after = versions(project), versions(joinpath(REPOSITORY, "test"))
        differing = sort!([name for name in keys(after)
            if haskey(before, name) && before[name] != after[name]])
        println(isempty(differing) ? "Dependency versions: identical on both sides." :
            "Dependency versions that differ at $reference: " *
            join(["$name $(before[name]) (working tree $(after[name]))" for name in differing], ", "))
        return f(project)
    finally
        run(pipeline(git("worktree", "remove", "--force", directory); stderr = devnull))
    end
end

# The revision's own copy of the file `relative` when it has one (containing `marker`
# when given), otherwise the working tree's copy. `own` tells which.
function own_copy(directory, relative; marker = nothing)
    path = joinpath(directory, relative)
    own = isfile(path) && (marker === nothing || occursin(marker, read(path, String)))
    return (; path = own ? path : joinpath(REPOSITORY, relative), own)
end
