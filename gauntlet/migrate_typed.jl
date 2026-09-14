# Finite export bridge for the pre-extraction checkpoint. Run in that checkout's
# test/gauntlet environment. This file never launches a numerical computation.
using LineCableModels, JLD2, SHA, TOML

length(ARGS)==3 || error("usage: migrate_typed.jl INPUT.jld2 OUTPUT.jld2 DECLARATION.toml")
const LEGACY_CHECKPOINT="2d694a2444f52942e9d14cd5ff835760295ff532"
repository=pkgdir(LineCableModels)
readchomp(`git -C $repository rev-parse HEAD`)==LEGACY_CHECKPOINT ||
    error("typed migration requires the pinned pre-extraction checkpoint $LEGACY_CHECKPOINT")
success(`git -C $repository diff --quiet HEAD -- src ext test/gauntlet Project.toml`) ||
    error("the legacy implementation must match its pinned checkpoint")
input, output, declaration_path=abspath.(ARGS)
ispath(output) && error("migration destination already exists: $output")
declaration=TOML.parsefile(declaration_path)
ports=String.(declaration["port_order"])
allunique(ports) || error("terminal labels must be unique and in stored matrix order")
# The legacy module supplies the original types; it is not installed as an alias
# in the new package. Published input files are never modified.
include(joinpath(repository, "test", "gauntlet", "runner.jl"))
document=JLD2.load(input)
result=document[get(declaration, "result_key", "result")]
result isa LineParameters ||
    error("the selected payload is not a legacy LineParameters result")
length(ports)==size(result.Z, 1) || error("terminal order does not span the matrices")
LineCableModels.domain(result)===PhaseDomain ||
    error("this finite bridge handles phase-domain matrices")
problem=document[get(declaration, "problem_key", "problem")]
digest=bytes2hex(open(sha256, input))
files=NamedTuple[]
mkpath(dirname(output))
legacy_path=joinpath(dirname(output), "files", "legacy", basename(input))
ispath(legacy_path) && error("legacy evidence destination is occupied")
mkpath(dirname(legacy_path));
cp(input, legacy_path)
push!(files, (path = relpath(legacy_path, dirname(output)), sha256 = digest))
sources=[(path, sha256 = bytes2hex(open(sha256, joinpath(repository, path))),
             source = read(joinpath(repository, path)))
         for path in
             split(readchomp(`git -C $repository ls-files src ext test/gauntlet Project.toml`), '\n')
         if endswith(path, ".jl") || endswith(path, ".toml")]
temporary=tempname(dirname(output))
try
    JLD2.jldsave(
        temporary; schema_version = 2, status = :complete, kind = :gauntlet_calculation,
        case_id = declaration["case_id"], backend = declaration["backend"],
        selection = (id = :legacy_import,), formulation = (
            legacy_checkpoint = LEGACY_CHECKPOINT,
            description = get(declaration, "description", "Retained legacy calculation; see original payload")),
        problem = LineCableModels.ImportExport.serialize_value(problem),
        frequencies = copy(result.f), basis = LineCableModels.basis(result), domain = :PhaseDomain,
        Z = copy(result.Z.values), Y = copy(result.Y.values), port_order = ports,
        comparison_unsupported = get(details(result), :comparison_unsupported, (;)),
        retained_files = files, implementation = LEGACY_CHECKPOINT, source_evidence = sources,
        computation_details = (legacy = (
            source_sha256 = digest, checkpoint = LEGACY_CHECKPOINT,
            declaration = TOML.parsefile(declaration_path)),))
    mv(temporary, output)
    write(output*".sha256", bytes2hex(open(sha256, output))*"  "*basename(output)*"\n")
finally
    isfile(temporary) && rm(temporary)
end
println("Converted ", output, "; original SHA-256 ", digest)
