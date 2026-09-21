# Shared input reading for the two manual inspectors. No solver, comparison,
# analysis writer, or compatibility conversion runs in this function.
using TOML, JLD2

function read_inspection_operands(campaign_directory, benchmark_id; previous=false)
    directory=joinpath(campaign_directory,string(benchmark_id))
    state=TOML.parsefile(joinpath(directory,"state.toml"))
    haskey(state,"current") || throw(ArgumentError("$benchmark_id has no completed saved attempt; a computation is required"))
    state["state"]=="complete" || previous || throw(ArgumentError(
        "$benchmark_id is $(state["state"]); set use_previous_complete=true explicitly to inspect its previous completed attempt"))
    attempt=joinpath(directory,state["current"])
    operands=map((:reference,:candidate)) do role
        path=joinpath(attempt,string(role),"calculation.jld2")
        kind=jldopen(file -> file["kind"],path,"r")
        kind===:gauntlet_moments && throw(ArgumentError(
            "$benchmark_id / $role retains only the removed gauntlet_moments representation. " *
            "The current UQ owner requires a retained scientific result; regenerate this UQ calculation. " *
            "Do not reconstruct correlated uncertainty from marginal means/stds."))
        Gauntlet.read_calculation(path;evidence=:numerical)
    end
    # The existing Gauntlet method attaches saved timing/performance measurements.
    return Gauntlet.read_benchmark((id=benchmark_id,reference=operands[1],
        candidate=operands[2],analyses=Dict{String,Any}[]))
end
