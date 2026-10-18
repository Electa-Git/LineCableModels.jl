"""
$(TYPEDEF)

Store terminal voltages and currents over time. Array axes are end × terminal × time,
with sending before receiving. Currents enter the line at each end.

$(TYPEDFIELDS)
"""
struct TransientResponse{T, D} <: AbstractCoreResult
    "Time axis \\[s\\]."
    t::Vector{T}
    "Terminal voltages \\[V\\]."
    voltage::Array{T, 3}
    "Terminal currents entering the line \\[A\\]."
    current::Array{T, 3}
    "Terminal labels in physical order."
    coordinates::Vector{String}
    "Supplemental computation record."
    details::D
end

"""
Select terminal voltage \\[V\\] from a transient response.
"""
function voltage end

"""
Select terminal current entering the line \\[A\\] from a transient response.
"""
function current end

details(response::TransientResponse) = response.details

function retain_gridpoint(source::TransientResponse, id; fields=(;))
    retained = ComputationDetails(merge(details(source).data,fields,(gridpoint=id,)))
    return TransientResponse(source.t,source.voltage,source.current,source.coordinates,retained)
end

function completed_inputs(problem::TransientAnalysisProblem)
    return (line=completed_inputs(problem.line),sending=problem.sending,
        receiving=problem.receiving,window=problem.window)
end
