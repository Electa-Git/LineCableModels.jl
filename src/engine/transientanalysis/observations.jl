_response_index(index::Integer, coordinates) = index
function _response_index(label::AbstractString, coordinates)
    index = findfirst(==(label), coordinates)
    index === nothing && throw(ArgumentError("unknown coordinate label $label"))
    return index
end

"""
$(TYPEDSIGNATURES)

Return terminal voltage \\[V\\] over the full time axis. Select the end by index
or label (`"sending"`, `"receiving"`) and the terminal by index or its coordinate label.
"""
function observe(response::TransientResponse, ::typeof(voltage), line_end, terminal, ::Colon)
    return response.voltage[_response_index(line_end, ("sending", "receiving")),
        _response_index(terminal, response.coordinates), :]
end

"""
$(TYPEDSIGNATURES)

Return terminal current entering the line \\[A\\] over the full time axis. Select the
end and terminal by index or label, as for voltage.
"""
function observe(response::TransientResponse, ::typeof(current), line_end, terminal, ::Colon)
    return response.current[_response_index(line_end, ("sending", "receiving")),
        _response_index(terminal, response.coordinates), :]
end

"""Return the terminal voltages \\[V\\] of both ends, indexed by end, terminal and time."""
observe(response::TransientResponse, ::typeof(voltage)) = response.voltage

"""Return the terminal currents entering the line \\[A\\], indexed by end, terminal and time."""
observe(response::TransientResponse, ::typeof(current)) = response.current

Units.quantity(::typeof(voltage)) = Units.Quantity{:terminal_voltage}()
Units.quantity(::typeof(current)) = Units.Quantity{:terminal_current}()

Commons.observables(::Type{<:TransientResponse}) = (voltage, current)
basis(::TransientResponse) = :total

function Commons.observation_quantity(response::TransientResponse, request;
        unit=nothing, clip=true, atol=nothing, frequencies=nothing)
    atol===nothing || throw(ArgumentError("transient responses have no resolution cutoffs"))
    frequencies===nothing || throw(ArgumentError("transient responses use their stored time axis"))
    resolved=Commons.observation_request(response,request)
    raw=observe(response,resolved.identity)
    indices=isempty(resolved.indices) ? (Colon(),Colon(),Colon()) : resolved.indices
    length(indices)==3 || throw(DimensionMismatch("transient requests require three selectors"))
    selected=map(Commons.observation_indices,indices,size(raw))
    coordinates=(kind=:matrix,indices,rows=selected[1],columns=selected[2],
        samples=selected[3],times=copy(response.t[selected[3]]),
        time_unit=Units.units(:base,:second),extent=size(response.voltage),
        labels=["sending","receiving"],domain=:LineEnd,
        column_labels=copy(response.coordinates),column_domain=:PhaseDomain)
    q=resolved.quantity
    native=Units.native_unit(q,basis(response))
    target=unit===nothing ? Units.display_unit(q,basis(response)) : unit
    return (request,quantity=q,family=:transient,statistic=:value,
        values=Commons.detach(raw[indices...],Units.scale_factor(native,target)),unit=target,
        basis=basis(response),coordinates,thresholds=nothing,available=nothing,
        engineering_zero=nothing,clipped=false,missing_reason=nothing)
end
