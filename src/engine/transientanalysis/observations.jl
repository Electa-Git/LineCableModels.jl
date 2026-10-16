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
