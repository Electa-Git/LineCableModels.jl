# Deterministic value encoding for calculation reuse and saved files.
# This file defines one encoding; consumers do not invent their own hash rules.

function _write_value(io::IO, value)
    if value === nothing
        print(io, "null;")
    elseif value isa Bool
        print(io, value ? "true;" : "false;")
    elseif value isa Integer
        print(io, "i", typeof(value), ':', value, ';')
    elseif value isa AbstractFloat
        print(io, "f", typeof(value), ':', repr(value), ';')
    elseif value isa Complex
        print(io, "c", typeof(value), '(')
        _write_value(io, real(value))
        _write_value(io, imag(value))
        print(io, ");")
    elseif value isa AbstractString
        print(io, "s", ncodeunits(value), ':', value, ';')
    elseif value isa Symbol
        _write_value(io, String(value))
    elseif value isa AbstractVector || value isa Tuple
        print(io, "a", length(value), '[')
        foreach(item -> _write_value(io, item), value)
        print(io, "];")
    elseif value isa AbstractArray
        # Retain shape for matrix/tensor products while preserving every
        # existing vector/scalar fingerprint byte for historical artifacts.
        print(io, "array{")
        _write_value(io, size(value))
        _write_value(io, vec(value))
        print(io, "};")
    elseif value isa NamedTuple
        print(io, "n", length(value), '{')
        for (name, item) in pairs(value)
            _write_value(io, String(name))
            _write_value(io, item)
        end
        print(io, "};")
    elseif value isa AbstractDict
        ordered = sort!(collect(keys(value)); by = string)
        print(io, "d", length(ordered), '{')
        for key in ordered
            _write_value(io, string(key))
            _write_value(io, value[key])
        end
        print(io, "};")
    else
        throw(ArgumentError(
            "semantic calculation records cannot encode $(typeof(value)); lower it to scalar records first",
        ))
    end
    return io
end

function semantic_sha256(value)
    io = IOBuffer()
    _write_value(io, value)
    return bytes2hex(sha256(take!(io)))
end

"""
    semantic_sha256(result, coordinates)

Hash the numerical observations and recorded terminal/parameter coordinates of
a completed calculation. Storage, reading and reporting use the same methods.
"""
function semantic_sha256(result::AbstractCoreResult, coordinates::NamedTuple)
    return semantic_sha256((Z=vec(observe(result, Z)), Y=vec(observe(result, Y)),
        frequencies=frequencies(result), port_order=coordinates.port_order, basis=basis(result)))
end

function semantic_sha256(result::MomentResult, coordinates::NamedTuple)
    return semantic_sha256(NamedTuple(result))
end

function semantic_sha256(result::AbstractParametricResult, coordinates::NamedTuple)
    return semantic_sha256((points=[semantic_sha256(point, coordinates) for point in result],
        axes=coordinates.axes))
end
