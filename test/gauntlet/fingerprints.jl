# Canonical value encoding shared by cache provenance and artifact transport.
# This file defines one encoding; consumers do not invent their own hash rules.

function _canonical_write(io::IO, value)
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
        _canonical_write(io, real(value))
        _canonical_write(io, imag(value))
        print(io, ");")
    elseif value isa AbstractString
        print(io, "s", ncodeunits(value), ':', value, ';')
    elseif value isa Symbol
        _canonical_write(io, String(value))
    elseif value isa AbstractVector || value isa Tuple
        print(io, "a", length(value), '[')
        foreach(item -> _canonical_write(io, item), value)
        print(io, "];")
    elseif value isa AbstractArray
        # Retain shape for matrix/tensor products while preserving every
        # existing vector/scalar fingerprint byte for historical artifacts.
        print(io, "array{")
        _canonical_write(io, size(value))
        _canonical_write(io, vec(value))
        print(io, "};")
    elseif value isa NamedTuple
        print(io, "n", length(value), '{')
        for (name, item) in pairs(value)
            _canonical_write(io, String(name))
            _canonical_write(io, item)
        end
        print(io, "};")
    elseif value isa AbstractDict
        ordered = sort!(collect(keys(value)); by = string)
        print(io, "d", length(ordered), '{')
        for key in ordered
            _canonical_write(io, string(key))
            _canonical_write(io, value[key])
        end
        print(io, "};")
    else
        throw(ArgumentError(
            "semantic provenance cannot encode $(typeof(value)); lower it to scalar records first",
        ))
    end
    return io
end

function semantic_sha256(value)
    io = IOBuffer()
    _canonical_write(io, value)
    return bytes2hex(sha256(take!(io)))
end
