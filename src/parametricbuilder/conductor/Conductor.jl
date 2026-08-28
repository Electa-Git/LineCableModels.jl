"""
    LineCableModels.ParametricBuilder.Conductor

Qualified conductor constructions for extension code.

These methods lower through the same `core`, `sheath`, `wires`, `tape`, and
`Group` functions used by the ordinary construction surface.
"""
module Conductor

import ..ParametricBuilder: Group, core, sheath, tape, terminal, wires
import ...DataModel

"""
    Conductor.Solid(tag, material; r, combine=:product)

Declare one solid circular conductor region of radius `r` [m].
"""
function Solid(tag, material; r, combine::Symbol = :product)
    return core(material; r, tag, combine)
end

"""
    Conductor.Shell(tag, material; t, combine=:product)

Declare one outward conductor shell of thickness `t` [m].
"""
function Shell(tag, material; t, combine::Symbol = :product)
    return sheath(material; t, tag, combine)
end

"""
    Conductor.Wires(tag, wire, material; n, r, lay=LayRatio(11), compact=nothing)

Declare one repeated conductor course on a ring of member-centre radius `r`
[m]. The result is the same `Group` produced by [`wires`](@ref).
"""
function Wires(
        tag,
        wire,
        material;
        n,
        r,
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    return wires(
        material;
        wire,
        n,
        r,
        lay,
        compact,
        tag,
        combine
    )
end

"""
    Conductor.Strip(tag, section, material; lay=LayRatio(11), compact=nothing)

Declare one conductive tape section and coalesce it into terminal `tag`.
"""
function Strip(
        tag,
        section,
        material;
        lay = DataModel.LayRatio(11),
        compact = nothing,
        combine::Symbol = :product
)
    part = tape(
        material;
        section,
        n = 1,
        lay,
        compact,
        tag,
        combine
    )
    return terminal(tag, part; combine)
end

"""
    Conductor.Tubular(tag, annulus, material; path=nothing, combine=:product)

Declare an annular conductor and retain an optional longitudinal path.
"""
function Tubular(
        tag,
        annulus,
        material;
        path = nothing,
        combine::Symbol = :product
)
    region = core(material, annulus; tag, combine)
    return Group(tag, region; path, combine)
end

end
