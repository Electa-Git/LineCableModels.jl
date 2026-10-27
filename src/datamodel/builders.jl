"""
    @cable id begin
        parts...
    end

Build the cable `id` from the parts of the block, from the inside out, as
[`build`](@ref)`(CableDesign, id, parts...)` does.

# Arguments

- `id`: name of the cable.
- `parts`: declared parts from the inside out, each on its own line of the block.

# Returns

- A [`CableDesign`](@ref).

# Examples

```jldoctest
copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
polyethylene = Material(:insulator, 1.0e14, 2.3)
wire = @cable "insulated_wire" begin
    @terminal :core begin
        core(copper; r = 4.0e-3) # [m]
    end
    insulation(polyethylene; t = 1.0e-3) # [m]
end
@assert r_ex(wire.outer) ≈ 5.0e-3
# output
```
"""
macro cable(id, block)
    parts = [esc(part) for part in block.args if !(part isa LineNumberNode)]
    return :($build($CableDesign, $(esc(id)), $(parts...)))
end

"""
    @terminal name begin
        parts...
    end

Assign the conductive regions of the parts of the block to the terminal `name`, as
[`Terminal`](@ref)`(name, Stack(parts...))` does.

# Arguments

- `name`: name of the terminal, a `Symbol`.
- `parts`: declared parts from the inside out, each on its own line of the block.

# Returns

- A [`Terminal`](@ref).
"""
macro terminal(name, block)
    parts = [esc(part) for part in block.args if !(part isa LineNumberNode)]
    return :($Terminal($(esc(name)), $Stack($(parts...))))
end

"""
$(TYPEDSIGNATURES)

Declare a solid round conductor of `material` at the center of the cable.

# Keywords

- `r`: radius \\[m\\].
- `tag`: name of the region. Default: `:core`.

# Returns

- A [`Region`](@ref) with a [`Circle`](@ref) of radius `r`.
"""
core(material::Material; r, tag = :core) = Region(tag, Circle(r), material)

"""
$(TYPEDSIGNATURES)

Declare an insulation layer of `material` around the part inside it.

# Keywords

- `t`: thickness \\[m\\].
- `tag`: name of the region. Default: `:insulation`.

# Returns

- A [`Region`](@ref) at the [`Offset`](@ref) `t` from the outline of the part inside it.
"""
insulation(material::Material; t, tag = :insulation) = Region(tag, Offset(t), material)

"""
$(TYPEDSIGNATURES)

Declare a screen layer of `material` around the part inside it, such as a semiconducting
screen. It takes the keywords of [`insulation`](@ref), with the default tag `:screen`.
"""
screen(material::Material; t, tag = :screen) = Region(tag, Offset(t), material)

"""
$(TYPEDSIGNATURES)

Declare a sheath layer of `material` around the part inside it. It takes the keywords of
[`insulation`](@ref), with the default tag `:sheath`.
"""
sheath(material::Material; t, tag = :sheath) = Region(tag, Offset(t), material)

"""
$(TYPEDSIGNATURES)

Declare a bedding layer of `material` around the part inside it. It takes the keywords of
[`insulation`](@ref), with the default tag `:bedding`.
"""
bedding(material::Material; t, tag = :bedding) = Region(tag, Offset(t), material)

"""
$(TYPEDSIGNATURES)

Declare a jacket layer of `material` around the part inside it. It takes the keywords of
[`insulation`](@ref), with the default tag `:jacket`.
"""
jacket(material::Material; t, tag = :jacket) = Region(tag, Offset(t), material)

"""
$(TYPEDSIGNATURES)

Declare a stranded conductor of `material` at the center of the cable. Courses of round wires
or rectangular strands lie around a center wire, inside the outline that `boundary`
declares. Each course lies on the course inside it, at no offset from its outline. A
circular pattern bends a rectangular strand into an annular sector of the same area on its
course. A fill takes `boundary` without the wires and sectors, and its region takes the name
`<tag>_fill`.

# Keywords

- `center`: cross-section of the center wire.
- `shape`: cross-section of one strand, a circle or a rectangle.
- `n`: number of members in each course, from the inside out.
- `lay`: [`LayRatio`](@ref) of every course.
- `boundary`: outline of the conductor, a [`Circle`](@ref) that the fill takes.
- `fill`: material of the fill. Default: an insulator with infinite resistivity and unit
  relative permittivity and permeability.
- `tag`: name of the regions of the center wire and the courses. Default: `:stranded`.

# Returns

- An [`Enclosure`](@ref) of the fill region around a [`Stack`](@ref) of the center wire and
  the [`Lay`](@ref) of each course.

# Errors

- Throws `DomainError` when a course does not fit on the course inside it, or when a wire or
  sector crosses `boundary`.
"""
function stranded(material::Material; center, shape,
        n::Union{Tuple{Vararg{Integer}}, AbstractVector{<:Integer}}, lay, boundary,
        fill = Material(:insulator, Inf), tag = :stranded)
    courses = [Lay(Helix(lay),
                   Repeat(CircularPattern(count; r = Offset(0)), Region(tag, shape, material)))
               for count in n]
    return Enclosure(Region(Symbol(tag, "_fill"), boundary, fill),
        Stack(Region(tag, center, material), courses...))
end

"""
$(TYPEDSIGNATURES)

Declare a course of round armor wires of `material` laid on the part inside it. The wires
lie at no offset from the outline of that part. A fill completes the layer that they occupy,
as thick as their diameter, and its region takes the name `<tag>_fill`.

# Keywords

- `shape`: cross-section of one wire, a [`Circle`](@ref).
- `n`: number of wires.
- `lay`: [`LayRatio`](@ref) of the course.
- `fill`: material of the fill. Default: an insulator with infinite resistivity and unit
  relative permittivity and permeability.
- `tag`: name of the regions of the wires. Default: `:armor`.

# Returns

- An [`Enclosure`](@ref) of the fill region around the [`Lay`](@ref) of the wires.

# Errors

- Throws `DomainError` when the wires do not fit on their course.
"""
function armor(material::Material; shape::Circle, n::Integer, lay,
        fill = Material(:insulator, Inf), tag = :armor)
    wires = Lay(Helix(lay),
        Repeat(CircularPattern(n; r = Offset(0)), Region(tag, shape, material)))
    return Enclosure(Region(Symbol(tag, "_fill"), Offset(2shape.r), fill), wires)
end
