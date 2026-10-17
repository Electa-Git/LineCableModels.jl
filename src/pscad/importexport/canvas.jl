# Layout constants: PowerImpedance.jl `src/PSCAD/plan.jl` lines 46-61 and 320-347.
# Interface ports: PowerImpedance.jl `src/PSCAD/emitters/cables.jl` lines 125-137, (±18, 54(k - 1)).
const _PSCAD_GRID = 18
const _PSCAD_INTERFACE_PORT_SPACING = 54

"""
    PSCADCanvas(schematic)

Own component identifiers and the next layout cell of an XML schematic.
Use one canvas per schematic. Identifiers start above those in the document.
Cells measure 432 by 216 pixels, six per row, starting at `(72, 72)`.
"""
mutable struct PSCADCanvas{N}
    schematic::N
    next_id::Int
    next_cell::Int
end

function PSCADCanvas(schematic)
    nodes = findall("//*[@id]", schematic)
    identifier = maximum(node -> parse(Int, node["id"]), nodes; init = 0) + 1
    return PSCADCanvas{typeof(schematic)}(schematic, identifier, 0)
end

"""
    allocate_cell!(canvas)

Reserve a cell and return its upper-left pixel coordinates on the PSCAD grid.
"""
function allocate_cell!(canvas::PSCADCanvas)
    index = canvas.next_cell
    canvas.next_cell += 1
    column, row = index % 6, index ÷ 6
    return (_PSCAD_GRID * (4 + 24column), _PSCAD_GRID * (4 + 12row))
end

function validate(position::NTuple{2, Int}, ::Type{PSCADCanvas})
    all(coordinate -> coordinate % _PSCAD_GRID == 0, position) ||
        throw(ArgumentError("PSCAD positions must be multiples of 18 pixels"))
    return position
end

"""
    component!(canvas, definition, position, parameters=Pair{String,String}[])

Place a component at pixel coordinates `position`, with string-valued parameters.
Unqualified definitions belong to `master`. Qualified bindings retain their prefix.
Return the component XML node.
"""
function component!(canvas::PSCADCanvas, definition::AbstractString,
        position::NTuple{2, Int}, parameters = Pair{String, String}[])
    validate(position, PSCADCanvas)
    binding = occursin(':', definition) ? definition : "master:$definition"
    component = _pscad_user!(canvas.schematic, canvas.next_id, binding, parameters;
        x = position[1], y = position[2])
    canvas.next_id += 1
    return component
end

"""
    node_label!(canvas, name, position)

Place an electrical node label at pixel coordinates `position`.
Labels with the same name connect electrical nodes within the schematic.
"""
function node_label!(canvas::PSCADCanvas, name::AbstractString, position::NTuple{2, Int})
    return component!(canvas, "master:nodelabel", position, ["Name" => String(name)])
end

"""
    wire!(canvas, a, b)

Join distinct, axis-aligned pixel positions `a` and `b` on the PSCAD grid.
Return a wire XML node with vertices relative to its first endpoint.
"""
function wire!(canvas::PSCADCanvas, a::NTuple{2, Int}, b::NTuple{2, Int})
    validate((a, b), PSCADCanvas)
    wire = addelement!(canvas.schematic, "Wire")
    # Native XML shape: mhi-xml 1.3.0, `pscad/component.py` lines 471-484.
    _pscad_attributes!(wire, (id = canvas.next_id, name = "", classid = "WireOrthogonal",
        x = a[1], y = a[2], w = abs(b[1] - a[1]), h = abs(b[2] - a[2]),
        orient = 0, disable = "false"))
    for (x, y) in ((0, 0), (b[1] - a[1], b[2] - a[2]))
        vertex = addelement!(wire, "vertex")
        _pscad_attributes!(vertex, (; x, y))
    end
    canvas.next_id += 1
    return wire
end

function validate(vertices::Tuple{NTuple{2, Int}, NTuple{2, Int}}, ::Type{PSCADCanvas})
    a, b = vertices
    validate(a, PSCADCanvas)
    validate(b, PSCADCanvas)
    a != b && (a[1] == b[1] || a[2] == b[2]) ||
        throw(ArgumentError("PSCAD wires must be nonzero and axis-aligned"))
    return vertices
end
