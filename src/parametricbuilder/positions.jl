# Position declarations and formation builders.
struct PositionDefinition{Kind, P <: Tuple, C <: Tuple}
    parameters::P
    connections::C
end

function PositionDefinition(::Val{Kind}, parameters::P, connections::C) where {
        Kind, P <: Tuple, C <: Tuple
}
    all(value -> value isa Real, parameters) || throw(ArgumentError(
        "position coordinates and spacing must resolve to real numbers",
    ))
    return PositionDefinition{Kind, P, C}(parameters, connections)
end

const _PhaseEntry = Union{
    Pair{Symbol, <:Any},
    Pair{String, <:Any},
    Tuple{Symbol, <:Any},
    Tuple{String, <:Any}
}

_phase_entries(entry::_PhaseEntry) = (entry,)
_phase_entries(entries::Tuple) = entries
_phase_entries(entries::AbstractVector) = Tuple(entries)
_phase_entries(::Nothing) = ()

function _phase_maps(phases, count::Int)
    maps = [Dict{String, Int}() for _ in 1:count]
    for entry in _phase_entries(phases)
        key_raw, value = entry isa Pair ?
                         (first(entry), last(entry)) : (entry[1], entry[2])
        name = String(key_raw)
        if value isa Integer
            for map in maps
                map[name] = Int(value)
            end
        elseif value isa Tuple || value isa AbstractVector
            length(value) == count || throw(DimensionMismatch(
                "phase '$name' requires $count assignments; got $(length(value))",
            ))
            for index in 1:count
                value[index] isa Integer || throw(ArgumentError(
                    "phase assignments must be integers",
                ))
                maps[index][name] = Int(value[index])
            end
        else
            throw(ArgumentError(
                "phase '$name' must be an integer or a collection of $count integers",
            ))
        end
    end
    return Tuple(maps)
end

_point_builder(x, y, connections) = PositionDefinition(Val(:point), (x, y), connections)
function _trifoil_builder(x, y, spacing, connections)
    PositionDefinition(Val(:trifoil), (x, y, spacing), connections)
end
function _hflat_builder(x, y, spacing, connections)
    PositionDefinition(Val(:hflat), (x, y, spacing), connections)
end
function _vflat_builder(x, y, spacing, connections)
    PositionDefinition(Val(:vflat), (x, y, spacing), connections)
end

"""
$(TYPEDSIGNATURES)

Describe one cable position for use by [`SystemBuilder`](@ref).

# Keywords

- `x`: Horizontal coordinate \\[m\\].
- `y`: Vertical coordinate \\[m\\].
- `phases=nothing`: Phase or conductor assignments. Supply pairs such as
  `:core => 1` and `:sheath => 2`.
- `combine=:product`: Local composition rule for direct varying inputs. Use
  `:zip` to pair compatible axes.

# Returns

- One position declaration for scalar inputs, or a `Gridspace` of position
  declarations when a direct input varies.
"""
function at(; x, y, phases = nothing, combine::Symbol = :product)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    connections = _phase_maps(phases, 1)
    values = (x, y, connections)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PositionDefinition}(_point_builder, grids; combine)
    end
    return _point_builder(values...)
end

"""
$(TYPEDSIGNATURES)

Describe three identical cables arranged in equilateral trifoil formation.
The formation is resolved when the parent [`SystemBuilder`](@ref) materialises
the cable system, so its spacing is checked against the cable outer diameter.

# Keywords

- `x=0.0`: Horizontal coordinate of the formation centre \\[m\\].
- `y`: Vertical coordinate of the formation centre \\[m\\].
- `spacing`: Centre-to-centre distance between adjacent cables \\[m\\].
- `phases`: Phase or conductor assignments. Each assignment may contain three
  integers, one for each cable, for example `:core => (1, 2, 3)`.
- `combine=:product`: Local composition rule for direct varying inputs. Use
  `:zip` to pair compatible axes.

# Returns

- One trifoil declaration for scalar inputs, or a `Gridspace` of declarations
  when a direct input varies.

# Errors

- System materialisation throws `ArgumentError` when `spacing` is nonpositive
  or smaller than the cable outer diameter.
- Throws `DimensionMismatch` when a phase assignment does not provide three
  values.
"""
function trifoil(;
        x = 0.0,
        y,
        spacing,
        phases,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    connections = _phase_maps(phases, 3)
    values = (x, y, spacing, connections)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PositionDefinition}(_trifoil_builder, grids; combine)
    end
    return _trifoil_builder(values...)
end

"""
$(TYPEDSIGNATURES)

Describe `n` identical cables in a horizontal flat formation. The first cable
is placed at `(x, y)` and subsequent cables are placed at increasing horizontal
coordinates.

# Keywords

- `x=0.0`: Horizontal coordinate of the first cable \\[m\\].
- `y=0.0`: Vertical coordinate shared by all cables \\[m\\].
- `spacing`: Centre-to-centre distance between adjacent cables \\[m\\].
- `phases`: Phase or conductor assignments. An assignment may provide one
  integer for every cable.
- `n=3`: Number of cables.
- `combine=:product`: Local composition rule for direct varying inputs. Use
  `:zip` to pair compatible axes.

# Returns

- One horizontal formation declaration for scalar inputs, or a `Gridspace` of
  declarations when a direct input varies.

# Errors

- Throws `ArgumentError` when `n` is nonpositive. System materialisation also
  rejects nonpositive spacing and overlapping cables.
- Throws `DimensionMismatch` when a phase assignment does not provide `n`
  values.
"""
function hflat(;
        x = 0.0,
        y = 0.0,
        spacing,
        phases,
        n::Int = 3,
        combine::Symbol = :product
)
    n > 0 || throw(ArgumentError("n must be positive"))
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    connections = _phase_maps(phases, n)
    values = (x, y, spacing, connections)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PositionDefinition}(_hflat_builder, grids; combine)
    end
    return _hflat_builder(values...)
end

"""
$(TYPEDSIGNATURES)

Describe `n` identical cables in a vertical flat formation. The first cable is
placed at `(x, y)` and subsequent cables are placed at decreasing vertical
coordinates.

# Keywords

- `x=0.0`: Horizontal coordinate shared by all cables \\[m\\].
- `y=0.0`: Vertical coordinate of the first cable \\[m\\].
- `spacing`: Centre-to-centre distance between adjacent cables \\[m\\].
- `phases`: Phase or conductor assignments. An assignment may provide one
  integer for every cable.
- `n=3`: Number of cables.
- `combine=:product`: Local composition rule for direct varying inputs. Use
  `:zip` to pair compatible axes.

# Returns

- One vertical formation declaration for scalar inputs, or a `Gridspace` of
  declarations when a direct input varies.

# Errors

- Throws `ArgumentError` when `n` is nonpositive. System materialisation also
  rejects nonpositive spacing and overlapping cables.
- Throws `DimensionMismatch` when a phase assignment does not provide `n`
  values.
"""
function vflat(;
        x = 0.0,
        y = 0.0,
        spacing,
        phases,
        n::Int = 3,
        combine::Symbol = :product
)
    n > 0 || throw(ArgumentError("n must be positive"))
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    connections = _phase_maps(phases, n)
    values = (x, y, spacing, connections)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PositionDefinition}(_vflat_builder, grids; combine)
    end
    return _vflat_builder(values...)
end
