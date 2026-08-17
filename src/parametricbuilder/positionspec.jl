struct PositionSpec{Kind, P <: Tuple, C <: Tuple}
    parameters::P
    connections::C
end

_position_kind(::PositionSpec{Kind}) where {Kind} = Kind

function PositionSpec(::Val{Kind}, parameters::P, connections::C) where {
        Kind, P <: Tuple, C <: Tuple
}
    all(value -> value isa Real, parameters) || throw(ArgumentError(
        "position coordinates and spacing must resolve to real numbers",
    ))
    return PositionSpec{Kind, P, C}(parameters, connections)
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
        key = String(key_raw)
        if value isa Integer
            for map in maps
                map[key] = Int(value)
            end
        elseif value isa Tuple || value isa AbstractVector
            length(value) == count || throw(DimensionMismatch(
                "phase '$key' requires $count assignments; got $(length(value))",
            ))
            for index in 1:count
                value[index] isa Integer || throw(ArgumentError(
                    "phase assignments must be integers",
                ))
                maps[index][key] = Int(value[index])
            end
        else
            throw(ArgumentError(
                "phase '$key' must be an integer or a collection of $count integers",
            ))
        end
    end
    return Tuple(maps)
end

_point_builder(x, y, connections) = PositionSpec(Val(:point), (x, y), connections)
function _trifoil_builder(x, y, spacing, connections)
    PositionSpec(Val(:trifoil), (x, y, spacing), connections)
end
function _hflat_builder(x, y, spacing, connections)
    PositionSpec(Val(:hflat), (x, y, spacing), connections)
end
function _vflat_builder(x, y, spacing, connections)
    PositionSpec(Val(:vflat), (x, y, spacing), connections)
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

- A [`Gridspace`](@ref) of position declarations. The parent system builder
  converts each resolved declaration into a cable position.
"""
function at(; x, y, phases = nothing, combine::Symbol = :product)
    connections = _phase_maps(phases, 1)
    return Gridspace{PositionSpec}(
        _point_builder,
        (_gridspace_axis(x), _gridspace_axis(y), _gridspace_axis(connections)),
        (:x, :y, :phases);
        combine
    )
end

"""
$(TYPEDSIGNATURES)

Describe three identical cables arranged in equilateral trifoil formation.
The formation is resolved when the parent [`SystemBuilder`](@ref) materializes
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

- A [`Gridspace`](@ref) of trifoil declarations. Each resolved declaration is
  one object-valued input to its parent system specification.

# Errors

- System materialization throws `ArgumentError` when `spacing` is nonpositive
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
    connections = _phase_maps(phases, 3)
    return Gridspace{PositionSpec}(
        _trifoil_builder,
        (_gridspace_axis(x), _gridspace_axis(y),
            _gridspace_axis(spacing), _gridspace_axis(connections)),
        (:x, :y, :spacing, :phases);
        combine
    )
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

- A [`Gridspace`](@ref) of horizontal flat-formation declarations.

# Errors

- Throws `ArgumentError` when `n` is nonpositive. System materialization also
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
    connections = _phase_maps(phases, n)
    return Gridspace{PositionSpec}(
        _hflat_builder,
        (_gridspace_axis(x), _gridspace_axis(y),
            _gridspace_axis(spacing), _gridspace_axis(connections)),
        (:x, :y, :spacing, :phases);
        combine
    )
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

- A [`Gridspace`](@ref) of vertical flat-formation declarations.

# Errors

- Throws `ArgumentError` when `n` is nonpositive. System materialization also
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
    connections = _phase_maps(phases, n)
    return Gridspace{PositionSpec}(
        _vflat_builder,
        (_gridspace_axis(x), _gridspace_axis(y),
            _gridspace_axis(spacing), _gridspace_axis(connections)),
        (:x, :y, :spacing, :phases);
        combine
    )
end
