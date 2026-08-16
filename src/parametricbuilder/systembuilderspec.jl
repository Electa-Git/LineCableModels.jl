@gridspace @relax struct EarthParameters{T<:Real}
    rho::T
    eps_r::T=1.0
    mu_r::T=1.0
    thickness::T=Inf
end

"""
$(TYPEDSIGNATURES)

Describe earth properties independently of analysis frequencies. The
frequency-dependent materialized `EarthModel` is constructed when the parent
problem is materialized.
"""
function Earth(;
    rho,
    eps_r=1.0,
    mu_r=1.0,
    thickness=Inf,
    combine::Symbol=:product,
)
    return EarthParameters(; rho, eps_r, mu_r, thickness, combine)
end

function _position_coordinates(
    position::PositionBuilder,
    design::DataModel.CableDesign,
)
    kind = _position_kind(position)
    if kind === :point
        x, y = position.parameters
        return ((x, y, only(position.connections)),)
    end

    x, y, spacing = position.parameters
    spacing > zero(spacing) ||
        throw(ArgumentError("formation spacing must be positive"))
    minimum_spacing = 2 * DataModel.get_outer_radius(design)
    spacing >= minimum_spacing || throw(ArgumentError(
        "formation spacing $spacing is smaller than the minimum non-overlap distance $minimum_spacing",
    ))

    coordinates = if kind === :trifoil
        length(position.connections) == 3 ||
            throw(DimensionMismatch("trifoil requires three phase maps"))
        values = DataModel.trifoil_formation(x, y, spacing / 2)
        ((values[1], values[2]), (values[3], values[4]), (values[5], values[6]))
    elseif kind === :hflat
        ntuple(
            index -> (x + (index - 1) * spacing, y),
            length(position.connections),
        )
    elseif kind === :vflat
        ntuple(
            index -> (x, y - (index - 1) * spacing),
            length(position.connections),
        )
    else
        throw(ArgumentError("unsupported position kind :$kind"))
    end
    return ntuple(
        index -> (
            coordinates[index][1],
            coordinates[index][2],
            position.connections[index],
        ),
        length(position.connections),
    )
end

function _materialize_earth(
    earth::EarthParameters,
    frequencies,
)
    return EarthProps.EarthModel(
        collect(frequencies),
        earth.rho,
        earth.eps_r,
        earth.mu_r;
        t=earth.thickness,
    )
end

function _materialize_earth(
    earth::EarthProps.EarthModel,
    ::AbstractVector,
)
    return earth
end

struct SystemMaterializer
    position_count::Int
end

function (materializer::SystemMaterializer)(identifier, design, values...)
    design isa DataModel.CableDesign ||
        throw(ArgumentError("SystemBuilder design must resolve to CableDesign"))
    expected = materializer.position_count + 4
    length(values) == expected || throw(DimensionMismatch(
        "SystemBuilder expected $expected resolved values after the design; got $(length(values))",
    ))

    positions = values[1:materializer.position_count]
    line_length = values[materializer.position_count + 1]
    temperature = values[materializer.position_count + 2]
    earth = values[materializer.position_count + 3]
    frequencies = values[materializer.position_count + 4]

    line_length isa Real && line_length > zero(line_length) ||
        throw(ArgumentError("line length must be a positive real number"))
    temperature isa Real ||
        throw(ArgumentError("temperature must be a real number"))
    frequencies isa AbstractVector ||
        throw(ArgumentError("frequencies must be an ordinary vector"))

    layouts = tuple((
        _position_coordinates(position, design)
        for position in positions
    )...)
    coordinates = tuple(Iterators.flatten(layouts)...)
    isempty(coordinates) &&
        throw(ArgumentError("SystemBuilder requires at least one cable position"))

    x, y, connections = first(coordinates)
    system = DataModel.LineCableSystem(
        String(identifier),
        line_length,
        DataModel.CablePosition(design, x, y, connections),
    )
    for (next_x, next_y, next_connections) in Iterators.drop(coordinates, 1)
        system = add!(system, design, next_x, next_y, next_connections)
    end

    earth_model = _materialize_earth(earth, frequencies)
    return Engine.LineParametersProblem(
        system;
        temperature,
        earth_props=earth_model,
        frequencies=collect(frequencies),
    )
end

struct SystemBuilderSpec{D,P<:Tuple,L,T,E,F,C} <:
       AbstractSpec{Engine.LineParametersProblem}
    identifier::String
    design::D
    positions::P
    line_length::L
    temperature::T
    earth::E
    frequencies::F
    combine::C
end

_flatten_positions(position::Gridspace{PositionBuilder}) = (position,)
function _flatten_positions(positions::Union{Tuple,AbstractVector})
    flattened = ()
    for position in positions
        flattened = (flattened..., _flatten_positions(position)...)
    end
    return flattened
end
_flatten_positions(value) = throw(ArgumentError(
    "SystemBuilder positions must be created by at/trifoil/hflat/vflat; got $(typeof(value))",
))

"""
$(TYPEDSIGNATURES)

Describe a line-cable analysis problem from a cable design, its spatial
arrangement, earth properties, and analysis frequencies. Iterating the
returned specification materializes
[`LineCableModels.Engine.LineParametersProblem`](@ref) objects without running
the numerical analysis.

# Arguments

- `identifier`: Name assigned to the materialized cable system.
- `design`: A materialized cable design or a [`CableBuilder`](@ref)
  specification.
- `positions`: One position declaration, or a tuple or vector of declarations
  created by [`at`](@ref), [`trifoil`](@ref), [`hflat`](@ref), or
  [`vflat`](@ref).

# Keywords

- `length=1000.0`: Line length \\[m\\].
- `temperature=20.0`: Conductor temperature \\[°C\\].
- `earth`: Earth declaration created by [`Earth`](@ref), or a materialized
  earth model.
- `frequencies`: Ordinary vector of analysis frequencies \\[Hz\\], or a `Grid`
  whose values are complete frequency vectors.
- `combine=:product`: Local composition rule for direct varying inputs. Use
  `:zip` to pair compatible axes.

# Returns

- A line-parameter problem specification. Pass it to
  [`LineCableModels.Engine.compute!`](@ref) with a
  [`LineCableModels.Engine.Formulation`](@ref) to evaluate its materialized
  configurations.

# Errors

- Throws `ArgumentError` for unsupported design, position, earth, frequency,
  or composition inputs.
- Materialization reports overlapping positions, invalid line lengths, and
  other invalid physical inputs from the materialized cable-system model.
"""
function SystemBuilder(
    identifier::AbstractString,
    design,
    positions;
    length=1000.0,
    temperature=20.0,
    earth,
    frequencies,
    combine::Symbol=:product,
)
    design isa Union{DataModel.CableDesign,AbstractSpec{DataModel.CableDesign}} ||
        throw(ArgumentError("design must be a materialized CableDesign or CableBuilder spec"))
    earth isa Union{EarthProps.EarthModel,AbstractSpec{EarthParameters}} ||
        throw(ArgumentError("earth must be a materialized EarthModel or Earth spec"))
    frequencies isa Union{AbstractVector,AbstractGrid} || throw(ArgumentError(
        "frequencies must be an ordinary vector or a Grid of complete vectors",
    ))
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    position_tuple = _flatten_positions(positions)
    isempty(position_tuple) &&
        throw(ArgumentError("SystemBuilder requires at least one position"))
    return SystemBuilderSpec(
        String(identifier),
        design,
        position_tuple,
        length,
        temperature,
        earth,
        frequencies,
        Val(combine),
    )
end

function gridspace(spec::SystemBuilderSpec)
    axes = (
        _gridspace_axis(spec.identifier),
        _gridspace_axis(spec.design),
        map(_gridspace_axis, spec.positions)...,
        _gridspace_axis(spec.line_length),
        _gridspace_axis(spec.temperature),
        _gridspace_axis(spec.earth),
        _gridspace_axis(spec.frequencies),
    )
    names = (
        :identifier,
        :design,
        ntuple(index -> Symbol(:position_, index), length(spec.positions))...,
        :length,
        :temperature,
        :earth,
        :frequencies,
    )
    return Gridspace{Engine.LineParametersProblem}(
        SystemMaterializer(length(spec.positions)),
        axes,
        names;
        combine=_valof(spec.combine),
    )
end
