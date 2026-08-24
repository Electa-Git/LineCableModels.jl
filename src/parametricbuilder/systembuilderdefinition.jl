struct EarthDefinition{R, E, M, H, C} <: _AbstractDefinition{EarthProps.EarthModel}
    rho::R
    eps_r::E
    mu_r::M
    thickness::H
    combine::C
end

function _earth_model(rho, eps_r, mu_r, thickness)
    EarthProps.EarthModel(rho, eps_r, mu_r; thickness)
end
_earth_model(rho, eps_r, mu_r, ::Nothing) = EarthProps.EarthModel(rho, eps_r, mu_r)

function Gridspace(spec::EarthDefinition)
    return Gridspace{EarthProps.EarthModel}(
        _earth_model,
        map(_gridspace_axis, (spec.rho, spec.eps_r, spec.mu_r, spec.thickness));
        combine = _valof(spec.combine)
    )
end

"""
$(TYPEDSIGNATURES)

Describe static earth properties independently of analysis frequencies. The
materialized `EarthModel` is constructed when the parent problem is materialized;
frequency-dependent evaluation remains an Engine formulation choice.
"""
function Earth(;
        rho,
        eps_r = 1.0,
        mu_r = 1.0,
        thickness = nothing,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    return Gridspace(EarthDefinition(rho, eps_r, mu_r, thickness, Val(combine)))
end

function _position_coordinates(
        position::PositionDefinition,
        design::DataModel.CableDesign
)
    kind = _position_kind(position)
    if kind === :point
        x, y = position.parameters
        return ((x, y, only(position.connections)),)
    end

    x, y, spacing = position.parameters
    spacing > zero(spacing) ||
        throw(ArgumentError("formation spacing must be positive"))
    minimum_spacing = 2 * DataModel.outer_radius(design)
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
            length(position.connections)
        )
    elseif kind === :vflat
        ntuple(
            index -> (x, y - (index - 1) * spacing),
            length(position.connections)
        )
    else
        throw(ArgumentError("unsupported position kind :$kind"))
    end
    return ntuple(
        index -> (
            coordinates[index][1],
            coordinates[index][2],
            position.connections[index]
        ),
        length(position.connections)
    )
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

    earth isa EarthProps.EarthModel || throw(ArgumentError(
        "earth must resolve to a static EarthModel",
    ))
    scalar_types = Any[
        eltype(design), eltype(earth), typeof(float(line_length)),
        typeof(float(temperature)), eltype(frequencies)
    ]
    for position in positions, value in position.parameters

        value isa Real && push!(scalar_types, typeof(float(value)))
    end
    T = promote_type(scalar_types...)
    design = convert(DataModel.CableDesign{T}, design)
    earth = convert(EarthProps.EarthModel{T}, earth)
    line_length = convert(T, float(line_length))
    temperature = convert(T, float(temperature))
    frequencies = T[convert(T, float(value)) for value in frequencies]

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
        DataModel.CablePosition(design, x, y, connections)
    )
    for (next_x, next_y, next_connections) in Iterators.drop(coordinates, 1)
        system = add!(system, design, next_x, next_y, next_connections)
    end

    return Engine.LineParametersProblem(
        system;
        temperature,
        earth_props = earth,
        frequencies
    )
end

struct SystemDefinition{D, P <: Tuple, L, T, E, F, C} <:
       _AbstractDefinition{Engine.LineParametersProblem}
    identifier::String
    design::D
    positions::P
    line_length::L
    temperature::T
    earth::E
    frequencies::F
    combine::C
end

_flatten_positions(position::Gridspace{PositionDefinition}) = (position,)
function _flatten_positions(positions::Union{Tuple, AbstractVector})
    tuple(Iterators.flatten(map(_flatten_positions, positions))...)
end
function _flatten_positions(value)
    throw(ArgumentError(
        "SystemBuilder positions must be created by at/trifoil/hflat/vflat; got $(typeof(value))",
    ))
end

"""
$(TYPEDSIGNATURES)

Describe a line-cable analysis problem from a cable design, its spatial
arrangement, earth properties, and analysis frequencies. Passing the returned
definition to [`Gridspace`](@ref) creates a lazy space of
[`LineCableModels.Engine.LineParametersProblem`](@ref) objects without running
the numerical analysis.

# Arguments

- `identifier`: Name assigned to the materialized cable system.
- `design`: A materialized cable design or a [`Gridspace`](@ref) returned by
  [`CableBuilder`](@ref).
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

- A line-parameter problem [`Gridspace`](@ref). Iterate it for primitive
  problems, or pass it to [`ParametricProblem`](@ref) with an explicit
  higher-order formulation.

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
        length = 1000.0,
        temperature = 20.0,
        earth,
        frequencies,
        combine::Symbol = :product
)
    design isa Union{DataModel.CableDesign, Gridspace{DataModel.CableDesign}} ||
        throw(ArgumentError("design must be a materialized CableDesign or CableBuilder result"))
    earth isa Union{EarthProps.EarthModel, Gridspace{EarthProps.EarthModel}} ||
        throw(ArgumentError("earth must be a materialized EarthModel or Earth result"))
    frequencies isa Union{AbstractVector, AbstractGrid} || throw(ArgumentError(
        "frequencies must be an ordinary vector or a Grid of complete vectors",
    ))
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    position_tuple = _flatten_positions(positions)
    isempty(position_tuple) &&
        throw(ArgumentError("SystemBuilder requires at least one position"))
    return Gridspace(SystemDefinition(
        String(identifier),
        design,
        position_tuple,
        length,
        temperature,
        earth,
        frequencies,
        Val(combine)
    ))
end

function Gridspace(spec::SystemDefinition)
    axes = (
        _gridspace_axis(spec.identifier),
        _gridspace_axis(spec.design),
        map(_gridspace_axis, spec.positions)...,
        _gridspace_axis(spec.line_length),
        _gridspace_axis(spec.temperature),
        _gridspace_axis(spec.earth),
        _gridspace_axis(
            spec.frequencies isa AbstractGrid ?
            spec.frequencies : Grid((spec.frequencies,))
        )
    )
    return Gridspace{Engine.LineParametersProblem}(
        SystemMaterializer(length(spec.positions)),
        axes;
        combine = _valof(spec.combine)
    )
end
