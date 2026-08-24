struct PartBuilder{Role, Part, Mode, C, D, A <: Tuple, M}
    component::C
    layers::Int
    dimension::D
    args::A
    material::M

    function PartBuilder(
            ::Val{Role},
            ::Val{Part},
            ::Val{Mode},
            component,
            layers,
            dimension,
            args::A,
            material
    ) where {Role, Part, Mode, A <: Tuple}
        component isa Symbol ||
            throw(ArgumentError("part component must be a Symbol"))
        layers isa Integer && layers > 0 ||
            throw(ArgumentError("part layers must be a positive integer"))
        Mode in (:radius, :solid, :thickness, :wire_radius) ||
            throw(ArgumentError("unsupported radial mode :$Mode"))
        Mode in (:radius, :solid) && layers != 1 &&
            throw(ArgumentError(
                "absolute radius is valid only for a single layer; use thickness for repeated stacking",
            ))
        dimension isa Real && dimension > zero(dimension) ||
            throw(ArgumentError("part radius/thickness must be a positive real number"))
        material isa Materials.Material ||
            throw(ArgumentError("part material must resolve to Materials.Material"))
        return new{
            Role,
            Part,
            Mode,
            typeof(component),
            typeof(dimension),
            A,
            typeof(material)
        }(component, Int(layers), dimension, args, material)
    end
end

_part_role(::PartBuilder{Role}) where {Role} = Role
_part_type(::PartBuilder{Role, Part}) where {Role, Part} = Part
_part_mode(::PartBuilder{Role, Part, Mode}) where {Role, Part, Mode} = Mode
function PartBuilder(
        role::Val,
        part::Val,
        mode::Val,
        component,
        layers,
        dimension,
        material
)
    return PartBuilder(role, part, mode, component, layers, dimension, (), material)
end

function PartBuilder(
        role::Val,
        part::Val{DataModel.CircStrands},
        mode::Val,
        component,
        layers,
        wire_radius,
        num_wires,
        lay_ratio,
        material
)
    num_wires isa Integer && num_wires > 0 ||
        throw(ArgumentError("num_wires must be a positive integer"))
    lay_ratio isa Real && lay_ratio >= zero(lay_ratio) ||
        throw(ArgumentError("lay_ratio must be a nonnegative real number"))
    return PartBuilder(
        role,
        part,
        mode,
        component,
        layers,
        wire_radius,
        (Int(num_wires), lay_ratio),
        material
    )
end

function PartBuilder(
        role::Val,
        part::Val{DataModel.Strip},
        mode::Val,
        component,
        layers,
        dimension,
        width,
        lay_ratio,
        material
)
    width isa Real && width > zero(width) ||
        throw(ArgumentError("strip width must be a positive real number"))
    lay_ratio isa Real && lay_ratio >= zero(lay_ratio) ||
        throw(ArgumentError("lay_ratio must be a nonnegative real number"))
    return PartBuilder(
        role,
        part,
        mode,
        component,
        layers,
        dimension,
        (width, lay_ratio),
        material
    )
end

function _radial_declaration(radius, thickness)
    xor(radius === nothing, thickness === nothing) || throw(ArgumentError(
        "provide exactly one of radius or thickness",
    ))
    radius !== nothing && return Val(:radius), radius
    return Val(:thickness), thickness
end

module Conductor

using ..ParametricBuilder:
                           AbstractGrid, Grid, Gridspace, PartBuilder,
                           _radial_declaration
import ...DataModel

function Solid(
        component::Symbol;
        radius,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    values = (
        Val(:conductor), Val(DataModel.Tubular), Val(:solid), component, 1,
        radius, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Tubular(
        component::Symbol;
        layers = 1,
        radius = nothing,
        thickness = nothing,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    mode, dimension = _radial_declaration(radius, thickness)
    values = (
        Val(:conductor), Val(DataModel.Tubular), mode, component, layers,
        dimension, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Wires(
        component::Symbol;
        layers = 1,
        wire_radius,
        num_wires,
        lay_ratio = 11.0,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    values = (
        Val(:conductor), Val(DataModel.CircStrands), Val(:wire_radius),
        component, layers, wire_radius, num_wires, lay_ratio, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Strip(
        component::Symbol;
        layers = 1,
        radius = nothing,
        thickness = nothing,
        width,
        lay_ratio = 0.0,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    mode, dimension = _radial_declaration(radius, thickness)
    values = (
        Val(:conductor), Val(DataModel.Strip), mode, component, layers,
        dimension, width, lay_ratio, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

function Stranded(
        component::Symbol;
        layers::Int,
        wire_radius,
        num_wires::Int = 6,
        lay_ratio = 11.0,
        material,
        combine::Symbol = :product
)
    layers >= 1 || throw(ArgumentError("layers must be at least one"))
    central = Wires(
        component;
        layers = 1,
        wire_radius,
        num_wires = 1,
        lay_ratio = 0.0,
        material,
        combine
    )
    layers == 1 && return (central,)
    rings = Wires(
        component;
        layers = layers - 1,
        wire_radius,
        num_wires,
        lay_ratio,
        material,
        combine
    )
    return (central, rings)
end

end

module Insulator

using ..ParametricBuilder:
                           AbstractGrid, Grid, Gridspace, PartBuilder,
                           _radial_declaration
import ...DataModel

function _annular(
        part,
        component::Symbol;
        layers = 1,
        radius = nothing,
        thickness = nothing,
        material,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    mode, dimension = _radial_declaration(radius, thickness)
    values = (
        Val(:insulator), Val(part), mode, component, layers, dimension, material
    )
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{PartBuilder}(PartBuilder, grids; combine)
    end
    return PartBuilder(values...)
end

Tubular(component::Symbol; kwargs...) = _annular(DataModel.Insulator, component; kwargs...)
Semicon(component::Symbol; kwargs...) = _annular(DataModel.Semicon, component; kwargs...)

end

function _resolved_outer_radius(builder::PartBuilder{<:Any, <:Any, :radius}, _)
    builder.dimension
end
_resolved_outer_radius(builder::PartBuilder{<:Any, <:Any, :solid}, _) = builder.dimension
function _resolved_outer_radius(builder::PartBuilder{<:Any, <:Any, :thickness}, r_in)
    r_in + builder.dimension
end

function (builder::PartBuilder{:conductor, DataModel.CircStrands})(
        r_in,
        layer::Int
)
    num_wires, lay_ratio = builder.args
    layer_count = num_wires == 1 ? 1 : layer * num_wires
    r_ex = layer_count == 1 ? builder.dimension : r_in + 2 * builder.dimension
    return DataModel.CircStrands(
        r_in, r_ex, builder.dimension, layer_count, lay_ratio, builder.material
    )
end

function (builder::PartBuilder{:conductor, DataModel.Strip})(
        r_in,
        ::Int
)
    width, lay_ratio = builder.args
    return DataModel.Strip(
        r_in, _resolved_outer_radius(builder, r_in), width, lay_ratio,
        builder.material
    )
end

function (builder::PartBuilder{:conductor, DataModel.Tubular, :solid})(
        r_in,
        ::Int
)
    iszero(r_in) || throw(ArgumentError("solid conductors must start at radius zero"))
    return DataModel.Tubular(r_in, builder.dimension, builder.material)
end

function (builder::PartBuilder{:conductor, DataModel.Tubular})(
        r_in,
        ::Int
)
    return DataModel.Tubular(
        r_in, _resolved_outer_radius(builder, r_in), builder.material
    )
end

function (builder::PartBuilder{:insulator, DataModel.Insulator})(
        r_in,
        ::Int
)
    return DataModel.Insulator(
        r_in, _resolved_outer_radius(builder, r_in), builder.material
    )
end

function (builder::PartBuilder{:insulator, DataModel.Semicon})(
        r_in,
        ::Int
)
    return DataModel.Semicon(
        r_in, _resolved_outer_radius(builder, r_in), builder.material
    )
end

function (builder::PartBuilder)(r_in, layer::Int)
    throw(ArgumentError(
        "unsupported part declaration $(_part_role(builder)) / " *
        "$(_part_type(builder)) / $(_part_mode(builder)) at layer $layer",
    ))
end

function _part_scalar(builder::PartBuilder)
    promote_type(
        typeof(float(builder.dimension)), eltype(builder.material),
        (typeof(float(value))
        for value in builder.args if value isa Real && !(value isa Integer))...
    )
end

function _convert_part(::Type{T}, builder::PartBuilder{
        Role, Part, Mode}) where {
        T <: Real, Role, Part, Mode
}
    args = map(builder.args) do value
        value isa Integer ? value : value isa Real ? convert(T, float(value)) : value
    end
    return PartBuilder(
        Val(Role), Val(Part), Val(Mode), builder.component, builder.layers,
        convert(T, float(builder.dimension)), args,
        convert(Materials.Material{T}, builder.material)
    )
end

function _promote_parts(builders::Tuple)
    T = promote_type((_part_scalar(builder) for builder in builders)...)
    return map(builder -> _convert_part(T, builder), builders), T
end

function _component_names(builders::Tuple)
    names = Symbol[]
    for builder in builders
        builder.component in names || push!(names, builder.component)
    end
    return names
end

function _build_component(component::Symbol, builders::Tuple, base_radius)
    conductor_group = nothing
    current_radius = base_radius
    for builder in builders
        builder.component === component && _part_role(builder) === :conductor ||
            continue
        for layer in 1:builder.layers
            part = builder(current_radius, layer)
            conductor_group = conductor_group === nothing ?
                              DataModel.ConductorGroup(part) :
                              add!(conductor_group, part)
            current_radius = conductor_group.r_ex
        end
    end
    conductor_group === nothing &&
        throw(ArgumentError("component :$component has no conductor layers"))

    insulator_group = nothing
    current_radius = conductor_group.r_ex
    for builder in builders
        builder.component === component && _part_role(builder) === :insulator ||
            continue
        for layer in 1:builder.layers
            part = builder(current_radius, layer)
            insulator_group = insulator_group === nothing ?
                              DataModel.InsulatorGroup(part) :
                              add!(insulator_group, part)
            current_radius = insulator_group.r_ex
        end
    end
    insulator_group === nothing &&
        throw(ArgumentError("component :$component has no insulator layers"))
    return DataModel.CableComponent(
        String(component),
        conductor_group,
        insulator_group
    )
end

struct DesignMaterializer{N}
    nominal::N
end

_nominal_data(::Nothing) = nothing
_nominal_data(nominal::DataModel.NominalData) = nominal
_nominal_data(nominal::NamedTuple) = DataModel.NominalData(; nominal...)
function _nominal_data(nominal)
    throw(ArgumentError(
        "nominal data must be a NamedTuple, NominalData, or nothing; got $(typeof(nominal))",
    ))
end

function (materializer::DesignMaterializer)(identifier, builders...)
    isempty(builders) && throw(ArgumentError("CableBuilder requires at least one part"))
    parts, T = _promote_parts(tuple(builders...))
    names = _component_names(parts)
    isempty(names) && throw(ArgumentError("CableBuilder has no components"))

    first_component = _build_component(first(names), parts, zero(T))
    nominal = _nominal_data(materializer.nominal)
    design = nominal === nothing ?
             DataModel.CableDesign(String(identifier), first_component) :
             DataModel.CableDesign(
        String(identifier),
        first_component;
        nominal_data = nominal
    )
    for name in Iterators.drop(names, 1)
        base_radius = design.components[end].insulator_group.r_ex
        add!(design, _build_component(name, parts, base_radius))
    end
    return design
end

_flatten_parts(part::PartBuilder) = (part,)
_flatten_parts(part::Gridspace{PartBuilder}) = (part,)
function _flatten_parts(parts::Union{Tuple, AbstractVector})
    tuple(Iterators.flatten(map(_flatten_parts, parts))...)
end
function _flatten_parts(value)
    throw(ArgumentError(
        "CableBuilder expects conductor or insulator part builders; got $(typeof(value))",
    ))
end

"""
$(TYPEDSIGNATURES)

Construct a cable design from radial parts declared by `Conductor` and
`Insulator`.

# Arguments

- `identifier`: Name assigned to each materialized cable design.
- `parts`: One or more part declarations, supplied in radial order. Nested
  tuples and vectors of declarations are flattened without changing that
  order.

# Keywords

- `nominal=nothing`: Optional reference data. A named tuple may contain the
  fields `designation_code`, `U0`, `U`, `conductor_cross_section`,
  `screen_cross_section`, `armor_cross_section`, `resistance`, `capacitance`,
  and `inductance`. Voltages are expressed in \\[kV\\], cross-sectional areas in
  \\[mm²\\], resistance in \\[Ω/km\\], capacitance in \\[μF/km\\], and inductance in
  \\[mH/km\\].
- `combine=:product`: Local composition rule for direct varying inputs. Use
  `:zip` to pair compatible axes.

# Returns

- A [`DataModel.CableDesign`](@ref) when every part is scalar, or a
  `Gridspace{CableDesign}` when a part varies explicitly.

# Errors

- Throws `ArgumentError` when no parts are supplied, a value is not a cable
  part declaration, or `combine` is neither `:product` nor `:zip`.
- Construction reports invalid radial geometry and other physical input errors
  from the cable model.
"""
function CableBuilder(
        identifier::AbstractString,
        parts...;
        nominal = nothing,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    flattened = _flatten_parts(parts)
    isempty(flattened) &&
        throw(ArgumentError("CableBuilder requires at least one part"))
    build = DesignMaterializer(nominal)
    values = (String(identifier), flattened...)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.CableDesign}(build, grids; combine)
    end
    return build(values...)
end
