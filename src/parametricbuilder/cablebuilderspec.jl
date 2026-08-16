struct PartBuilder{Role,Part,Mode,C,D,A<:Tuple,M}
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
        material,
    ) where {Role,Part,Mode,A<:Tuple}
        component isa Symbol ||
            throw(ArgumentError("part component must be a Symbol"))
        layers isa Integer && layers > 0 ||
            throw(ArgumentError("part layers must be a positive integer"))
        Mode in (:radius, :solid) && layers != 1 && throw(ArgumentError(
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
            typeof(material),
        }(component, Int(layers), dimension, args, material)
    end
end

_part_role(::PartBuilder{Role}) where {Role} = Role
_part_type(::PartBuilder{Role,Part}) where {Role,Part} = Part
_part_mode(::PartBuilder{Role,Part,Mode}) where {Role,Part,Mode} = Mode
_valof(::Val{Value}) where {Value} = Value

function _simple_part_builder(
    role::Val,
    part::Val,
    mode::Val,
    component,
    layers,
    dimension,
    material,
)
    return PartBuilder(role, part, mode, component, layers, dimension, (), material)
end

function _wire_part_builder(
    role::Val,
    part::Val,
    mode::Val,
    component,
    layers,
    wire_radius,
    num_wires,
    lay_ratio,
    material,
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
        material,
    )
end

function _strip_part_builder(
    role::Val,
    part::Val,
    mode::Val,
    component,
    layers,
    dimension,
    width,
    lay_ratio,
    material,
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
        material,
    )
end

function _part_space(target, axes::Tuple, names::Tuple; combine::Symbol=:product)
    return Gridspace{PartBuilder}(
        target,
        map(_gridspace_axis, axes),
        names;
        combine,
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
    Gridspace, PartBuilder, _part_space, _radial_declaration,
    _simple_part_builder, _wire_part_builder, _strip_part_builder, _valof
import ...DataModel

function Solid(
    component::Symbol;
    radius,
    material,
    combine::Symbol=:product,
)
    return _part_space(
        _simple_part_builder,
        (
            Val(:conductor),
            Val(DataModel.Tubular),
            Val(:solid),
            component,
            1,
            radius,
            material,
        ),
        (:role, :part, :mode, :component, :layers, :radius, :material);
        combine,
    )
end

function Tubular(
    component::Symbol;
    layers=1,
    radius=nothing,
    thickness=nothing,
    material,
    combine::Symbol=:product,
)
    mode, dimension = _radial_declaration(radius, thickness)
    return _part_space(
        _simple_part_builder,
        (
            Val(:conductor),
            Val(DataModel.Tubular),
            mode,
            component,
            layers,
            dimension,
            material,
        ),
        (:role, :part, :mode, :component, :layers, _valof(mode), :material);
        combine,
    )
end

function Wires(
    component::Symbol;
    layers=1,
    wire_radius,
    num_wires,
    lay_ratio=11.0,
    material,
    combine::Symbol=:product,
)
    return _part_space(
        _wire_part_builder,
        (
            Val(:conductor),
            Val(DataModel.CircStrands),
            Val(:wire_radius),
            component,
            layers,
            wire_radius,
            num_wires,
            lay_ratio,
            material,
        ),
        (
            :role,
            :part,
            :mode,
            :component,
            :layers,
            :wire_radius,
            :num_wires,
            :lay_ratio,
            :material,
        );
        combine,
    )
end

function Strip(
    component::Symbol;
    layers=1,
    radius=nothing,
    thickness=nothing,
    width,
    lay_ratio=0.0,
    material,
    combine::Symbol=:product,
)
    mode, dimension = _radial_declaration(radius, thickness)
    return _part_space(
        _strip_part_builder,
        (
            Val(:conductor),
            Val(DataModel.Strip),
            mode,
            component,
            layers,
            dimension,
            width,
            lay_ratio,
            material,
        ),
        (
            :role,
            :part,
            :mode,
            :component,
            :layers,
            _valof(mode),
            :width,
            :lay_ratio,
            :material,
        );
        combine,
    )
end

function Stranded(
    component::Symbol;
    layers::Int,
    wire_radius,
    num_wires::Int=6,
    lay_ratio=11.0,
    material,
    combine::Symbol=:product,
)
    layers >= 1 || throw(ArgumentError("layers must be at least one"))
    central = Wires(
        component;
        layers=1,
        wire_radius,
        num_wires=1,
        lay_ratio=0.0,
        material,
        combine,
    )
    layers == 1 && return (central,)
    rings = Wires(
        component;
        layers=layers - 1,
        wire_radius,
        num_wires,
        lay_ratio,
        material,
        combine,
    )
    return (central, rings)
end

end

module Insulator

using ..ParametricBuilder:
    _part_space, _radial_declaration, _simple_part_builder, _valof
import ...DataModel

function _annular(
    part,
    component::Symbol;
    layers=1,
    radius=nothing,
    thickness=nothing,
    material,
    combine::Symbol=:product,
)
    mode, dimension = _radial_declaration(radius, thickness)
    return _part_space(
        _simple_part_builder,
        (
            Val(:insulator),
            Val(part),
            mode,
            component,
            layers,
            dimension,
            material,
        ),
        (:role, :part, :mode, :component, :layers, _valof(mode), :material);
        combine,
    )
end

Tubular(component::Symbol; kwargs...) =
    _annular(DataModel.Insulator, component; kwargs...)
Semicon(component::Symbol; kwargs...) =
    _annular(DataModel.Semicon, component; kwargs...)

end

function _resolved_outer_radius(builder::PartBuilder, r_in)
    mode = _part_mode(builder)
    mode === :radius && return builder.dimension
    mode === :solid && return builder.dimension
    mode === :thickness && return r_in + builder.dimension
    throw(ArgumentError("part mode :$mode does not define an annular outer radius"))
end

function _materialize_part(builder::PartBuilder, r_in, layer::Int)
    role = _part_role(builder)
    part = _part_type(builder)
    mode = _part_mode(builder)

    if part === DataModel.CircStrands
        num_wires, lay_ratio = builder.args
        layer_count = num_wires == 1 ? 1 : layer * num_wires
        return DataModel.CircStrands(
            r_in,
            builder.dimension,
            layer_count,
            lay_ratio,
            builder.material,
        )
    elseif part === DataModel.Strip
        width, lay_ratio = builder.args
        r_ex = _resolved_outer_radius(builder, r_in)
        return DataModel.Strip(
            r_in,
            r_ex,
            width,
            lay_ratio,
            builder.material,
        )
    elseif part === DataModel.Tubular
        mode === :solid && !iszero(r_in) && throw(ArgumentError(
            "solid conductors must start at radius zero",
        ))
        r_ex = _resolved_outer_radius(builder, r_in)
        return DataModel.Tubular(r_in, r_ex, builder.material)
    elseif part === DataModel.Insulator
        r_ex = _resolved_outer_radius(builder, r_in)
        return DataModel.Insulator(r_in, r_ex, builder.material)
    elseif part === DataModel.Semicon
        r_ex = _resolved_outer_radius(builder, r_in)
        return DataModel.Semicon(r_in, r_ex, builder.material)
    end
    throw(ArgumentError("unsupported materialized part type $part for role :$role"))
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
            part = _materialize_part(builder, current_radius, layer)
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
            part = _materialize_part(builder, current_radius, layer)
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
        insulator_group,
    )
end

struct DesignMaterializer{N}
    nominal::N
end

_nominal_data(::Nothing) = nothing
_nominal_data(nominal::DataModel.NominalData) = nominal
_nominal_data(nominal::NamedTuple) = DataModel.NominalData(; nominal...)
_nominal_data(nominal) = throw(ArgumentError(
    "nominal data must be a NamedTuple, NominalData, or nothing; got $(typeof(nominal))",
))

function (materializer::DesignMaterializer)(identifier, builders...)
    isempty(builders) && throw(ArgumentError("CableBuilder requires at least one part"))
    parts = tuple(builders...)
    names = _component_names(parts)
    isempty(names) && throw(ArgumentError("CableBuilder has no components"))

    first_component = _build_component(first(names), parts, 0.0)
    nominal = _nominal_data(materializer.nominal)
    design = nominal === nothing ?
        DataModel.CableDesign(String(identifier), first_component) :
        DataModel.CableDesign(
            String(identifier),
            first_component;
            nominal_data=nominal,
        )
    for name in Iterators.drop(names, 1)
        base_radius = design.components[end].insulator_group.r_ex
        add!(design, _build_component(name, parts, base_radius))
    end
    return design
end

struct CableDesignSpec{P<:Tuple,N,C} <: AbstractSpec{DataModel.CableDesign}
    identifier::String
    parts::P
    nominal::N
    combine::C
end

_flatten_parts(part::Gridspace{PartBuilder}) = (part,)
function _flatten_parts(parts::Union{Tuple,AbstractVector})
    flattened = ()
    for part in parts
        flattened = (flattened..., _flatten_parts(part)...)
    end
    return flattened
end
_flatten_parts(value) = throw(ArgumentError(
    "CableBuilder expects part Gridspaces; got $(typeof(value))",
))

"""
$(TYPEDSIGNATURES)

Describe a cable design using the radial parts declared by `Conductor` and
`Insulator`. Iterating the returned specification materializes
[`LineCableModels.DataModel.CableDesign`](@ref) objects through the existing
cable-part, group, and component calculations.

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

- A cable-design specification. Ordinary iteration enumerates deterministic
  configurations; `rand` draws a realization from its uncertainty-bearing
  configuration.

# Errors

- Throws `ArgumentError` when no parts are supplied, a value is not a cable
  part declaration, or `combine` is neither `:product` nor `:zip`.
- Materialization reports invalid radial geometry and other physical input
  errors from the materialized cable model.
"""
function CableBuilder(
    identifier::AbstractString,
    parts...;
    nominal=nothing,
    combine::Symbol=:product,
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    flattened = _flatten_parts(parts)
    isempty(flattened) &&
        throw(ArgumentError("CableBuilder requires at least one part"))
    return CableDesignSpec(
        String(identifier),
        flattened,
        nominal,
        Val(combine),
    )
end

function gridspace(spec::CableDesignSpec)
    axes = (_gridspace_axis(spec.identifier), map(_gridspace_axis, spec.parts)...)
    names = (:identifier, ntuple(index -> Symbol(:part_, index), length(spec.parts))...)
    return Gridspace{DataModel.CableDesign}(
        DesignMaterializer(spec.nominal),
        axes,
        names;
        combine=_valof(spec.combine),
    )
end
