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
