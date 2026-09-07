"""
$(TYPEDEF)

Retain a circular common pipe and the coaxial assemblies in its cavity.
Conductor and assembly indices refer to the enclosing numerical payload.

$(TYPEDFIELDS)
"""
struct PipeAssembly{T <: Real}
    "Innermost conductor row of the enclosing pipe assembly."
    conductor::Int
    "Directly enclosed coaxial assembly indices."
    children::Vector{Int}
    "Physical homogeneous nonmagnetic cavity material."
    material::Material{T}
end
