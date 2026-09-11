
"""
$(TYPEDSIGNATURES)

**Identification.** Static earth material.

**Expression.** Preserve the supplied static soil resistivity, relative permittivity, and
relative permeability at every positive evaluation frequency.

**Reference.** Package policy with no frequency dependence or fitted parameters.
"""
description(::Formula{:default}) = "Static earth material"

function earth_material(
        ::Val{:default}, material::EarthMaterial, frequency::Real,
        values::NamedTuple, options::NamedTuple, workspace
)
    return material
end

computation_options(::FormulaMethod{:default, typeof(earth_material)}) = (;)

:default
