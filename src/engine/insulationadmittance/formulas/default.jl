
"""
$(TYPEDSIGNATURES)

`:default` is the package routing alias for `:lossless`; it has no separate
constitutive equation.
"""
function description(::Type{<:Formula{:default}}; compact::Bool=false)
    compact ? "Default" : "Default routing alias to lossless cable-insulation admittivity"
end

"""
$(TYPEDSIGNATURES)

Route the default cable-insulation selection to the explicit `:lossless`
constitutive relation.
"""
@inline function insulation_material(
        ::Val{:default}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    insulation_material(
        Val(:lossless), material, frequency, temperature, values, options, workspace
    )
end

computation_options(::FormulaMethod{:default, typeof(insulation_material)}) =
    computation_options(FormulaMethod(Val(:lossless), insulation_material))

:default
