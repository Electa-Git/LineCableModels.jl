
"""
$(TYPEDSIGNATURES)

`:default` is the package routing alias for `:lossless`; it has no separate
constitutive equation.
"""
function description(::Type{<:Formula{:default}}; compact::Bool=false)
    compact ? "Default" : "Default routing alias to lossless semiconducting-screen admittivity"
end

"""
$(TYPEDSIGNATURES)

Route the default semiconducting-screen selection to the explicit `:lossless`
constitutive relation.
"""
@inline function semicon_material(
        ::Val{:default}, material::Material{T}, frequency::T,
        temperature::T, values::NamedTuple, options::NamedTuple, workspace
) where {T <: Real}
    semicon_material(
        Val(:lossless), material, frequency, temperature, values, options, workspace
    )
end

computation_options(::FormulaMethod{:default, typeof(semicon_material)}) =
    computation_options(FormulaMethod(Val(:lossless), semicon_material))

:default
