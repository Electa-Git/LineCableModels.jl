"""
$(TYPEDSIGNATURES)

Construct cable blueprints for each selected formulation. Identical local
selections share the completed blueprints; equivalent lossless domains share
their coefficient matrices. Earth-return choices do not enter this calculation.
The returned outer vector follows formulation order, and each inner vector
follows design order. Sharing is confined to this construction call.
"""
function flatten(engine::LineCableModelsCoaxial, designs::AbstractVector,
        ::Type{T}, formulations::AbstractVector{<:AbstractFormulation}) where {T <: Real}
    solutions = NamedTuple[]
    selections = NamedTuple[]
    blueprints = Vector{CableBlueprint{T}}[]
    for formulation in formulations
        methods = formulation.methods
        selected = blueprint_dependencies(methods.shunt_model, methods)
        previous = findfirst(value -> isequal(value, selected), selections)
        current = previous === nothing ?
                  CableBlueprint{T}[flatten(engine, design, T, selected, solutions, index)
                                    for (index, design) in pairs(designs)] :
                  blueprints[previous]
        push!(selections, selected)
        push!(blueprints, current)
    end
    return blueprints
end

# Install completed terminal coefficients; no boundary equation is evaluated.
function _shunt_potential!(destination, blocks::AbstractVector{<:InternalShuntBlock})
    for block in blocks
        inner, reference = first(block.terminals), last(block.terminals)
        # Charge on the inner anchor includes all conductors shielded inside it.
        @inbounds for j in first(block.assembly):(reference - 1),
            i in first(block.assembly):(reference - 1)

            destination[i, j] += block.P[max(1, i-inner+1), max(1, j-inner+1)]
        end
    end
    return destination
end

function _shunt_admittance!(destination, blocks::AbstractVector{<:InternalShuntBlock}, s)
    for block in blocks
        first_index, reference = first(block.terminals), last(block.terminals)
        @inbounds for j in axes(block.C, 2), i in axes(block.C, 1)

            value = s*block.C[i, j]
            row, column = first_index+i-1, first_index+j-1
            destination[row, column] += value
            destination[row, reference] -= value
            destination[reference, column] -= value
            destination[reference, reference] += value
        end
    end
    return destination
end
