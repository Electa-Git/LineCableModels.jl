"""
$(TYPEDSIGNATURES)

Complete the selected shunt model of each formulation for the flattened designs of one
system. Solves of the `:boundary` model share their coefficient matrices when their
physical inputs and numerical controls agree, within this call only.

# Arguments

- `designs`: completed physical cable designs.
- `blueprints`: the blueprints of `designs`, in the same order.
- `formulations`: selected formulations.

# Returns

- One result per formulation, in formulation order. `blocks` holds the completed
  coefficient blocks in system conductor indices, `covered` marks the radial intervals
  that they replace, and `details` holds the local-model outcomes in system terminal
  order.
"""
function internal_shunt_response(designs::AbstractVector,
        blueprints::AbstractVector{<:CableBlueprint{T}},
        formulations::AbstractVector{<:AbstractFormulation}) where {T <: Real}
    Response = NamedTuple{(:blocks, :covered, :details),
        Tuple{Vector{InternalShuntBlock{T}}, BitVector, NamedTuple}}
    solutions = NamedTuple[]
    return Response[internal_shunt_response(formulation.methods.shunt_model, designs,
                        blueprints, formulation.methods, solutions)
                    for formulation in formulations]
end

# Complete one formulation's shunt model for every design, and move each design's
# response from cable-local to system indices.
function internal_shunt_response(selected::ShuntModelFormulation, designs::AbstractVector,
        blueprints::AbstractVector{<:CableBlueprint{T}}, methods::NamedTuple,
        solutions::Vector) where {T <: Real}
    parts = [internal_shunt_response(selected, design, blueprint, T, methods, solutions,
                 design_index)
             for (design_index, (design, blueprint)) in enumerate(zip(designs, blueprints))]
    blocks = InternalShuntBlock{T}[]
    covered = falses(sum(length, blueprints))
    reports = similar(first(parts).details.domains, 0)
    diagnostics = similar(first(parts).details.diagnostics, 0)
    solved = Base.IdSet{Matrix{T}}()
    offset = 0
    for (design_index, (blueprint, part)) in enumerate(zip(blueprints, parts))
        for block in part.blocks
            block.assembly in blueprint.assembly_ranges || throw(ArgumentError(
                "a shunt block must belong to a conductor assembly"))
            first(block.terminals) in block.assembly &&
            last(block.terminals) in block.assembly ||
                throw(ArgumentError("shunt block terminals must lie in their assembly"))
            n = length(block.terminals) - 1
            size(block.C) == size(block.P) == (n, n) || throw(DimensionMismatch(
                "shunt block coefficients must match the nonreference terminals"))
            terminals = block.terminals .+ offset
            push!(blocks, InternalShuntBlock(block.assembly .+ offset, terminals,
                block.C, block.P))
            covered[first(terminals):(last(terminals) - 1)] .= true
            push!(solved, block.C)
        end
        append!(diagnostics, part.details.diagnostics)
        for report in part.details.domains
            push!(reports, merge(report,
                (; design = design_index, terminals = report.terminals .+ offset)))
        end
        offset += length(blueprint)
    end
    details = (requested = first(parts).details.requested,
        effective = isempty(reports) ? :equivalent :
                    all(r -> r.effective === :boundary, reports) ? :boundary :
                    all(r -> r.effective === :equivalent, reports) ? :equivalent : :mixed,
        solves = length(solved), domains = reports, diagnostics)
    return (; blocks, covered, details)
end

# Install completed terminal coefficients. No boundary equation is evaluated.
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
