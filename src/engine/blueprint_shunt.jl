"""
$(TYPEDSIGNATURES)

Complete the internal potential coefficients of each formulation for the flattened designs
of one system: the capacitance and potential-coefficient blocks of the cable interiors, which
pair with the exterior potential coefficients of the earth. A dielectric group of the
`:boundary` model reuses an earlier block of the same call, also of another formulation,
when its physical state and numerical controls agree.

# Arguments

- `designs`: completed physical cable designs.
- `blueprints`: the blueprints of `designs`, in the same order.
- `formulations`: the selected formulations in a vector or a tuple. A tuple gives a tuple
  of results.

# Returns

- One result per formulation, in formulation order. `blocks` holds the completed
  internal capacitances in system conductor indices, `covered` marks the radial intervals
  that they replace, and `details` holds the local-model outcomes in system terminal
  order.
"""
function internal_potential_coefficients(designs::AbstractVector,
        blueprints::AbstractVector{<:CableBlueprint{T}},
        formulations::AbstractVector{<:AbstractFormulation}) where {T <: Real}
    Response = NamedTuple{(:blocks, :covered, :details),
        Tuple{Vector{InternalCapacitance{T}}, BitVector, NamedTuple}}
    responses = Vector{Response}(undef, length(formulations))
    # A block holds system indices and matrices, not the group and controls that reuse
    # compares. Each formulation with groups keeps them here, beside its blocks.
    groups = Any[]
    for (index, formulation) in pairs(formulations)
        responses[index] = internal_potential_coefficients(
            formulation.formulas.shunt_model, designs, blueprints, formulation.formulas,
            formulations, responses, groups, index)
    end
    return responses
end

# A tuple of formulations gives a tuple of results.
function internal_potential_coefficients(designs::AbstractVector,
        blueprints::AbstractVector{<:CableBlueprint},
        formulations::Tuple{Vararg{AbstractFormulation, N}}) where {N}
    responses = internal_potential_coefficients(designs, blueprints, collect(formulations))
    return ntuple(index -> responses[index], Val(N))
end

# Complete the internal capacitances of the formulation at `index` for every design. Each
# group reuses a block of the earlier `responses`, whose groups `groups` keeps, or of this
# formulation, when the physical state and controls agree. Any other group is solved once.
# The formulation's groups then join `groups`. The blocks move from cable-local to system
# indices.
function internal_potential_coefficients(selected::ShuntModelFormulation,
        designs::AbstractVector, blueprints::AbstractVector{<:CableBlueprint{T}},
        formulas::NamedTuple, formulations::AbstractVector, responses::AbstractVector,
        groups::AbstractVector, index::Int) where {T <: Real}
    parts = [Expression(selected, ShuntModel.dielectric_groups)(
                 design, blueprint, formulas, design_index)
             for (design_index, (design, blueprint)) in enumerate(zip(designs, blueprints))]
    # A solve reads the group's descriptors, its layer count and the terminal and layer of
    # each open conductor.
    same(a, b) = length(a.layers) === length(b.layers) &&
                 length(a.wires) === length(b.wires) && length(a.tapes) === length(b.tapes) &&
                 all(i -> a.wires[i].terminal === b.wires[i].terminal &&
                     a.wires[i].layer === b.wires[i].layer, eachindex(a.wires)) &&
                 all(i -> a.tapes[i].terminal === b.tapes[i].terminal &&
                     a.tapes[i].layer === b.tapes[i].layer, eachindex(a.tapes)) &&
                 all(pair -> first(pair) === last(pair), zip(Vector(a), Vector(b)))
    # The block that the group at `position` reuses comes from an earlier formulation with the
    # same controls or from an earlier group of this formulation. Without one, it is solved.
    function reused(group, position)
        for (formulation, solved) in groups,
            (block, other) in enumerate(Iterators.flatten(part.groups for part in solved))

            formulations[formulation].formulas.shunt_model === selected &&
                same(other, group) && return (formulation, block)
        end
        for (block, other) in enumerate(Iterators.flatten(part.groups for part in parts))
            block < position || break
            same(other, group) && return (index, block)
        end
        return nothing
    end
    pending = similar(first(parts).groups, 0)
    position = 0
    for part in parts, group in part.groups
        position += 1
        reused(group, position) === nothing && push!(pending, group)
    end
    solved = isempty(pending) ? nothing : _internal_capacitances(selected, pending)
    blocks = InternalCapacitance{T}[]
    covered = falses(sum(length, blueprints))
    reports = similar(first(parts).reports, 0)
    diagnostics = ShuntModel.InternalShuntDiagnostic[]
    position = 0
    next = 0
    offset = 0
    for (design_index, (blueprint, part)) in enumerate(zip(blueprints, parts))
        for group in part.groups
            position += 1
            source = reused(group, position)
            C, P, diagnostic = if source === nothing
                next += 1
                solved[next]
            elseif first(source) < index
                response = responses[first(source)]
                block = response.blocks[last(source)]
                (block.C, block.P, response.details.diagnostics[last(source)])
            else
                (blocks[last(source)].C, blocks[last(source)].P, diagnostics[last(source)])
            end
            block = validate(InternalCapacitance(group.assembly, group.terminals,
                group.reference, C, P), blueprint)
            terminals = block.terminals .+ offset
            push!(blocks, InternalCapacitance(block.assembly .+ offset, terminals,
                block.reference, block.C, block.P))
            covered[first(terminals):(first(terminals)+size(block.C, 1)-1)] .= true
            push!(diagnostics, diagnostic)
        end
        for report in part.reports
            push!(reports, merge(report,
                (; design = design_index, terminals = report.terminals .+ offset)))
        end
        offset += length(blueprint)
    end
    solves = count(i -> findfirst(other -> other.C === blocks[i].C, blocks) == i,
        eachindex(blocks))
    details = (requested = formula_id(selected), solves, domains = reports, diagnostics)
    any(part -> !isempty(part.groups), parts) && push!(groups, (index, parts))
    return (; blocks, covered, details)
end

# Solve each new dielectric group once for its internal capacitance and potential
# coefficients. The element type of the results is the type of the first result.
function _internal_capacitances(selected::ShuntModelFormulation, pending::AbstractVector)
    first_solved = internal_potential_coefficients(selected, first(pending))
    solved = Vector{typeof(first_solved)}(undef, length(pending))
    solved[1] = first_solved
    return _internal_capacitances(solved, selected, pending)
end

# Solve the groups after the first into `solved`.
function _internal_capacitances(solved::AbstractVector, selected::ShuntModelFormulation,
        pending::AbstractVector)
    for index in 2:length(pending)
        solved[index] = internal_potential_coefficients(selected, pending[index])
    end
    return solved
end

# Solve one dielectric group for its internal capacitance, its potential coefficients and its
# numerical record.
function internal_potential_coefficients(selected::ShuntModelFormulation,
        group::ShuntModel.DielectricGroup{T}) where {T}
    result = Expression(selected, ShuntModel.internal_capacitance)(group)
    C = Matrix{T}(result.C)
    P = lu(C) \ Matrix{T}(I, size(C, 1), size(C, 1))
    return (C, P, result.diagnostic)
end

function validate(block::InternalCapacitance, blueprint::CableBlueprint)
    block.assembly in blueprint.assembly_ranges || throw(ArgumentError(
        "a shunt block must belong to a conductor assembly"))
    first(block.terminals) in block.assembly &&
    last(block.terminals) in block.assembly ||
        throw(ArgumentError("shunt block terminals must lie in their assembly"))
    validate(block.reference, blueprint, block)
    return block
end

function validate(reference::ShuntShield, blueprint::CableBlueprint, block::InternalCapacitance)
    n = length(block.terminals) - 1
    size(block.C) == size(block.P) == (n, n) || throw(DimensionMismatch(
        "shield-relative coefficients must match the nonreference terminals"))
    reference.radius == blueprint.conductors[last(block.terminals)].r_in ||
        throw(ArgumentError("the shunt reference must match the shield's inner radius"))
    return reference
end

function validate(reference::ShuntJacket, blueprint::CableBlueprint, block::InternalCapacitance)
    n = length(block.terminals)
    size(block.C) == size(block.P) == (n, n) || throw(DimensionMismatch(
        "jacket-relative coefficients must match every real terminal"))
    last(block.terminals) == last(block.assembly) || throw(ArgumentError(
        "a jacket reference must enclose the assembly's outermost conductor"))
    layers = blueprint[last(block.assembly)].layers
    !isempty(layers) && reference.radius == last(layers).r_ex || throw(ArgumentError(
        "the shunt reference must match the jacket's outer radius"))
    return reference
end

# Install completed terminal coefficients. No boundary equation is evaluated.
function _shunt_potential!(destination, blocks::AbstractVector{<:InternalCapacitance})
    for block in blocks
        inner = first(block.terminals)
        outer = inner+size(block.P, 1)-1
        # Charge on the inner anchor includes all conductors shielded inside it.
        @inbounds for j in first(block.assembly):outer,
            i in first(block.assembly):outer

            destination[i, j] += block.P[max(1, i-inner+1), max(1, j-inner+1)]
        end
    end
    return destination
end

function _shunt_admittance!(destination, blocks::AbstractVector{<:InternalCapacitance}, s)
    for block in blocks
        _shunt_admittance!(block.reference, destination, block, s)
    end
    return destination
end

function _shunt_admittance!(::ShuntShield, destination, block, s)
    inner, reference = first(block.terminals), last(block.terminals)
    @inbounds for j in axes(block.C, 2), i in axes(block.C, 1)
        value = s*block.C[i, j]
        row, column = inner+i-1, inner+j-1
        destination[row, column] += value
        destination[row, reference] -= value
        destination[reference, column] -= value
        destination[reference, reference] += value
    end
    return destination
end

function _shunt_admittance!(::ShuntJacket, destination, block, s)
    inner = first(block.terminals)
    @inbounds for j in axes(block.C, 2), i in axes(block.C, 1)
        destination[inner+i-1, inner+j-1] += s*block.C[i, j]
    end
    return destination
end
