function earth!(workspace::LineParametersWorkspace, frequency::Int,
        calculations::Tuple = workspace.invariants.earth_calculations,
        materials::Tuple = workspace.buffers.earth_materials)
    foreach(calculations, materials) do calculation, material
        earth!(calculation, material, workspace, frequency)
    end
    return workspace
end

# Construct formula state, evaluate indexed coefficients,
# convert the complete matrix when required, then select physical outputs.
function earth!(binding::NamedTuple, materials::NamedTuple, workspace, frequency::Int)
    calculation = binding.selection(materials, binding, workspace, frequency)
    earth!(calculation.coefficients, binding, calculation.state, materials, workspace)
    physical = earth!(binding.selection, calculation, workspace)
    if !isempty(binding.impedance_indices) && physical.impedance !== workspace.buffers.Zearth
        destination = workspace.buffers.Zearth
        for index in binding.impedance_indices
            pair = binding.interactions[index].pair
            destination[pair.row, pair.column] = physical.impedance[pair.row, pair.column]
        end
    end
    if !isempty(binding.potential_indices) && physical.potential !== workspace.buffers.Pearth
        destination = workspace.buffers.Pearth
        for index in binding.potential_indices
            pair = binding.interactions[index].pair
            destination[pair.row, pair.column] = physical.potential[pair.row, pair.column]
        end
    end
    return workspace
end

"""
$(TYPEDSIGNATURES)

Allocate representative indices and diagnostic ranges for `count` earth
interactions. Each traversal clears its warning records and overwrites the
representatives. Numerical contributions stay in the coefficient matrices.
"""
function initialize_buffers(::typeof(earth!), count::Integer)
    return (representatives = zeros(Int, count),
        integral_ranges = Vector{UnitRange{Int}}(undef, count),
        warning_ranges = Vector{UnitRange{Int}}(undef, count), warnings = NamedTuple[])
end

"""
$(TYPEDSIGNATURES)

Evaluate bound indexed earth equations and distribute their scalar coefficients
into aligned matrices. The binding retains `EarthPair` geometry, source-target
layer dispatch, earlier interactions with matching inputs, and selected equation controls.

Interactions share values only when all declared invariant inputs and current
material values agree under `same_physical_state`. Default bindings include
destination indices, preserving equations that use an index numerically.
Every logical integral and warning retains its receiving row and source column.
"""
function earth!(destinations::Tuple{Vararg{AbstractMatrix}}, binding::NamedTuple,
        state::NamedTuple, materials::NamedTuple, workspace)
    work = workspace.buffers.earth_interactions
    length(work.representatives) >= length(binding.interactions) ||
        throw(DimensionMismatch("earth interaction scratch is too small"))
    empty!(work.warnings)
    foreach(binding.equations) do group
        earth!(destinations, binding.selection, group, binding, state, materials, workspace)
    end
    empty!(work.warnings)
    return destinations
end

function earth!(destinations, selection, group, binding, state, materials, workspace)
    work = workspace.buffers.earth_interactions
    observations = workspace.buffers.observations
    for index in group.indices
        interaction = binding.interactions[index]
        pair = interaction.pair
        previous_interaction = binding.previous[index]
        while previous_interaction != 0
            same = let previous_interaction = previous_interaction
                all(
                    values -> same_physical_state(@view(values[:, index]),
                        @view(values[:, previous_interaction])),
                    (materials.rho, materials.epsilon, materials.mu))
            end
            same && break
            previous_interaction = binding.previous[previous_interaction]
        end
        if previous_interaction == 0
            work.representatives[index] = index
            first_integral = observations === nothing ? 1 : length(observations)+1
            first_warning = length(work.warnings)+1
            functor = selection(state, interaction, group.declaration)
            result = functor(workspace)
            values = result isa Number ? (result,) : result
            length(values) == length(destinations) ||
                throw(DimensionMismatch("earth coefficients must match their destinations"))
            foreach(destinations, values) do destination, value
                destination[pair.row, pair.column] = value
            end
            work.integral_ranges[index] = first_integral:(observations === nothing ? 0 :
                                                          length(observations))
            work.warning_ranges[index] = first_warning:length(work.warnings)
        else
            representative = work.representatives[previous_interaction]
            work.representatives[index] = representative
            previous_pair = binding.interactions[representative].pair
            for destination in destinations
                destination[pair.row, pair.column] = destination[previous_pair.row, previous_pair.column]
            end
            context = (receiver = pair.row, source = pair.column)
            if observations !== nothing
                for position in work.integral_ranges[representative]
                    record = observations[position]
                    integral_context = record.context isa NamedTuple ?
                                       merge(record.context, context) : record.context
                    push!(observations, merge(record, (context = integral_context,)))
                end
            end
            for position in work.warning_ranges[representative]
                record = work.warnings[position]
                integral_context = record.context isa NamedTuple ?
                                   merge(record.context, context) : record.context
                record_integral!(nothing, nothing, record.value, record.estimated_error,
                    record.controls, integral_context)
            end
        end
    end
    return destinations
end
