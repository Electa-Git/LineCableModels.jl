function earth!(workspace::LineParametersWorkspace, frequency::Int)
    bindings = workspace.invariants.earth_bindings
    materials = workspace.buffers.earth_materials
    Z, P = workspace.buffers.Zearth, workspace.buffers.Pearth
    for (index, binding) in pairs(bindings.earth_impedance.cases)
        partner = binding.partner == 0 ? nothing :
                  bindings.earth_admittance.cases[binding.partner]
        selection = partner === nothing ? nothing : partner.selection
        earth!(Z, P, binding.selection, selection, binding, partner,
            materials.earth_impedance[index], workspace, frequency)
    end
    for (index, binding) in pairs(bindings.earth_admittance.cases)
        binding.partner == 0 || continue
        earth!(Z, P, nothing, binding.selection, nothing, binding,
            materials.earth_admittance[index], workspace, frequency)
    end
    return workspace
end

function earth!(Z, P, selection::EarthImpedanceFormulation, ::Nothing,
        binding, ::Nothing, materials, workspace, frequency)
    return earth!(Z, selection, binding, materials, workspace.input.jω[frequency],
        workspace, materials.thickness)
end

function earth!(Z, P, ::Nothing, selection::EarthAdmittanceFormulation,
        ::Nothing, binding, materials, workspace, frequency)
    return earth!(P, selection, binding, materials, workspace.input.jω[frequency],
        workspace, materials.thickness)
end

# Ordinary formulas evaluate their bound indexed equations directly.
function earth!(destination, selection, binding, materials, jω, workspace, thickness)
    thickness = media(selection) === Val(:stratified) ? thickness : nothing
    for group in binding.equations
        earth!(destination, selection, group, binding.interactions,
            materials, jω, workspace, thickness)
    end
    return destination
end

function earth!(destination, selection, group::NamedTuple{(:declaration, :indices)},
        interactions, materials, jω, workspace, thickness)
    for index in group.indices
        interaction = interactions[index]
        pair = interaction.pair
        functor = selection(@view(materials.rho[:, index]),
            @view(materials.epsilon[:, index]), @view(materials.mu[:, index]),
            jω, pair, group.declaration; thickness,
            physical_pair = interaction.physical_pair)
        destination[pair.row, pair.column] = functor(workspace)
    end
    return destination
end
