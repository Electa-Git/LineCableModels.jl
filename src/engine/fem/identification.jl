"""
Entity identification functions for the FEMTools.jl module.
These functions handle the identification of entities after boolean operations.
"""

"""
$(TYPEDSIGNATURES)

Perform boolean fragmentation on all entities in the model.

# Arguments

- `workspace`: The [`FEMWorkspace`](@ref) containing the entities to fragment.

# Returns

- Nothing. Modifies the Gmsh model in place.

# Examples

```julia
$(FUNCTIONNAME)(workspace)
```

# Notes

This function performs boolean fragmentation on all surfaces and curves in the model.
After fragmentation, the original entities are replaced with new entities that respect
the intersections between them. The original entity tags are no longer valid after
this operation.
"""
function process_fragments(workspace::FEMWorkspace)

    # Get all entities
    surfaces = gmsh.model.get_entities(2)
    curves = gmsh.model.get_entities(1)
    points = gmsh.model.get_entities(0)

    @debug "Initial counts: $(length(surfaces)) surfaces, $(length(curves)) curves, $(length(points)) points"

    # Fragment points onto curves
    if !isempty(curves) && !isempty(points)
        @debug "Fragmenting points onto curves..."
        gmsh.model.occ.fragment(curves, points)
        gmsh.model.occ.synchronize()
    end

    # Get updated entities after first fragmentation
    updated_curves = gmsh.model.get_entities(1)
    updated_points = gmsh.model.get_entities(0)

    @debug "After fragmenting points onto curves: $(length(updated_curves)) curves, $(length(updated_points)) points"

    # Fragment curves onto surfaces
    if !isempty(surfaces) && !isempty(updated_curves)
        @debug "Fragmenting curves onto surfaces..."
        gmsh.model.occ.fragment(surfaces, updated_curves)
        gmsh.model.occ.synchronize()
    end

    # Remove duplicates
    @debug "Removing duplicate entities..."
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # Final counts
    final_surfaces = gmsh.model.get_entities(2)
    final_curves = gmsh.model.get_entities(1)
    final_points = gmsh.model.get_entities(0)

    @info "Boolean fragmentation completed"
    @debug "Before: $(length(surfaces)) surfaces, $(length(curves)) curves, $(length(points)) points"
    @debug "After: $(length(final_surfaces)) surfaces, $(length(final_curves)) curves, $(length(final_points)) points"
    @debug "Unique markers in workspace: $(length(workspace.unassigned_entities)) markers"

end

function identify_by_marker(workspace::FEMWorkspace)

    # Get all surfaces after fragmentation
    all_surfaces = gmsh.model.get_entities(2)

    # Track statistics
    total_entities = length(workspace.unassigned_entities)
    identified_count = 0

    # Copy keys to avoid modifying dict during iteration
    markers = collect(keys(workspace.unassigned_entities))

    # For each marker, find which surface contains it
    for marker in markers
        entity_data = workspace.unassigned_entities[marker]
        physical_group_tag = entity_data.core.physical_group_tag
        elementary_name = entity_data.core.elementary_name

        for (dim, tag) in all_surfaces
            if !(entity_data isa CurveEntity)
                # Check if marker is inside this surface
                if gmsh.model.is_inside(dim, tag, marker) == 1
                    fem_entity = GmshObject(tag, entity_data)

                    # Place in appropriate container
                    if entity_data isa CablePartEntity
                        if entity_data.cable_part isa AbstractConductorPart
                            push!(workspace.conductors, fem_entity)
                        elseif entity_data.cable_part isa AbstractInsulatorPart
                            push!(workspace.insulators, fem_entity)
                        end
                    elseif entity_data isa SurfaceEntity
                        push!(workspace.space_regions, fem_entity)
                    end

                    delete!(workspace.unassigned_entities, marker)
                    identified_count += 1
                    @debug "Marker at $(marker) identified entity $(tag) as $(elementary_name) (tag: $(physical_group_tag))"
                    break
                end
            end
        end
    end

    # Get all remaining curves after fragmentation
    all_curves = gmsh.model.get_entities(1)

    # Update keys to avoid modifying dict during iteration
    markers = collect(keys(workspace.unassigned_entities))

    # For each marker, find which surface contains it
    for marker in markers
        entity_data = workspace.unassigned_entities[marker]
        physical_group_tag = entity_data.core.physical_group_tag
        elementary_name = entity_data.core.elementary_name

        for (dim, tag) in all_curves
            # Check if marker is inside this curve
            if gmsh.model.is_inside(dim, tag, marker) == 1
                # Found match - create GmshObject and add to appropriate container
                fem_entity = GmshObject(tag, entity_data)

                # Place in appropriate container
                if entity_data isa CurveEntity
                    push!(workspace.boundaries, fem_entity)
                end

                delete!(workspace.unassigned_entities, marker)
                identified_count += 1
                @debug "Marker at $(marker) identified entity $(tag) as $(elementary_name) (tag: $(physical_group_tag))"
                break
            end
        end
    end

    # Report identification stats
    @info "Entity identification completed: $(identified_count)/$(total_entities) entities identified"

    if !isempty(workspace.unassigned_entities)
        @warn "$(length(workspace.unassigned_entities))/$(total_entities) markers could not be matched to entities"
    end
end
function assign_physical_groups(workspace::FEMWorkspace)
    # Group entities by physical tag and dimension
    entities_by_physical_group_tag = Dict{Tuple{Int,Int},Vector{Int}}()

    # Process all entity containers
    for container in [workspace.conductors, workspace.insulators, workspace.space_regions, workspace.boundaries]
        for entity in container
            physical_group_tag = entity.data.core.physical_group_tag
            elementary_name = entity.data.core.elementary_name
            dim = entity.data isa CurveEntity ? 1 : 2

            # Key is now a tuple of (physical_group_tag, dimension)
            group_key = (physical_group_tag, dim)

            if !haskey(entities_by_physical_group_tag, group_key)
                entities_by_physical_group_tag[group_key] = Int[]
            end

            # Add this entity to the collection for this physical tag
            current_physical_group = gmsh.model.get_physical_groups_for_entity(dim, entity.tag)
            if !isempty(current_physical_group)
                @debug "Entity $(entity.tag) already has physical group: $(current_physical_group)"
            end
            push!(entities_by_physical_group_tag[group_key], entity.tag)

            if !isempty(elementary_name)
                # Append the complete name to the shape
                current_name = gmsh.model.get_entity_name(dim, entity.tag)
                if !isempty(current_name)
                    @debug "Entity $(entity.tag) already has elementary name: $(current_name)"
                end
                gmsh.model.set_entity_name(dim, entity.tag, elementary_name)
            end
        end
    end

    # Create physical groups for each physical tag
    successful_groups = 0
    failed_groups = 0

    for ((physical_group_tag, dim), entity_tags) in entities_by_physical_group_tag
        try
            physical_group_name = create_physical_group_name(workspace, physical_group_tag)
            @debug "Creating physical group $(physical_group_name) (tag: $(physical_group_tag), dim: $(dim)) with $(length(entity_tags)) entities"

            # Use the correct dimension when creating the physical group
            gmsh.model.add_physical_group(dim, entity_tags, physical_group_tag, physical_group_name)
            successful_groups += 1
        catch e
            @warn "Failed to create physical group tag: $(physical_group_tag), dim: $(dim): $(e)"
            failed_groups += 1
        end
    end

    @info "Physical groups assigned: $(successful_groups) successful, $(failed_groups) failed out of $(length(entities_by_physical_group_tag)) total"
end
