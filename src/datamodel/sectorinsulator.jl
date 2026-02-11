"""
$(TYPEDEF)

Represents an insulating layer surrounding a sector-shaped conductor.

$(TYPEDFIELDS)
"""
struct SectorInsulator{T<:REALSCALAR} <: AbstractInsulatorPart{T}
    "Inner radius (not applicable, defined by inner sector) \\[m\\]."
    radius_in::T
    "Outer radius (equivalent back radius of outer boundary) \\[m\\]."
    radius_ext::T
    "The inner sector conductor that this insulator surrounds."
    inner_sector::Sector{T}
    "The thickness of the insulating layer \\[m\\]."
    thickness::T
    "Material properties of the insulator."
    material_props::Material{T}
    "Operating temperature of the insulator \\[°C\\]."
    temperature::T
    "Cross-sectional area of the insulating layer \\[m²\\]."
    cross_section::T
    "Shunt capacitance (approximated) \\[F/m\\]."
    shunt_capacitance::T
    "Shunt conductance (approximated) \\[S·m\\]."
    shunt_conductance::T
    "Calculated vertices of the outer boundary polygon."
    outer_vertices::Vector{Point{2, T}}
end


function SectorInsulator(
    inner_sector::Sector{T},
    thickness::T,
    material_props::Material{T};
    temperature::T=T₀,
) where {T<:REALSCALAR}

    # 1. Calculate the outer vertices by offsetting the inner sector's geometry
    outer_vertices = _calculate_offset_polygon(inner_sector.vertices, thickness)

    # 2. Calculate areas
    inner_area = inner_sector.cross_section

    # tol = 1e-9
    # if !isempty(outer_vertices)
    #     firstp = outer_vertices[1]
    #     lastp  = outer_vertices[end]
    #     if !(isapprox(firstp[1], lastp[1]; atol=tol, rtol=0.0) &&
    #          isapprox(firstp[2], lastp[2]; atol=tol, rtol=0.0))
    #         push!(outer_vertices, firstp)
    #         @debug "(Sector) outer_vertices not closed — appended first point to close polygon."
    #     end
    # end



    #outer_area = PolygonOps.area(outer_vertices)
    outer_area = _shoelace_area(outer_vertices)
    @debug "SectorInsulator inner area: $(inner_area*1e6) mm²"
    @debug "SectorInsulator outer area: $(outer_area*1e6) mm²"
    cross_section = outer_area - inner_area

    # 3. Approximate capacitance and conductance using equivalent coaxial circles
    # This is more consistent with the package's approach than a parallel plate model.
    r_eq_in = sqrt(inner_area / π)
    r_eq_ext = sqrt(outer_area / π)
    shunt_capacitance = calc_shunt_capacitance(r_eq_in, r_eq_ext, material_props.eps_r)
    shunt_conductance = calc_shunt_conductance(r_eq_in, r_eq_ext, material_props.rho)

    # 4. Determine outer radius from the new params
    outer_r_back = inner_sector.params.r_back + thickness

    return SectorInsulator{T}(
        inner_sector.radius_ext,
        outer_r_back,
        inner_sector,
        thickness,
        material_props,
        temperature,
        cross_section,
        shunt_capacitance,
        shunt_conductance,
        outer_vertices
    )
end

# --- Geometric Helper Functions (internal) --- 
# REVISED: This function now takes vertices and thickness to compute a geometric offset.
function _calculate_offset_polygon(vertices::Vector{Point{2, T}}, thickness::T) where {T<:REALSCALAR}
    num_vertices = length(vertices)
    if num_vertices < 3
        error("Polygon must have at least 3 vertices.")
    end

    new_vertices = similar(vertices)

    for i in 1:num_vertices
        p_prev = vertices[i == 1 ? num_vertices : i - 1]
        p_curr = vertices[i]
        p_next = vertices[i == num_vertices ? 1 : i + 1]

        v1 = p_curr - p_prev
        v2 = p_next - p_curr

        # Normalize the vectors
        v1_norm = v1 / norm(v1)
        v2_norm = v2 / norm(v2)

        # Normal vectors (rotated 90 degrees clockwise for outward direction)
        n1 = Point(v1_norm[2], -v1_norm[1])
        n2 = Point(v2_norm[2], -v2_norm[1])

        # Bisector of the normals
        bisector = (n1 + n2) / norm(n1 + n2)

        # Angle between the two vectors to calculate the correct offset distance
        angle = acos(clamp(v1_norm ⋅ v2_norm, -1.0, 1.0))
        
        # Miter length
        offset_distance = thickness / sin((π - angle) / 2)

        if isinf(offset_distance) # The vectors are parallel
            new_vertices[i] = p_curr + n1 * thickness
        else
            new_vertices[i] = p_curr + bisector * offset_distance
        end
    end

    return new_vertices
end