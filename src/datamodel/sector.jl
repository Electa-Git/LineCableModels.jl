"""
$(TYPEDEF)

Holds the geometric parameters that define the shape of a sector conductor.

$(TYPEDFIELDS)
"""
struct SectorParams{T<:REALSCALAR}
    "Number of sectors in the full cable (e.g., 3 or 4)."
    n_sectors::Int
    "Back radius of the sector (outermost curve) \\[m\\]."
    r_back::T
    "Depth of the sector from the back to the base \\[m\\]."
    d_sector::T
    "Corner radius for rounding sharp edges \\[m\\]."
    r_corner::T
    "Angular width of the conductor's flat base/sides [degrees]."
    theta_cond_deg::T
    "Insulation thickness \\[m\\]."
    d_insulation::T  # needed to correct for the offset
end

"""
$(TYPEDEF)

Represents a single sector-shaped conductor with defined geometric and material properties.

$(TYPEDFIELDS)
"""
struct Sector{T<:REALSCALAR} <: AbstractConductorPart{T}
    "Inner radius (not applicable, typically 0 for the central point) \\[m\\]."
    radius_in::T
    "Outer radius (equivalent back radius) \\[m\\]."
    radius_ext::T
    "Geometric parameters defining the sector's shape."
    params::SectorParams{T}
    "Rotation angle of this specific sector around the cable's center [degrees]."
    rotation_angle_deg::T
    "Material properties of the conductor."
    material_props::Material{T}
    "Operating temperature of the conductor \\[°C\\]."
    temperature::T
    "Cross-sectional area of the sector \\[m²\\]."
    cross_section::T
    "Electrical resistance of the sector \\[Ω/m\\]."
    resistance::T
    "Geometric mean radius (GMR) of the sector (approximated) \\[m\\]."
    gmr::T
    "Calculated vertices defining the polygon shape."
    vertices::Vector{Point{2,T}}
    "Geometric centroid of the sector shape."
    centroid::Point{2,T}
end


function Sector(
    params::SectorParams{T},
    rotation_angle_deg::T,
    material_props::Material{T};
    temperature::T=T₀,
) where {T<:REALSCALAR}

    # 1. Calculate the geometry and vertices for a base (unrotated) sector
    base_vertices = _calculate_sector_polygon_points(params)

    # 2. Rotate the vertices to the specified angle
    rotation_angle_rad = deg2rad(rotation_angle_deg)
    rotated_vertices = [_rotate_point(p, rotation_angle_rad) for p in base_vertices]

    # Ensure the polygon is closed: first point == last point (within tolerance)
    # tol = 1e-9
    # if !isempty(rotated_vertices)
    #     firstp = rotated_vertices[1]
    #     lastp  = rotated_vertices[end]
    #     if !(isapprox(firstp[1], lastp[1]; atol=tol, rtol=0.0) &&
    #          isapprox(firstp[2], lastp[2]; atol=tol, rtol=0.0))
    #         push!(rotated_vertices, firstp)
    #         @debug "(Sector) rotated_vertices not closed — appended first point to close polygon."
    #     end
    # end
    # 3. Calculate cross-sectional area using the Shoelace formula
    #cross_section = PolygonOps.area(rotated_vertices)
    cross_section = _shoelace_area(rotated_vertices)
    #@debug "Sector cross-sectional area: $(cross_section*1e6) mm²"
    @debug "Sector cross-sectional area (Shoelace): $(cross_section*1e6) mm²"
    
    # Calculate centroid
    centroid = _calculate_polygon_centroid(rotated_vertices)
    @debug "Sector numerically calculated centroid point is: $(centroid)"
    # 4. Calculate DC resistance
    rho_eff = calc_temperature_correction(material_props.alpha, temperature, material_props.T0) * material_props.rho
    resistance = rho_eff / cross_section

    # 5. Approximate GMR based on a circle of equivalent area
    r_equiv = sqrt(cross_section / π)
    gmr = r_equiv * exp(-0.25 * material_props.mu_r) # GMR of an equivalent solid round conductor

    return Sector{T}(
        zero(T),
        params.r_back,
        params,
        rotation_angle_deg,
        material_props,
        temperature,
        cross_section,
        resistance,
        gmr,
        rotated_vertices,
        centroid,
    )
end



# --- Geometric Helper Functions (internal) ---

function _calculate_sector_geometry(p::SectorParams)
    phi_deg = 90.0 - p.theta_cond_deg / 2.0
    phi_rad = deg2rad(phi_deg) # ϕ
    @debug "(Sector) phi_deg: $phi_deg"
    @debug "(Sector) phi_rad: $phi_rad rad"

    if abs(cos(phi_rad)) < 1e-9
        error("theta_cond_deg is too close to 180, leading to division by zero. Check parameters.")
    end

    d_base_corner = p.r_corner * (1.0 / cos(phi_rad) - 1.0) # D_B 
    @debug "(Sector) d_base_corner (D_B) : $d_base_corner m"
    d_offset = p.r_back - p.d_sector - d_base_corner  # D_O
    @debug "(Sector) d_offset (D_O) : $d_offset m"
    d_insulation_offset = (p.d_insulation / (cos((pi - ((2 * pi) / p.n_sectors)) / 2.0))) - d_offset # D_I
    @debug "(Sector) d_insulation_offset (D_I) : $d_insulation_offset m"

    # trick test (works though different that Urquhart's)
    #d_offset = (p.d_insulation / (cos((pi - ((2 * pi) / p.n_sectors)) / 2.0))) # D_I
    #d_insulation_offset = 0.0 

    x_base_corner = p.r_corner * sin(phi_rad) # X_A = -X_F
    y_base_corner = x_base_corner * tan(phi_rad) + d_offset # Y_A = Y_F
    node_A = Point2f(x_base_corner, y_base_corner + d_insulation_offset)

    k = p.r_corner / cos(phi_rad) + d_offset
    qa = 1.0 + tan(phi_rad)^2
    qb = 2.0 * k * tan(phi_rad)
    qc = k^2 - (p.r_back - p.r_corner)^2

    discriminant = qb^2 - 4.0 * qa * qc
    if discriminant < 0
        error("Cannot calculate side corner center: negative discriminant ($discriminant). Check parameters.")
    end

    x_side_center = (-qb + sqrt(discriminant)) / (2.0 * qa) # X_N
    @debug "(Sector) x_side_center (X_N): $(x_side_center*1e3) mm"
    @debug "(Sector) r_corner: $(p.r_corner*1e3) mm"
    y_side_center = x_side_center * tan(phi_rad) + k # Y_N

    # Sector width (for checks):
    # w = 2 * X_N + 2 * r_corner
    w_sector = 2.0 * x_side_center + 2.0 * p.r_corner
    @debug "(Sector) sector width (computed): $(w_sector) m ($(w_sector*1e3) mm)"

    x_side_lower = x_side_center + p.r_corner * sin(phi_rad)  # X_B
    y_side_lower = y_side_center - p.r_corner * cos(phi_rad)  # Y_B
    node_B = Point2f(x_side_lower, y_side_lower + d_insulation_offset)

    dist_origin_side_center = sqrt(x_side_center^2 + y_side_center^2)
    if dist_origin_side_center < 1e-9
        error("Side corner center is at the origin, cannot determine upper point direction.")
    end
    x_side_upper = x_side_center * p.r_back / dist_origin_side_center # X_C
    y_side_upper = y_side_center * p.r_back / dist_origin_side_center # Y_C
    node_C = Point2f(x_side_upper, y_side_upper + d_insulation_offset)

    nodes = (
        A=node_A,
        B=node_B,
        C=node_C,
        D=Point2f(-node_C[1], node_C[2]),
        E=Point2f(-node_B[1], node_B[2]),
        F=Point2f(-node_A[1], node_A[2])
    )
    @debug "(Sector) Nodes: $nodes"
    centers = (
        Back=Point2f(0, 0+ d_insulation_offset),
        Base=Point2f(0, d_offset + p.r_corner / cos(phi_rad) + d_insulation_offset),
        RightSide=Point2f(x_side_center, y_side_center+ d_insulation_offset),
        LeftSide=Point2f(-x_side_center, y_side_center+ d_insulation_offset)
    )
    @debug "(Sector) Centers: $centers"
    return (Nodes=nodes, Centers=centers, Params=p)
end

function _generate_arc_points(center, radius, start_angle, end_angle, num_points)
    while end_angle < start_angle
        @debug "(Sector) end_angle < start_angle: $(rad2deg(end_angle)) < $(rad2deg(start_angle))"
        end_angle += 2pi
        @debug "(Sector)  Adjusted end_angle to be greater than start_angle: $(rad2deg(end_angle)) > $(rad2deg(start_angle))"
    end
    while end_angle - start_angle > pi
        @debug "(Sector) overshoot! end_angle - start_angle > π: $(rad2deg(end_angle - start_angle)) > 180"
        end_angle -= 2pi
        @debug "(Sector) Adjusted end_angle to be less than start_angle: $(rad2deg(end_angle)) < $(rad2deg(start_angle))"
    end
    angle_range = range(start_angle, stop=end_angle, length=num_points)
    @debug "(Sector) ______________________________ ∠ $(rad2deg(start_angle-end_angle))."
    return [Point2f(center[1] + radius * cos(a), center[2] + radius * sin(a)) for a in angle_range]
end

function _calculate_sector_polygon_points(params; num_arc_points=30) # increase `num_arc_points` for higher accuracy
    geom = _calculate_sector_geometry(params)
    nodes, centers = geom.Nodes, geom.Centers

    poly_points = Point2f[]
    get_angle(p1, p2) = atan((p1[2] - p2[2]), (p1[1] - p2[1]))
    # Start at F, go to A (Base)
    push!(poly_points, nodes.F)
    if params.r_corner > 1e-9
        start_angle = get_angle(nodes.F, centers.Base)
        @debug "Arc from F to A: start_angle=$(rad2deg(start_angle))"
        end_angle = get_angle(nodes.A, centers.Base)
        @debug "Arc from F to A: end_angle=$(rad2deg(end_angle))"
        append!(poly_points,
            _generate_arc_points(centers.Base, params.r_corner, start_angle, end_angle, num_arc_points)[2:end])
    else
        push!(poly_points, Point2f(0, params.r_back - params.d_sector), nodes.A)
    end

    # Line A to B
    push!(poly_points, nodes.B)

    # Arc B to C (Right Side)
    if params.r_corner > 1e-9
        start_angle = get_angle(nodes.B, centers.RightSide)
        @debug "Arc from B to C: start_angle=$(rad2deg(start_angle))"
        end_angle = get_angle(nodes.C, centers.RightSide)
        @debug "Arc from B to C: end_angle=$(rad2deg(end_angle))"
        append!(poly_points, _generate_arc_points(centers.RightSide, params.r_corner, start_angle, end_angle, num_arc_points)[2:end])
    else
        push!(poly_points, nodes.C)
    end

    # Arc C to D (Back)
    start_angle = get_angle(nodes.C, centers.Back)
    @debug "Arc from C to D: start_angle=$(rad2deg(start_angle))"
    end_angle = get_angle(nodes.D, centers.Back)
    @debug "Arc from C to D: end_angle=$(rad2deg(end_angle))"
    append!(poly_points, _generate_arc_points(centers.Back, params.r_back, start_angle, end_angle, num_arc_points)[2:end])

    # Arc D to E (Left Side)
    if params.r_corner > 1e-9
        start_angle = get_angle(nodes.D, centers.LeftSide)
        @debug "Arc from D to E: start_angle=$(rad2deg(start_angle))"
        end_angle = get_angle(nodes.E, centers.LeftSide)
        @debug "Arc from D to E: end_angle=$(rad2deg(end_angle))"
        append!(poly_points, _generate_arc_points(centers.LeftSide, params.r_corner, start_angle, end_angle, num_arc_points)[2:end])
    else
        push!(poly_points, nodes.E)
    end

    return poly_points
end

function _rotate_point(p::Point2f, angle_rad::Real)
    cos_a, sin_a = cos(angle_rad), sin(angle_rad)
    return Point2f(p[1] * cos_a - p[2] * sin_a, p[1] * sin_a + p[2] * cos_a)
end


function _shoelace_area(vertices::AbstractVector{Point{2,T}}) where {T<:Real}
    n::Int = length(vertices)
    if n < 3
        @warn "Polygon must have at least 3 vertices to compute area. Returning 0 area"
        return zero(T)
    end
    area::T = zero(T)
    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[mod1(i + 1, n)]
        area += (p1[1] * p2[2] - p2[1] * p1[2])
    end
    return abs(area) / T(2)
end

function _calculate_polygon_centroid(vertices::AbstractVector{Point{2,T}}) where {T<:Real}
    n = length(vertices)
    if n < 3
        return Point2f(0, 0)
    end

    
    Cx = zero(T)
    Cy = zero(T)
    signed_area = zero(T)

    for i in 1:n
        p1 = vertices[i]
        p2 = vertices[mod1(i + 1, n)]
        cross_prod = (p1[1] * p2[2] - p2[1] * p1[2])
        signed_area += (p1[1] * p2[2] - p2[1] * p1[2])
        Cx += (p1[1] + p2[1]) * cross_prod
        Cy += (p1[2] + p2[2]) * cross_prod
    end
    signed_area /= 2

    if abs(signed_area) < 1e-12
        return Point2f(0, 0)
    end

    Cx /= (6 * signed_area)
    Cy /= (6 * signed_area)

    return Point2f(Cx, Cy)
end