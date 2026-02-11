"""
Primitive drawing functions for the FEMTools.jl module.
These functions handle the creation of geometric entities in Gmsh.
"""

function add_mesh_points(;
	radius_ext::Number,
	theta_0::Number,
	theta_1::Number,
	mesh_size::Number,
	radius_in::Number = 0.0,
	num_points_ang::Integer = 8,
	num_points_rad::Integer = 0,
	C::Tuple{Number, Number} = (0.0, 0.0),
	theta_offset::Number = 0.0)

	point_tags = Vector{Int}()
	center_x, center_y = C

	# Handle special cases
	if num_points_ang <= 0 && num_points_rad <= 0
		# Single point at center C
		point_tag = gmsh.model.occ.add_point(center_x, center_y, 0.0, mesh_size)
		gmsh.model.set_entity_name(
			0,
			point_tag,
			"mesh_size_$(round(mesh_size, sigdigits=6))",
		)
		return [point_tag]
	end

	# Circular arc (default case or when num_points_rad=0)
	if num_points_rad == 0
		r = radius_ext  # Use external radius as default
		np_ang = max(2, num_points_ang)  # At least 2 points for an arc

		for i in 0:(np_ang-1)
			t_ang = i / (np_ang - 1)
			theta = theta_0 + t_ang * (theta_1 - theta_0) + theta_offset

			x = center_x + r * cos(theta)
			y = center_y + r * sin(theta)

			point_tag = gmsh.model.occ.add_point(x, y, 0.0, mesh_size)
			gmsh.model.set_entity_name(
				0,
				point_tag,
				"mesh_size_$(round(mesh_size, sigdigits=6))",
			)
			push!(point_tags, point_tag)
		end

		return point_tags
	end

	# Radial line (when theta_0 == theta_1)
	if theta_0 == theta_1
		theta = theta_0 + theta_offset
		np_rad = max(2, num_points_rad)

		for j in 0:(np_rad-1)
			t_rad = j / (np_rad - 1)
			r = radius_in + t_rad * (radius_ext - radius_in)

			x = center_x + r * cos(theta)
			y = center_y + r * sin(theta)

			point_tag = gmsh.model.occ.add_point(x, y, 0.0, mesh_size)
			gmsh.model.set_entity_name(
				0,
				point_tag,
				"mesh_size_$(round(mesh_size, sigdigits=6))",
			)
			push!(point_tags, point_tag)
		end

		return point_tags
	end

	# 2D array of points (both radial and angular)
	np_rad = max(2, num_points_rad)
	np_ang = max(2, num_points_ang)

	for j in 0:(np_rad-1)
		t_rad = j / (np_rad - 1)
		r = radius_in + t_rad * (radius_ext - radius_in)

		for i in 0:(np_ang-1)
			t_ang = i / (np_ang - 1)
			theta = theta_0 + t_ang * (theta_1 - theta_0) + theta_offset
			*
			x = center_x + r * cos(theta)
			y = center_y + r * sin(theta)

			point_tag = gmsh.model.occ.add_point(x, y, 0.0, mesh_size)
			gmsh.model.set_entity_name(
				0,
				point_tag,
				"mesh_size_$(round(mesh_size, sigdigits=6))",
			)
			push!(point_tags, point_tag)
		end
	end

	return point_tags
end


"""
$(TYPEDSIGNATURES)

Draw a point with specified coordinates and mesh size.

# Arguments

- `x`: X-coordinate \\[m\\].
- `y`: Y-coordinate \\[m\\].
- `z`: Z-coordinate \\[m\\].

# Returns

- Gmsh point tag \\[dimensionless\\].

# Examples

```julia
point_tag = $(FUNCTIONNAME)(0.0, 0.0, 0.0, 0.01)
```
"""
function draw_point(x::Number, y::Number, z::Number)
	return gmsh.model.occ.add_point(x, y, z)
end

"""
$(TYPEDSIGNATURES)

Draw a line between two points.

# Arguments

- `x1`: X-coordinate of the first point \\[m\\].
- `y1`: Y-coordinate of the first point \\[m\\].
- `x2`: X-coordinate of the second point \\[m\\].
- `y2`: Y-coordinate of the second point \\[m\\].

# Returns

- Gmsh line tag \\[dimensionless\\].

# Examples

```julia
line_tag = $(FUNCTIONNAME)(0.0, 0.0, 1.0, 0.0, 0.01)
```
"""
function draw_line(
	x1::Number,
	y1::Number,
	x2::Number,
	y2::Number,
	mesh_size::Number,
	num_points::Number,
)

	# Calculate line parameters
	line_length = sqrt((x2 - x1)^2 + (y2 - y1)^2)
	x_center = (x1 + x2) / 2
	y_center = (y1 + y2) / 2

	# Calculate angle in polar coordinates (in radians)
	theta = atan(y2 - y1, x2 - x1)

	# Use the distance as a "domain radius" for placing mesh points
	radius = line_length / 2

	# Create a unique marker for this line
	marker = [x_center, y_center, 0.0]  # Center of the line

	marker_tag = gmsh.model.occ.add_point(marker[1], marker[2], marker[3], mesh_size)
	gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_size, sigdigits=6))")

	mesh_points = add_mesh_points(
		radius_in = -radius,
		radius_ext = radius,
		theta_0 = theta,
		theta_1 = theta,
		mesh_size = mesh_size,
		num_points_ang = 0,
		num_points_rad = num_points, # Not strictly a circumference, but the trick works
		C = (x_center, y_center),
	)

	tag = gmsh.model.occ.add_line(mesh_points[1], mesh_points[end])

	# Add midpoint markers between each pair of mesh points
	segment_markers = Vector{Vector{Float64}}()
	push!(segment_markers, marker)

	if length(mesh_points) >= 2
		# Iterate through adjacent pairs of mesh points
		for i in 1:(num_points-1)
			t = (i - 0.5) / (num_points - 1)  # Parametric coordinate (0.5 between points)
			mid_x = x1 + t * (x2 - x1)
			mid_y = y1 + t * (y2 - y1)
			# Create marker at midpoint
			mid_marker = [mid_x, mid_y, 0.0]
			push!(segment_markers, mid_marker)
		end
	end

	return tag, mesh_points, segment_markers
end

"""
$(TYPEDSIGNATURES)

Draw a circular disk with specified center and radius.

# Arguments

- `x`: X-coordinate of the center \\[m\\].
- `y`: Y-coordinate of the center \\[m\\].
- `radius`: Radius of the disk \\[m\\].

# Returns

- Gmsh surface tag \\[dimensionless\\].

# Examples

```julia
disk_tag = $(FUNCTIONNAME)(0.0, 0.0, 0.5, 0.01)
```
"""
function draw_disk(
	x::Number,
	y::Number,
	radius::Number,
	mesh_size::Number,
	num_points::Number,
)

	tag = gmsh.model.occ.add_disk(x, y, 0.0, radius, radius)

	mesh_points = add_mesh_points(
		radius_in = radius,
		radius_ext = radius,
		theta_0 = 0,
		theta_1 = 2 * pi,
		mesh_size = mesh_size,
		num_points_ang = num_points,
		C = (x, y),
		theta_offset = 0, #pi / 15
	)


	marker = [x, y + 0.99 * radius, 0.0] # A very small offset inwards the circle
	marker_tag = gmsh.model.occ.add_point(marker[1], marker[2], marker[3], mesh_size)
	gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_size, sigdigits=6))")

	# Add midpoint markers between each pair of mesh points
	arc_markers = Vector{Vector{Float64}}()

	if num_points >= 2
		# Calculate the angular step between mesh points
		theta_step = 2 * pi / num_points

		# Add a midpoint marker for each arc segment
		for i in 1:num_points
			# Calculate midpoint theta (angle)
			theta_mid = (i - 0.5) * theta_step

			# Calculate midpoint coordinates
			mid_x = x + radius * cos(theta_mid)
			mid_y = y + radius * sin(theta_mid)

			# Create marker at the midpoint
			mid_marker = [mid_x, mid_y, 0.0]

			push!(arc_markers, mid_marker)
		end
	end

	return tag, mesh_points, marker, arc_markers
end


"""
$(TYPEDSIGNATURES)

Draw an annular (ring) shape with specified center, inner radius, and outer radius.

# Arguments

- `x`: X-coordinate of the center \\[m\\].
- `y`: Y-coordinate of the center \\[m\\].
- `radius_in`: Inner radius of the annular shape \\[m\\].
- `radius_ext`: Outer radius of the annular shape \\[m\\].

# Returns

- Gmsh surface tag \\[dimensionless\\].

# Examples

```julia
annular_tag = $(FUNCTIONNAME)(0.0, 0.0, 0.3, 0.5, 0.01)
```
"""
function draw_annular(
	x::Number,
	y::Number,
	radius_in::Number,
	radius_ext::Number,
	mesh_size::Number,
	num_points::Number;
	inner_points::Bool = false,
)
	# Create outer disk
	outer_disk = gmsh.model.occ.add_disk(x, y, 0.0, radius_ext, radius_ext)

	# Create inner disk
	inner_disk = gmsh.model.occ.add_disk(x, y, 0.0, radius_in, radius_in)

	# Cut inner disk from outer disk to create annular shape
	annular_obj, _ = gmsh.model.occ.cut([(2, outer_disk)], [(2, inner_disk)])

	# Return the tag of the resulting surface
	if length(annular_obj) > 0
		tag = annular_obj[1][2]
	else
		Base.error("Failed to create annular shape.")
	end

	mesh_points = add_mesh_points(
		radius_in = radius_ext,
		radius_ext = radius_ext,
		theta_0 = 0,
		theta_1 = 2 * pi,
		mesh_size = mesh_size,
		num_points_ang = num_points,
		C = (x, y),
		theta_offset = 0, #pi / 15
	)

	if inner_points
		mesh_points = add_mesh_points(
			radius_in = radius_in,
			radius_ext = radius_in,
			theta_0 = 0,
			theta_1 = 2 * pi,
			mesh_size = mesh_size,
			num_points_ang = num_points,
			C = (x, y),
			theta_offset = pi / 3,
		)
	end

	marker = [x, y + (radius_in + 0.99 * (radius_ext - radius_in)), 0.0]
	marker_tag = gmsh.model.occ.add_point(marker[1], marker[2], marker[3], mesh_size)
	gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_size, sigdigits=6))")

	# Add midpoint markers between each pair of mesh points
	arc_markers = Vector{Vector{Float64}}()

	if num_points >= 2
		# Calculate the angular step between mesh points
		theta_step = 2 * pi / num_points

		# Add a midpoint marker for each arc segment
		for i in 1:num_points
			# Calculate midpoint theta (angle)
			theta_mid = (i - 0.5) * theta_step

			# Calculate midpoint coordinates
			mid_x = x + radius_ext * cos(theta_mid)
			mid_y = y + radius_ext * sin(theta_mid)

			# Create marker at the midpoint
			mid_marker = [mid_x, mid_y, 0.0]

			push!(arc_markers, mid_marker)
		end
	end

	return tag, mesh_points, marker, arc_markers
end


"""
$(TYPEDSIGNATURES)

Draw a rectangle with specified center, width, and height.

# Arguments

- `x`: X-coordinate of the center \\[m\\].
- `y`: Y-coordinate of the center \\[m\\].
- `width`: Width of the rectangle \\[m\\].
- `height`: Height of the rectangle \\[m\\].

# Returns

- Gmsh surface tag \\[dimensionless\\].

# Examples

```julia
rect_tag = $(FUNCTIONNAME)(0.0, 0.0, 1.0, 0.5, 0.01)
```
"""
function draw_rectangle(x::Number, y::Number, width::Number, height::Number)
	# Calculate corner coordinates
	x1 = x - width / 2
	y1 = y - height / 2
	x2 = x + width / 2
	y2 = y + height / 2

	# Create rectangle
	return gmsh.model.occ.add_rectangle(x1, y1, 0.0, width, height)
end

"""
$(TYPEDSIGNATURES)

Draw a circular arc between two points with a specified center.

# Arguments

- `x1`: X-coordinate of the first point \\[m\\].
- `y1`: Y-coordinate of the first point \\[m\\].
- `x2`: X-coordinate of the second point \\[m\\].
- `y2`: Y-coordinate of the second point \\[m\\].
- `xc`: X-coordinate of the center \\[m\\].
- `yc`: Y-coordinate of the center \\[m\\].

# Returns

- Gmsh curve tag \\[dimensionless\\].

# Examples

```julia
arc_tag = $(FUNCTIONNAME)(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.01)
```
"""
function draw_arc(x1::Number, y1::Number, x2::Number, y2::Number, xc::Number, yc::Number)
	p1 = gmsh.model.occ.add_point(x1, y1, 0.0)
	p2 = gmsh.model.occ.add_point(x2, y2, 0.0)
	pc = gmsh.model.occ.add_point(xc, yc, 0.0)

	return gmsh.model.occ.add_circle_arc(p1, pc, p2)
end

"""
$(TYPEDSIGNATURES)

Draw a circle with specified center and radius.

# Arguments

- `x`: X-coordinate of the center \\[m\\].
- `y`: Y-coordinate of the center \\[m\\].
- `radius`: Radius of the circle \\[m\\].

# Returns

- Gmsh curve tag \\[dimensionless\\].

# Examples

```julia
circle_tag = $(FUNCTIONNAME)(0.0, 0.0, 0.5, 0.01)
```
"""
function draw_circle(x::Number, y::Number, radius::Number)
	return gmsh.model.occ.add_circle(x, y, 0.0, radius)
end

"""
$(TYPEDSIGNATURES)

Draw a polygon with specified vertices.

# Arguments

- `vertices`: Array of (x,y) coordinates for the vertices \\[m\\].

# Returns

- Gmsh surface tag \\[dimensionless\\].

# Examples

```julia
vertices = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
polygon_tag = $(FUNCTIONNAME)(vertices, 0.01)
```
"""
function draw_polygon(vertices::Vector{<:Tuple{<:Number, <:Number}})
	# Create points
	points = Vector{Int}()
	for (x, y) in vertices
		push!(points, gmsh.model.occ.add_point(x, y, 0.0))
	end

	# Create lines
	lines = Vector{Int}()
	for i in 1:length(points)
		next_i = i % length(points) + 1
		push!(lines, gmsh.model.occ.add_line(points[i], points[next_i]))
	end

	# Create curve loop
	curve_loop = gmsh.model.occ.add_curve_loop(lines)

	# Create surface
	return gmsh.model.occ.add_plane_surface([curve_loop])
end

function draw_transition_region(
	x::Number,
	y::Number,
	radii::Vector{<:Number},
	mesh_sizes::Vector{<:Number},
	num_points::Number,
)
	# Validate inputs
	if length(radii) != length(mesh_sizes)
		Base.error("Radii and mesh_sizes vectors must have the same length")
	end

	n_regions = length(radii)
	if n_regions < 1
		Base.error("At least one radius must be provided")
	end

	# Sort radii in ascending order if not already sorted
	if !issorted(radii)
		p = sortperm(radii)
		radii = radii[p]
		mesh_sizes = mesh_sizes[p]
	end

	# Note to future self: the reason why I did this is because in the edge case when the bounding box coincides with the cable outermost radius (i.e. when you want the transition region around 1 single cable and not several), the earth marker ends up inside the cable region, which present me does not need to explain to future me why it's bad.

	rad_buffer = 0.001
	radii[1] += rad_buffer  # Ensure the innermost radius is slightly larger than zero to avoid ambiguous regions
	tags = Int[]
	all_mesh_points = Int[]
	markers = Vector{Vector{Float64}}()

	# Create all disks
	disk_tags = Int[]
	for i in 1:n_regions
		disk_tag = gmsh.model.occ.add_disk(x, y, 0.0, radii[i], radii[i])
		gmsh.model.occ.synchronize()
		push!(disk_tags, disk_tag)
	end

	# Add the innermost disk to output
	push!(tags, disk_tags[1])

	# Add mesh points for innermost disk
	inner_mesh_points = add_mesh_points(
		radius_in = radii[1],
		radius_ext = radii[1],
		theta_0 = 0,
		theta_1 = 2 * pi,
		mesh_size = mesh_sizes[1],
		num_points_ang = num_points,
		C = (x, y),
		theta_offset = 0,
	)
	append!(all_mesh_points, inner_mesh_points)

	# Create marker at the midpoint between the real radius and the added buffer
	# Feel free to implement a less dumb way to do this
	radius_inner_marker = (radii[1] + radii[1] - rad_buffer) / 2
	inner_marker = [x, y + radius_inner_marker, 0.0]

	marker_tag = gmsh.model.occ.add_point(
		inner_marker[1],
		inner_marker[2],
		inner_marker[3],
		mesh_sizes[1],
	)
	gmsh.model.set_entity_name(0, marker_tag, "marker_$(round(mesh_sizes[1], sigdigits=6))")
	push!(markers, inner_marker)

	# Synchronize the model
	gmsh.model.occ.synchronize()

	# Create annular regions for the rest
	for i in 2:n_regions
		# Cut the inner disk from the outer disk
		annular_obj, _ =
			gmsh.model.occ.cut([(2, disk_tags[i])], [(2, disk_tags[i-1])], false, false)

		# Get the resulting surface tag
		if length(annular_obj) > 0
			annular_tag = annular_obj[1][2]
			push!(tags, annular_tag)

			# Add mesh points on the boundary
			boundary_points = add_mesh_points(
				radius_in = radii[i],
				radius_ext = radii[i],
				theta_0 = 0,
				theta_1 = 2 * pi,
				mesh_size = mesh_sizes[i],
				num_points_ang = num_points,
				C = (x, y),
				theta_offset = 0,
			)
			append!(all_mesh_points, boundary_points)

			# Create marker at 99% of the way from inner to outer radius
			radius_marker = radii[i-1] + 0.99 * (radii[i] - radii[i-1])
			annular_marker = [x, y + radius_marker, 0.0]
			marker_tag = gmsh.model.occ.add_point(
				annular_marker[1],
				annular_marker[2],
				annular_marker[3],
				mesh_sizes[i],
			)
			gmsh.model.set_entity_name(
				0,
				marker_tag,
				"marker_$(round(mesh_sizes[i], sigdigits=6))",
			)
			push!(markers, annular_marker)
		else
			Base.error(
				"Failed to create annular region for radii $(radii[i-1]) and $(radii[i])",
			)
		end
	end

	return tags, all_mesh_points, markers
end

function draw_polygon(vertices::Vector{<:Point})
	@debug "Drawing a polygon for a Point vertices"
    # Create points
    points = [gmsh.model.occ.add_point(v[1], v[2], 0.0) for v in vertices]

    # Create lines
    lines = [gmsh.model.occ.add_line(points[i], points[i % length(points) + 1]) for i in 1:length(points)]

    # Create curve loop and surface
    curve_loop = gmsh.model.occ.add_curve_loop(lines)
    surface_tag = gmsh.model.occ.add_plane_surface([curve_loop])
    
    # Synchronize to make the new entity available for calculations
    gmsh.model.occ.synchronize()

    # calculate the centroid of the vertices as a marker
    if isempty(vertices)
        error("Cannot calculate centroid of empty vertex list.")
    end
    avg_x = sum(v[1] for v in vertices) / length(vertices)
    avg_y = sum(v[2] for v in vertices) / length(vertices)
    marker = [avg_x, avg_y, 0.0]

    return surface_tag, marker
end

"""
$(TYPEDSIGNATURES)

Draw a polygon with a hole.

# Arguments
- `outer_vertices`: A vector of (x,y) coordinates for the outer boundary.
- `inner_vertices`: A vector of (x,y) coordinates for the inner boundary (the hole).

# Returns
- A tuple containing the Gmsh surface tag and a marker point `[x, y, z]`.
"""
function _densify_vertices(vertices::Vector{<:Point}, max_len::Number)
    new_vertices = Point[]
    if isempty(vertices)
        return new_vertices
    end
    for i in 1:length(vertices)
        p1 = vertices[i]
        p2 = vertices[i % length(vertices) + 1]
        
        push!(new_vertices, p1)
        
        edge_vec = p2 - p1
        edge_len = norm(edge_vec)
        
        if edge_len > max_len
            num_segments = ceil(Int, edge_len / max_len)
            for j in 1:(num_segments - 1)
                intermediate_point = p1 + (j / num_segments) * edge_vec
                push!(new_vertices, intermediate_point)
            end
        end
    end
    return new_vertices
end

function draw_polygon_with_hole(outer_vertices::Vector{<:Point}, inner_vertices::Vector{<:Point}, max_edge_length::Number)
    
    new_outer_vertices = _densify_vertices(outer_vertices, max_edge_length)
    new_inner_vertices = _densify_vertices(inner_vertices, max_edge_length)

    # Create outer boundary
    outer_points = [gmsh.model.occ.add_point(v[1], v[2], 0.0) for v in new_outer_vertices]
    outer_lines = [gmsh.model.occ.add_line(outer_points[i], outer_points[i % length(outer_points) + 1]) for i in 1:length(outer_points)]
    outer_loop = gmsh.model.occ.add_curve_loop(outer_lines)

    # Create inner boundary (hole)
    inner_points = [gmsh.model.occ.add_point(v[1], v[2], 0.0) for v in new_inner_vertices]
    inner_lines = [gmsh.model.occ.add_line(inner_points[i], inner_points[i % length(inner_points) + 1]) for i in 1:length(inner_points)]
    inner_loop = gmsh.model.occ.add_curve_loop(inner_lines)

    # Create surface with hole
    surface_tag = gmsh.model.occ.add_plane_surface([outer_loop, inner_loop])
    
    # Synchronize to make the new entity available for calculations
    gmsh.model.occ.synchronize()

    # A marker point must be inside the insulator, but outside the conductor.
    # A point halfway between the inner and outer boundaries along one of the vertices should work.
    marker_point = (outer_vertices[1] + inner_vertices[1]) / 2
    marker = [marker_point[1], marker_point[2], 0.0]

    return surface_tag, marker
end

function get_system_centroid(cable_system::LineCableSystem, cable_idx::Vector{<:Integer})
	# Check if cable_idx is empty
	if isempty(cable_idx)
		Base.error("Cable index vector cannot be empty")
	end

	# Check if any index is out of bounds
	if any(idx -> idx < 1 || idx > length(cable_system.cables), cable_idx)
		Base.error("Cable index out of bounds")
	end

	# Extract coordinates
	horz_coords = [cable_system.cables[idx].horz for idx in cable_idx]
	vert_coords = [cable_system.cables[idx].vert for idx in cable_idx]

	# Calculate centroid
	centroid_x = sum(horz_coords) / length(horz_coords)
	centroid_y = sum(vert_coords) / length(vert_coords)

	# Find the maximum distance from centroid to any cable's edge
	max_distance = 0.0
	characteristic_len = Inf

	for idx in cable_idx
		cable_position = cable_system.cables[idx]

		# Calculate distance from centroid to cable center
		distance_to_center = sqrt(
			(cable_position.horz - centroid_x)^2 + (cable_position.vert - centroid_y)^2,
		)

		# Get the outermost component (last component in the vector)
		if !isempty(cable_position.design_data.components)
			last_component = cable_position.design_data.components[end]

			outer_radius = last_component.insulator_group.radius_ext

			insulator_radius_in = last_component.insulator_group.layers[end].radius_in
			last_layer_thickness = outer_radius - insulator_radius_in


			# Add cable radius to get distance to edge
			total_distance = distance_to_center + outer_radius
			max_distance = max(max_distance, total_distance)
			characteristic_len = min(characteristic_len, last_layer_thickness)
		end
	end

	return (centroid_x, centroid_y, max_distance, characteristic_len)
end


"""
$(TYPEDSIGNATURES)

Calculate the coordinates of air gaps in a wire array.

# Arguments

- `num_wires`: Number of wires in the array \\[dimensionless\\].
- `radius_wire`: Radius of each wire \\[m\\].
- `radius_in`: Inner radius of the wire array \\[m\\].

# Returns

- Vector of marker positions (3D coordinates) for air gaps \\[m\\].

# Notes

This function calculates positions for markers that are guaranteed to be in the air gaps
between wires in a wire array. These markers are used to identify the air regions after
boolean fragmentation operations.

# Examples

```julia
markers = $(FUNCTIONNAME)(7, 0.002, 0.01)
```
"""
function get_air_gap_markers(num_wires::Int, radius_wire::Number, radius_in::Number)
	markers = Vector{Vector{Float64}}()

	lay_radius = radius_in + radius_wire

	num_angular_markers = num_wires == 1 ? 6 : num_wires
	# For multiple wires, place markers between adjacent wires
	angle_step = 2π / num_angular_markers
	for i in 0:(num_angular_markers-1)
		angle = i * angle_step + (angle_step / 2)  # Midway between wires
		r = lay_radius + (radius_wire / 2)  # Slightly outward
		x = r * cos(angle)
		y = r * sin(angle)
		push!(markers, [x, y, 0.0])
	end
	return markers
end

function draw_polygon(vertices::Vector{<:Point}, max_edge_length::Number)
    if isempty(vertices)
        error("Cannot draw a polygon with no vertices.")
    end

    new_vertices = _densify_vertices(vertices, max_edge_length)

    # Create points
    points = [gmsh.model.occ.add_point(v[1], v[2], 0.0) for v in new_vertices]

    # Create lines
    lines = [gmsh.model.occ.add_line(points[i], points[i % length(points) + 1]) for i in 1:length(points)]

    # Create curve loop and surface
    curve_loop = gmsh.model.occ.add_curve_loop(lines)
    surface_tag = gmsh.model.occ.add_plane_surface([curve_loop])
    
    # Synchronize to make the new entity available for calculations
    gmsh.model.occ.synchronize()

    # calculate the centroid of the original vertices as a marker
    if isempty(vertices)
        error("Cannot calculate centroid of empty vertex list.")
    end
    avg_x = sum(v[1] for v in vertices) / length(vertices)
    avg_y = sum(v[2] for v in vertices) / length(vertices)
    marker = [avg_x, avg_y, 0.0]

    return surface_tag, marker
end


