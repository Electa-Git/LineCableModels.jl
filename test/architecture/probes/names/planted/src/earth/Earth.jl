module Earth
import ..Commons
struct Ring end
struct Layer
    resistivity::Float64
    mean_free_path::Float64
end
struct Coefficients
    a::Float64
    ab::Float64
    bc::Float64
    c::Float64
end
struct BuildPlan end
struct AbstractDisk end
struct Formula <: Commons.AbstractFormulation end
struct Plan end
earth_layer(x) = x
ring_poses(ring::Ring) = ring
radial_half_extent(x) = x
mesh_build(x) = x
get_radius(x) = x
frobnicate(x) = x
abc(x) = x
spread(x; outer_layer_count = 1) = x
length(::Layer) = 1
compute(::Layer) = 1
end
