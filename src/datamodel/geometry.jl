"""
Return the outer radius of a materialized cable design.
"""
function outer_radius(design::CableDesign)
    return support(boundary(design.geometry))
end

