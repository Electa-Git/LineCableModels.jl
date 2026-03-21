# Global accessors
@inline r_in(s::AbstractShape) = s.r_in
@inline r_ex(s::AbstractShape) = s.r_ex



include("solidcore.jl")
include("tubular.jl")
# include("enclosure.jl")
# include("wires.jl")
# include("helical.jl")
# include("stranded.jl")

# ---------------------------------------------------------
# The Fuzzy Characteristic Length Trait
# ---------------------------------------------------------
# Returns the characteristic dimension of the primitive. 
# If someone writes a Rectangle primitive and defines char_len as the diagonal, 
# the stacking engine will blindly build overlapping garbage.
@inline char_len(s::SolidCore) = 2 * r_ex(s)
# @inline char_len(s::TubularLayer) = r_ex(s) - r_in(s)
# @inline char_len(s::CircularWire) = 2 * s.r
# @inline char_len(s::RectangularWire) = s.h

# ==========================================
# THE TOPOLOGICAL FOOTPRINT
# Extracts the absolute boundary as a concrete primitive.
# ==========================================
@inline boundary(p::AbstractCablePart) = boundary(p.shape)

# A SolidCore's footprint is a circle at r_ex.
@inline boundary(s::SolidCore) = Circular(r_ex(s))

# A Tubular's footprint is a circle at r_ex.
@inline boundary(s::TubularLayer) = Circular(r_ex(s))

# A PartGroup's bounding footprint is currently a circumscribed circle.
@inline boundary(g::PartGroup) = Circular(r_ex(g))