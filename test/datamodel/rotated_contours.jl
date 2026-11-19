using InteractiveUtils: subtypes
using LineCableModels.BaseGeometry: Contour

include("shape_gallery.jl")

@testset "a rotated copy of a shape has the same contour" begin
    gallery = shape_gallery()
    @testset "the gallery holds every concrete shape type" begin
        family = Any[LineCableModels.AbstractShape]
        for type in family
            append!(family, subtypes(type))
        end
        @test Set(Base.typename(typeof(shape)).wrapper for shape in gallery) ==
              Set(filter(!isabstracttype, family))
    end
    # A placement rounds each coordinate in two products, two sums, a sine and a cosine.
    # Each magnitude is below 1/16 m, with error at most 4eps(1/16). Either order places
    # each point twice and rounds the sum of two rotations below 8 rad. That angular
    # rounding moves a point within 0.01 m of its center by less than 4eps(1/16).
    # The copies of each point differ by at most 24eps(1/16).
    bound = 24eps(1 / 16)
    @testset "$(nameof(typeof(shape))) turned by $α rad" for shape in gallery,
        α in (0.3, 2.0, -1.1)

        at = Placement(0.0, 0.0, α)
        placed = Contour(LineCableModels.resolve(at, shape))
        turned = LineCableModels.resolve(at, Contour(shape))
        paths = [placed.outer, placed.holes...]
        turned_paths = [turned.outer, turned.holes...]
        @test length.(paths) == length.(turned_paths)
        @test [typeof.(path) for path in paths] == [typeof.(path) for path in turned_paths]
        @test maximum(maximum(abs, curve(t) .- turned_curve(t))
        for (path, turned_path) in zip(paths, turned_paths)
        for (curve, turned_curve) in zip(path, turned_path) for t in (0, 1 / 2, 1)) <=
              bound
    end
end

# No check of the ends runs where a cross-section builds its contour. Each such method
# closes its paths, the sector's in each of its forms.
@testset "the contour of each shape is closed" begin
    at = Placement(0.012, -0.007, 0.4)
    sectors = [Sector(ri, 10.0e-3, span, apex, fillet, at)
               for (ri, span, apex, fillet) in (
        (0, 2π, 0, 0), (4.0e-3, 2π, 0, 0), (0, 2π / 3, 0, 0), (4.0e-3, 2π / 3, 0, 0),
        (0, 2π / 3, 1.0e-3, 0), (4.0e-3, 2π / 3, 1.0e-3, 0),
        (0, 2π / 3, 1.0e-3, 1.0e-3), (4.0e-3, 2π / 3, 1.0e-3, 1.0e-3))]
    @testset "$(nameof(typeof(shape)))" for shape in
                                            [shape_gallery(); sectors;
                                             Annulus(0, 3.0e-3, at)]
        contour = Contour(shape)
        @test LineCableModels.validate(contour) === contour
    end
end
