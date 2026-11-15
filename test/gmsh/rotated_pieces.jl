using Gmsh: gmsh
using LineCableModels.BaseGeometry: Arc, Contour, dimension, draw
using LineCableModels.Commons: ComputationOptions
using LineCableModels.Engine: GmshBackend

include("../datamodel/shape_gallery.jl")

# Cut `shape` and the contacts of its circles in a fresh geometry of `model`, as `build`
# does. Draw it by `draw!`. Return each curve, its registered ends and its chord, in the
# order drawn.
function drawn_pieces(draw!, model, options, shape)
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    gmsh.model.add("pieces")
    geometry = extension.Geometry(model, options)
    extension.cut!(geometry, shape)
    extension.cut!(geometry)
    tags = draw!(geometry, shape)
    ends = [(geometry.positions[p], geometry.positions[q])
            for (p, q) in (geometry.ends[abs(tag)] for tag in tags)]
    pieces = (curves = [geometry.curves[abs(tag)] for tag in tags], ends = ends,
        chords = [hypot((p .- q)...) for (p, q) in ends])
    gmsh.model.remove()
    return pieces
end

# Draw the outer path and then each hole of the contour of `shape`, as `faces!` does.
function contour_tags(geometry, shape)
    contour = Contour(shape)
    return reduce(vcat, [draw(geometry, path).curves
                         for path in (contour.outer, contour.holes...)])
end

# Draw the strips of `shape` at `depth`, and return the outer side, the inner side and the
# spokes of each strip.
function strip_tags(geometry, shape, depth)
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    extension.strips!(geometry, shape, 1, depth)
    return reduce(vcat, [vcat(strip.curves..., collect(strip.spokes))
                         for strip in geometry.strips])
end

# Draw the faces of `shape`, and return every curve that they drew, in drawing order.
function face_tags(geometry, shape)
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    extension.faces!(geometry, shape)
    return sort!(collect(keys(geometry.curves)))
end

# Compare the pieces with those of the copy `turned` by `at`. The point copies differ
# by at most `point` in each coordinate. Registration rounds an end to 15 significant
# digits, half of 1e-16 m below 0.1 m. Rotation by `at` adds 4eps(1/16).
#
# The chords or radii of the copies differ by at most 2√2(1e-16 m + point).
# An end moved by √2 point turns an arc of radius r by √2 point / r.
# Its break angles include at most eight roundings below 4π.
function test_same_pieces(pieces, turned, at, point)
    lengths = 2sqrt(2) * (1.0e-16 + point)
    @test length(pieces.curves) == length(turned.curves)
    @test typeof.(pieces.curves) == typeof.(turned.curves)
    @test maximum(
        maximum(abs, at(p) .- q)
        for (ends, turned_ends) in zip(pieces.ends, turned.ends)
        for (p, q) in zip(ends, turned_ends);
        init = 0.0) <= 1.0e-16 + point + 4eps(1 / 16)
    @test maximum(abs(a - b) for (a, b) in zip(pieces.chords, turned.chords); init = 0.0) <=
          lengths
    arcs = [curve for curve in pieces.curves if curve isa Arc]
    turned_arcs = [curve for curve in turned.curves if curve isa Arc]
    @test maximum(abs(a.r - b.r) for (a, b) in zip(arcs, turned_arcs); init = 0.0) <=
          lengths
    @test all(abs(abs(a.span) - abs(b.span)) <= 2sqrt(2) * point / a.r + 8eps(4π)
    for (a, b) in zip(arcs, turned_arcs))
    # The extension matches a break within 1e-12 rad.
    @test all(abs(arc.span) <= π / 2 + 1.0e-12 for arc in [arcs; turned_arcs])
end

@testset "a rotated copy of a shape is cut into the same pieces" begin
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    model = extension.FieldModel(last(first(REFERENCES)), FEM())
    options = ComputationOptions(GmshBackend)
    gallery = shape_gallery()
    # A placement rounds each coordinate in two products, two sums, a sine and a cosine.
    # Each magnitude is below 1/16 m, with error at most 4eps(1/16). Either order places
    # each point twice and rounds the sum of two rotations below 8 rad. That angular
    # rounding moves a point within 0.01 m of its center by less than 4eps(1/16).
    # The copies of each point differ by at most 24eps(1/16).
    point = 24eps(1 / 16)
    angles = (0.3, 2.0, -1.1)
    # A strip lies at 0.45 of the distance across, the depth of a strip section. The wires
    # of a course touch the inner circle of its annulus: it has no distance across.
    stripped = filter(shape -> dimension(shape) > 0, gallery)
    perforated = filter(shape -> shape isa PerforatedShape, gallery)
    gmsh.initialize()
    try
        @testset "$(nameof(typeof(shape))) contour, $α rad" for shape in gallery,
            α in angles

            at = Placement(0.0, 0.0, α)
            turned = LineCableModels.resolve(at, shape)
            test_same_pieces(drawn_pieces(contour_tags, model, options, shape),
                drawn_pieces(contour_tags, model, options, turned), at, point)
        end
        @testset "$(nameof(typeof(shape))) strips, $α rad" for shape in stripped,
            α in angles

            at = Placement(0.0, 0.0, α)
            turned = LineCableModels.resolve(at, shape)
            depth = 0.45 * dimension(shape)
            test_same_pieces(
                drawn_pieces((geometry, member) -> strip_tags(geometry, member, depth),
                    model, options, shape),
                drawn_pieces((geometry, member) -> strip_tags(geometry, member, depth),
                    model, options, turned),
                at, point)
        end
        @testset "$(nameof(typeof(shape))) faces, $α rad" for shape in perforated,
            α in angles

            at = Placement(0.0, 0.0, α)
            turned = LineCableModels.resolve(at, shape)
            test_same_pieces(drawn_pieces(face_tags, model, options, shape),
                drawn_pieces(face_tags, model, options, turned), at, point)
        end
    finally
        gmsh.finalize()
    end
end
