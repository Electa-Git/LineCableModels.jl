struct PreviewPolygon{G, C, S}
    geometry::G
    label::Union{Nothing, String}
    group::Symbol
    color::C
    stroke::S
    width::Float64
end

struct PreviewReference{V, C}
    values::V
    group::Symbol
    color::C
    width::Float64
end
