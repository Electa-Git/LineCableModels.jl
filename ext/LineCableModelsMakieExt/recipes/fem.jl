function _fem_planar(coordinates)
    all(isfinite, coordinates) || throw(ArgumentError("mesh coordinates must be finite"))
    z = view(coordinates, 3, :)
    isempty(z) && throw(ArgumentError("mesh has no coordinates"))
    all(==(first(z)), z) || throw(ArgumentError(
        "spatial plots require a mesh in a plane of constant axial coordinate"))
    return nothing
end

_fem_mesh_input(source::Engine.FEMMesh) = source
_fem_mesh_input(source::AbstractString) = ImportExport.import_data(:msh, source)

include("fem_mesh.jl")

function _fem_field_samples(field; component, part, step)
    part in (:real, :imag, :magnitude, :phase) || throw(ArgumentError(
        "field part must be :real, :imag, :magnitude, or :phase"))
    step isa Integer && !(step isa Bool) && step > 0 || throw(ArgumentError(
        "field step must be a positive one-based index"))
    field.representation === :complex && step != 1 &&
        throw(ArgumentError(
            "a complex field contains one phasor, not two independent steps"))
    field.representation === :real && part in (:imag, :phase) &&
        throw(ArgumentError(
            "imaginary part and phase require representation=:complex when importing"))
    samples = Matrix{Float64}[]
    for block in field.blocks
        ncomponents, nodes, steps, elements = size(block.values)
        component === nothing ||
            (component isa Integer && !(component isa Bool) &&
             1 <= component <= ncomponents) ||
            throw(ArgumentError("field component is out of range"))
        ncomponents > 1 && component === nothing && part !== :magnitude &&
            throw(ArgumentError("select a field component for $part"))
        step <= steps || throw(ArgumentError("field step is out of range"))
        chosen = component === nothing ? (1:ncomponents) : (component:component)
        values = Matrix{Float64}(undef, nodes, elements)
        for element in 1:elements, node in 1:nodes

            if part === :magnitude
                values[node, element] = sqrt(sum(chosen) do c
                    abs2(block.values[c, node, step, element]) +
                    (field.representation === :complex ?
                     abs2(block.values[c, node, 2, element]) : 0.0)
                end)
            else
                c = first(chosen)
                re = block.values[c, node, step, element]
                im = field.representation === :complex ? block.values[c, node, 2, element] :
                     0.0
                values[node, element] = part === :real ? re :
                                        part === :imag ? im :
                                        iszero(re) && iszero(im) ? NaN : atan(im, re)
            end
        end
        push!(samples, values)
    end
    return samples
end

function _fem_layer!(
        axis, field::Engine.FEMFieldMap; component = nothing, part = :real, step = 1,
        colormap = :viridis, colorscale = identity, colorrange = nothing,
        arrows = false, arrow_stride = 1, arrow_attributes = (;))
    arrow_stride isa Integer && !(arrow_stride isa Bool) && arrow_stride > 0 ||
        throw(ArgumentError("arrow_stride must be a positive integer"))
    arrows && part ∉ (:real, :imag) &&
        throw(ArgumentError(
            "vector arrows require a real or imaginary component field"))
    samples = _fem_field_samples(field; component, part, step)
    points = GeometryBasics.Point2d[]
    faces = GeometryBasics.TriangleFace{Int}[]
    values = Float64[]
    origins = GeometryBasics.Point2d[]
    directions = GeometryBasics.Vec2d[]
    element_index = 0
    for (block, data) in zip(field.blocks, samples)
        block.element_type in (2, 3) || throw(ArgumentError(
            "field plotting supports triangle and quadrangle samples; found element type $(block.element_type)"))
        _fem_planar(reshape(block.coordinates, 3, :))
        nodes = size(block.coordinates, 2)
        for element in axes(block.coordinates, 3)
            offset = length(points)
            for node in 1:nodes
                push!(points,
                    GeometryBasics.Point2d(block.coordinates[1, node, element],
                        block.coordinates[2, node, element]))
                push!(values, data[node, element])
            end
            push!(faces, GeometryBasics.TriangleFace{Int}(offset+1, offset+2, offset+3))
            nodes == 4 &&
                push!(faces, GeometryBasics.TriangleFace{Int}(offset+1, offset+3, offset+4))
            element_index += 1
            if arrows && mod(element_index-1, arrow_stride) == 0
                size(block.values, 1) == 3 ||
                    throw(ArgumentError("arrows require a vector field"))
                selected = part === :imag ? 2 : step
                push!(origins,
                    GeometryBasics.Point2d(
                        sum(block.coordinates[1, :, element])/nodes,
                        sum(block.coordinates[2, :, element])/nodes))
                push!(directions,
                    GeometryBasics.Vec2d(
                        sum(block.values[1, :, selected, element])/nodes,
                        sum(block.values[2, :, selected, element])/nodes))
            end
        end
    end
    finite = filter(isfinite, values)
    isempty(finite) &&
        throw(ArgumentError("field contains no defined finite color samples"))
    logarithmic = colorscale in (log, log2, log10)
    excluded = count(value -> !isfinite(value) || (logarithmic && value <= 0), values)
    if logarithmic
        filter!(>(0), finite)
        isempty(finite) &&
            throw(ArgumentError("logarithmic field colors require positive samples"))
    end
    values = [isfinite(value) && (!logarithmic || value > 0) ? value : NaN
              for value in values]
    limits = colorrange === nothing ? extrema(finite) : colorrange
    length(limits) == 2 && all(isfinite, limits) && first(limits) <= last(limits) ||
        throw(ArgumentError("colorrange must contain two finite ordered values"))
    logarithmic && first(limits) <= 0 &&
        throw(ArgumentError("logarithmic colorrange must be positive"))
    if first(limits) == last(limits)
        value = first(limits)
        padding = max(abs(value)*0.01, eps(Float64))
        limits = logarithmic ? (value/1.01, value*1.01) : (value-padding, value+padding)
    end
    plot = Makie.mesh!(axis, points, faces; color = values, colormap, colorscale,
        colorrange = limits, shading = Makie.NoShading)
    layers = Any[plot]
    arrows && push!(layers, Makie.arrows2d!(axis, origins, directions;
        merge((; color = :black), arrow_attributes)...))
    selection = component === nothing ? "" : "; component=$component"
    unit_label = part === :phase ? "Phase [rad]" : string(part)
    label = "$unit_label$selection\n$(replace(field.label, "; " => "\n"))"
    field.representation === :real && (label *= "; step=$step")
    scale = (; colormap, limits, ticks = Makie.automatic, label)
    return (; layers, scales = (scale,), excluded,
        colorscale, group = :field, label = field.label)
end

function _addon_fem_plot(source, attributes;
        geometry = nothing, mesh = nothing, title = nothing, figure_title = nothing,
        title_attributes = (;), size = (1000, 700), backend = nothing,
        display_plot = true, controls = true, export_theme = :default, open_export = true,
        legend_position = nothing, legend_attributes = (;), colorbar_position = :right,
        colorbar_attributes = (;), mesh_color = (:black, 0.25), mesh_linewidth = 0.4, kwargs...)
    _addon_activate_backend(backend)
    overlay = mesh === nothing ? nothing : _fem_mesh_input(mesh)
    display_title = title === nothing ? basename(source.source) : String(title)
    return with_theme(_addon_theme(; export_theme)) do
        shell = _addon_shell(; size, controls, kwargs...)
        panel = _addon_panel!(shell, (1, 1))
        axis = Axis(panel.content;
            merge(
                (; xlabel = "y [m]", ylabel = "z [m]",
                    title = display_title, aspect = DataAspect(), tellwidth = false, tellheight = false),
                shell.axis_attributes)...)
        rendered = _fem_layer!(axis, source; attributes...)
        sidebar = source isa Engine.FEMMesh ?
            _fem_mesh_sidebar!(shell, panel, rendered.mesh_view; controls) : nothing
        groups = Dict{Symbol, Vector{Any}}(rendered.group=>rendered.layers)
        labels = Dict{Symbol, Any}(rendered.group=>rendered.label)
        order = [rendered.group]
        scales = rendered.scales
        scale_attributes = merge((; scale = rendered.colorscale), colorbar_attributes)
        if overlay !== nothing
            groups[:mesh_overlay] = _fem_layer!(axis, overlay; color = mesh_color, linewidth = mesh_linewidth).layers
            push!(order, :mesh_overlay)
            labels[:mesh_overlay] = "Mesh overlay"
        end
        if geometry !== nothing
            geometry isa DataModel.LineCableSystem ||
                throw(ArgumentError("geometry must be a LineCableSystem"))
            polygons, references = _native_system_shapes(geometry, false; display_dielectric_pattern = false)
            groups[:geometry] = Any[poly!(axis, polygon.geometry; color = :transparent,
                                        strokecolor = :black, strokewidth = 0.8)
                                    for polygon in polygons]
            append!(groups[:geometry],
                [hlines!(axis, reference.values; color = reference.color,
                     linewidth = reference.width, xautolimits = false, yautolimits = false)
                 for reference in references])
            push!(order, :geometry)
            labels[:geometry] = "Geometry"
        end
        p = _addon_finish!(shell, Any[axis], Function[_addon_reset!(axis)],
            groups, order, labels; title = display_title, figure_title, title_attributes,
            legend_position, legend_attributes, panels = (panel,), color_scales = scales,
            colorbar_position = isempty(scales) ? nothing : colorbar_position,
            colorbar_attributes = scale_attributes, controls, display_plot=false,
            export_name = splitext(basename(source.source))[1], export_theme, open_export)
        p.addon_state = merge(p.addon_state, (; spatial_data = source))
        sidebar === nothing || _fem_mesh_controls!(p, rendered.mesh_view, sidebar)
        rendered.excluded > 0 &&
            (p.status[] = "$(rendered.excluded) undefined or nonpositive color samples are not colored")
        display_plot && _addon_display!(p.figure, display_title)
        p
    end
end

function plot(mesh::Engine.FEMMesh; color = (:black, 0.45), linewidth = 0.5,
        color_by=:uniform, inspect=:none, kwargs...)
    _addon_fem_plot(mesh, (; color, linewidth, color_by, inspect); kwargs...)
end

function plot(field::Engine.FEMFieldMap; component = nothing, part = :real, step = 1,
        colormap = :viridis, colorscale = identity, colorrange = nothing,
        arrows = false, arrow_stride = 1, arrow_attributes = (;), kwargs...)
    _addon_fem_plot(field,
        (; component, part, step, colormap, colorscale,
            colorrange, arrows, arrow_stride, arrow_attributes);
        kwargs...)
end

function Makie.plot(source::Union{Engine.FEMMesh, Engine.FEMFieldMap}; kwargs...)
    plot(source; kwargs...)
end

function plot(path::AbstractString; view = nothing,
        representation = :auto, coordinate_scale = 1, kwargs...)
    extension = lowercase(splitext(path)[2])
    if extension == ".msh"
        view === nothing && representation === :auto || throw(ArgumentError(
            "view and representation apply to field files"))
        return plot(ImportExport.import_data(:msh, path; coordinate_scale); kwargs...)
    elseif extension == ".pos"
        field = ImportExport.import_data(:pos, path; view, representation, coordinate_scale)
        field isa AbstractVector && throw(ArgumentError(
            "file contains multiple field views; select one with view=..."))
        return plot(field; kwargs...)
    end
    throw(ArgumentError("spatial plot input must be a .msh or .pos file"))
end
