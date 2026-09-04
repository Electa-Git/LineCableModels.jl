@enum FEMRunState begin
    created
    geometry_ready
    mesh_ready
    running
    completed
    failed
    cancelled
    not_executed
end

mutable struct FEMRun
    path::String
    state::FEMRunState
    message::String
    mesh_source::Symbol
    mesh_fingerprint::String
    getdp_invocations::Int
end

function FEMRun(path, state, message, mesh_source, mesh_fingerprint)
    FEMRun(
        path, state, message, mesh_source, mesh_fingerprint, 0
    )
end

struct FEMMaterialPlan{T <: Real}
    object_id::String
    field::Symbol
    kind::Symbol
    rho::T
    eps_r::T
    mu_r::T
    tan_delta::T
    physical_tag::Int
    physical_name::String
end

struct FEMRegionPlan{T <: Real, S}
    object_id::String
    cable_index::Int
    region_index::Int
    terminal_index::Int
    material_index::Int
    shape::S
    mesh_size::T
end

struct FEMMeshPlan{T <: Real}
    frequency_index::Int
    frequency::T
    domain_radius::T
    shell_outer_radius::T
    domain_mesh_size::T
    infinite_mesh_size::T
    interface_mesh_size::T
    cable_interface_mesh_sizes::Vector{T}
end

struct FEMResolvedModel{T <: Real, P <: LineParametersProblem}
    problem::P
    terminal_ids::Vector{String}
    terminal_names::Vector{String}
    region_plans::Vector{FEMRegionPlan}
    material_plans::Vector{FEMMaterialPlan{T}}
    cable_boundaries::Vector{Any}
    cable_hosts::Vector{Symbol}
    tags::NamedTuple
    centre::Tuple{T, T}
    domain_radius::T
    shell_outer_radius::T
    fine_mesh_size::T
    coarse_mesh_size::T
    maximum_frequency::T
    cable_outer_mesh_sizes::Vector{T}
    mesh_growth_factor::T
    mesh_plans::Vector{FEMMeshPlan{T}}
end

struct FEMGeometry
    model_name::String
    terminal_surfaces::Vector{Vector{Int}}
    terminal_curves::Vector{Vector{Int}}
    material_surfaces::Vector{Vector{Int}}
    cable_curves::Vector{Vector{Int}}
    air_surfaces::Vector{Int}
    earth_surfaces::Vector{Int}
    infinite_surfaces::Vector{Int}
    outer_curves::Vector{Int}
    outer_air_curves::Vector{Int}
    outer_earth_curves::Vector{Int}
    inner_shell_curves::Vector{Int}
    interface_curves::Vector{Int}
end

struct FEMScan{T <: Real}
    Z::Array{Complex{T}, 3}
    P::Array{Complex{T}, 3}
    map_paths::Vector{String}
end

struct FEMRunRecord
    state::FEMRunState
    run_directory::Union{Nothing, String}
    mesh_source::Symbol
    mesh_fingerprint::String
    getdp_invocations::Int
    map_paths::Vector{String}
end

function _fem_error(
        category::Symbol,
        object_id,
        field::Symbol,
        message::AbstractString;
        run_directory = nothing
)
    throw(LineCableModelsFEMError(
        category,
        object_id,
        field,
        message;
        run_directory
    ))
end

const FEM_FLOAT_TAGS = ("Float16", "Float32", "Float64", "BigFloat")

function _fem_float64_scalar(value::AbstractDict)
    scalar = ImportExport.deserialize_value(value)
    converted = Float64(scalar)
    isfinite(scalar) && !isfinite(converted) &&
        throw(OverflowError(
            "finite value $scalar is outside the Float64 range"
        ))
    return ImportExport.serialize_value(converted)
end

function _fem_float64_document(value::AbstractVector)
    return [_fem_float64_document(item) for item in value]
end

function _fem_float64_document(value::AbstractDict)
    if haskey(value, "__type__")
        marker = String(value["__type__"])
        marker in FEM_FLOAT_TAGS && return _fem_float64_scalar(value)
        marker == "Measurement" && return _fem_float64_document(
            get(value, "value") do
            throw(ArgumentError("serialized Measurement has no nominal value"))
        end
        )
        marker == "Complex" && return Dict(
            "__type__" => "Complex",
            "re" => _fem_float64_document(get(value, "re") do
                throw(ArgumentError("serialized Complex has no real component"))
            end),
            "im" => _fem_float64_document(get(value, "im") do
                throw(ArgumentError("serialized Complex has no imaginary component"))
            end)
        )
    end
    return Dict(
        String(key) => _fem_float64_document(item) for (key, item) in value
    )
end

_fem_float64_document(value) = value

"""
Rebuild one FEM problem from uncertainty-free `Float64` declarations.

The conversion is deliberately performed before adaptation or runtime setup.
Integer-valued topology is retained as integer data; all serialized floating
values, including the nominal component of `Measurements.Measurement`, become
`Float64`.
"""
function _preflight_fem_problem(problem::LineParametersProblem)
    object_id = problem.system.system_id
    document = try
        _fem_float64_document(ImportExport.serialize_value(problem))
    catch exception
        _fem_error(
            :preflight,
            object_id,
            :numeric_type,
            "could not normalize FEM inputs to nominal Float64: " *
            sprint(showerror, exception)
        )
    end
    normalized = try
        system = ImportExport.deserialize_value(document["system"])
        temperature = ImportExport.deserialize_value(document["temperature"])
        earth_props = ImportExport.deserialize_value(document["earth_props"])
        frequencies = ImportExport.deserialize_value(document["frequencies"])
        propagation = ImportExport.deserialize_value(document["Gamma"])
        LineParametersProblem(
            system;
            temperature,
            earth_props,
            frequencies,
            Γ = propagation
        )
    catch exception
        _fem_error(
            :preflight,
            object_id,
            :numeric_type,
            "could not rebuild the nominal Float64 FEM problem: " *
            sprint(showerror, exception)
        )
    end
    eltype(normalized) === Float64 || _fem_error(
        :preflight,
        object_id,
        :numeric_type,
        "FEM preflight produced $(eltype(normalized)); expected Float64"
    )
    return normalized
end

function _validate_material(material, object_id::String)
    material.kind === :conductor && !isfinite(material.rho) &&
        _fem_error(
            :adaptation,
            object_id,
            :rho,
            "a conductor requires finite electrical resistivity"
        )
    isfinite(material.eps_r) || _fem_error(
        :adaptation, object_id, :eps_r, "relative permittivity must be finite"
    )
    isfinite(material.mu_r) || _fem_error(
        :adaptation, object_id, :mu_r, "relative permeability must be finite"
    )
    isfinite(material.tan_delta) || _fem_error(
        :adaptation, object_id, :tan_delta, "loss tangent must be finite"
    )
    return nothing
end

function _validate_fem_shape(
        shape::Union{
            DataModel.Disk,
            DataModel.Rectangle,
            DataModel.Ellipse,
            DataModel.Annulus,
            DataModel.Polygon,
            DataModel.BentStrip,
            DataModel.SectorShape},
        object_id::String
)
    return nothing
end

function _validate_fem_shape(shape::DataModel.ShellShape, object_id::String)
    _validate_fem_shape(shape.inner, object_id)
    _validate_fem_shape(shape.outer, object_id)
    return nothing
end

function _validate_fem_shape(shape::DataModel.DifferenceShape, object_id::String)
    _validate_fem_shape(shape.outer, object_id)
    foreach(hole -> _validate_fem_shape(hole, object_id), shape.holes)
    return nothing
end

function _validate_fem_shape(shape::DataModel.AssemblyShape, object_id::String)
    foreach(member -> _validate_fem_shape(member, object_id), shape.members)
    return nothing
end

function _validate_fem_shape(shape, object_id::String)
    _fem_error(
        :unsupported,
        object_id,
        :primitive,
        "resolved shape $(typeof(shape)) has no built-in geo-kernel adaptation"
    )
end

function _temperature_resistivity(material, problem, formulation)
    formulation.options.temperature_correction || return material.rho
    isfinite(material.rho) || return material.rho
    return material.rho * (
        one(material.rho) + material.alpha * (problem.temperature - material.T0)
    )
end

function _validate_material_partition(design)
    envelope_area = LineCableModels.area(design.geometry.outer)
    declared_area = sum(
        region -> LineCableModels.area(region.primitive),
        design.geometry.regions;
        init = zero(envelope_area)
    )
    tolerance = max(abs(envelope_area), one(envelope_area)) * 1e-10
    isapprox(
        declared_area,
        envelope_area;
        rtol = 1e-10,
        atol = tolerance
    ) || _fem_error(
        :adaptation,
        design.cable_id,
        :material_partition,
        "the resolved cable cross-section is not a complete material " *
        "partition: declared area $(declared_area), envelope area " *
        "$(envelope_area); represent interstitial media explicitly with " *
        "Enclosure"
    )
    return nothing
end

function _fem_region_mesh_size(region::DataModel.PlacedRegion)
    source = region.source.primitive
    shape = region.primitive
    repeated = !isempty(region.placement.patterns)
    if source isa DataModel.Disk
        return repeated ? source.r : source.r / 5
    end
    if source isa DataModel.Annulus
        return (source.ro - source.ri) / 2
    end
    if source isa DataModel.Rectangle
        return min(source.w, source.h) / (repeated ? 2 : 5)
    end
    if source isa DataModel.Ellipse
        return min(source.a, source.b) / (repeated ? 1 : 5)
    end
    if shape isa DataModel.Annulus
        return (shape.ro - shape.ri) / 2
    end
    region_area = LineCableModels.area(shape)
    region_perimeter = DataModel.perimeter(shape)
    scale = region_area / region_perimeter
    isfinite(scale) && scale > zero(scale) || _fem_error(
        :adaptation,
        String(region.source.tag),
        :mesh_size,
        "could not derive a positive local characteristic length"
    )
    return scale
end

function _fem_cable_outer_mesh_sizes(
        region_plans::Vector{FEMRegionPlan}, cable_boundaries, cable_count::Int
)
    sizes = Vector{Float64}(undef, cable_count)
    for cable_index in 1:cable_count
        indices = findall(plan -> plan.cable_index == cable_index, region_plans)
        isempty(indices) && _fem_error(
            :adaptation,
            "cable_$cable_index",
            :mesh_size,
            "a cable has no resolved material regions"
        )
        boundary = cable_boundaries[cable_index]
        sizes[cable_index] = if boundary isa DataModel.AssemblyShape
            minimum(region_plans[index].mesh_size for index in indices)
        else
            region_plans[last(indices)].mesh_size
        end
    end
    return sizes
end

function _fem_mesh_plans(
        problem::LineParametersProblem{T},
        centre_x::T,
        layout_radius::T,
        cable_outer_mesh_sizes::Vector{T},
        growth_factor::T
) where {T <: Real}
    earth = problem.earth_props.layers[2]
    plans = FEMMeshPlan{T}[]
    for (frequency_index, frequency) in enumerate(problem.frequencies)
        earth_skin_depth = sqrt(
            earth.rho /
            (convert(T, π) * frequency * earth.mu_r * convert(T, 4π * 1e-7))
        )
        domain_radius = max(layout_radius, convert(T, earth_skin_depth))
        shell_outer_radius = convert(T, 1.25) * domain_radius
        domain_mesh_size = domain_radius / 20
        infinite_mesh_size = 2domain_mesh_size
        cable_interface_mesh_sizes = T[]
        for (cable_index, (design, position)) in enumerate(zip(
            problem.system.designs, problem.system.positions
        ))
            distance = max(
                zero(T), abs(position.y) - LineCableModels.outer_radius(design)
            )
            push!(cable_interface_mesh_sizes,
                min(
                    domain_mesh_size,
                    cable_outer_mesh_sizes[cable_index] +
                    (growth_factor - one(T)) * distance
                ))
        end
        interface_mesh_size = minimum([domain_mesh_size; cable_interface_mesh_sizes])
        push!(plans,
            FEMMeshPlan(
                frequency_index,
                frequency,
                domain_radius,
                shell_outer_radius,
                domain_mesh_size,
                infinite_mesh_size,
                interface_mesh_size,
                cable_interface_mesh_sizes
            ))
    end
    return plans
end

function _formations(regions, terminal_map, object_id)
    formations = NamedTuple[]
    assigned = falses(length(regions))
    for (region_index, region) in pairs(regions)
        assigned[region_index] && continue
        entry_index = findfirst(
            entry -> entry.pattern isa DataModel.BoundedPlacement,
            region.placement.patterns
        )
        entry_index === nothing && continue
        entry = region.placement.patterns[entry_index]
        enclosing = region.placement.patterns[(entry_index + 1):end]
        members = Int[]
        member_ids = Int[]
        for (peer_index, peer) in pairs(regions)
            peer_entry_index = findfirst(
                candidate -> candidate.pattern isa DataModel.BoundedPlacement,
                peer.placement.patterns
            )
            peer_entry_index === nothing && continue
            peer_entry = peer.placement.patterns[peer_entry_index]
            isequal(peer_entry.pattern.boundary, entry.pattern.boundary) || continue
            peer.terminal === region.terminal || continue
            isequal(
                peer.placement.patterns[(peer_entry_index + 1):end], enclosing
            ) || continue
            push!(members, peer_index)
            push!(member_ids, peer_entry.member)
        end
        order = sortperm(member_ids)
        members = members[order]
        member_ids = member_ids[order]
        member_ids == collect(1:length(member_ids)) || _fem_error(
            :adaptation,
            object_id,
            :material_partition,
            "bounded-formation member identities must be contiguous from one"
        )
        assigned[members] .= true

        boundary = entry.pattern.boundary
        occupied_area = sum(
            index -> LineCableModels.area(regions[index].source.primitive),
            members
        )
        boundary_area = LineCableModels.area(boundary)
        complete = isapprox(
            occupied_area,
            boundary_area;
            rtol = 5.0e-6,
            atol = zero(boundary_area)
        )
        member_shapes = [regions[index].primitive for index in members]
        if complete
            material = regions[first(members)].source.material
            all(
                index -> regions[index].source.material == material,
                members
            ) || _fem_error(
                :unsupported,
                object_id,
                :material_partition,
                "a complete bounded formation can collapse only when every " *
                "member has the same material"
            )
            terminal = terminal_map[first(members)]
            all(index -> terminal_map[index] == terminal, members) || _fem_error(
                :unsupported,
                object_id,
                :terminal,
                "a complete bounded formation can collapse only when every " *
                "member has the same terminal owner"
            )
        else
            filled = any(regions) do candidate
                shape = candidate.primitive
                candidate.source.material.kind === :conductor && return false
                shape isa DataModel.DifferenceShape || return false
                all(member_shapes) do member_shape
                    any(hole -> isequal(hole, member_shape), shape.holes)
                end
            end
            filled || _fem_error(
                :adaptation,
                object_id,
                :material_partition,
                "an incomplete bounded formation leaves unassigned cross-sectional " *
                "area; contain it in Enclosure with an explicit fill material"
            )
        end
        push!(formations, (;
            members,
            member_shapes,
            boundary,
            complete,
            mesh_size = minimum(
                index -> _fem_region_mesh_size(regions[index]), members
            )
        ))
    end
    return formations
end

function _coalesce(shape::DataModel.DifferenceShape, formations, object_id)
    holes = Any[shape.holes...]
    for formation in formations
        formation.complete || continue
        matched = findall(eachindex(holes)) do hole_index
            any(
                member_shape -> isequal(holes[hole_index], member_shape),
                formation.member_shapes
            )
        end
        isempty(matched) && continue
        length(matched) == length(formation.member_shapes) || _fem_error(
            :adaptation,
            object_id,
            :material_partition,
            "an enclosing material excludes only part of a collapsed bounded formation"
        )
        insertion = first(matched)
        retained = Any[]
        for hole_index in eachindex(holes)
            hole_index == insertion && push!(retained, formation.boundary)
            hole_index in matched || push!(retained, holes[hole_index])
        end
        holes = retained
    end
    return DataModel.DifferenceShape(shape.outer, Tuple(holes))
end

_coalesce(shape, ::Any, ::Any) = shape

function _resolved_fem_model(
        problem::LineParametersProblem{T},
        formulation::LineCableModelsFEM
) where {T <: Real}
    try
        LineCableModels.validate(problem)
    catch exception
        _fem_error(
            :adaptation,
            problem.system.system_id,
            :problem,
            sprint(showerror, exception)
        )
    end
    problem.Γ === nothing || _fem_error(
        :unsupported,
        problem.system.system_id,
        :Γ,
        "LineCableModelsFEM uses its fixed quasi-TEM propagation constant; " *
        "problem-level propagation constants are unsupported"
    )
    earth = problem.earth_props
    earth.vertical_layers && _fem_error(
        :unsupported,
        problem.system.system_id,
        :vertical_layers,
        "vertical earth layers are not supported by the two-dimensional FEM domain"
    )
    length(earth.layers) == 2 || _fem_error(
        :unsupported,
        problem.system.system_id,
        :earth_props,
        "the FEM backend currently supports one homogeneous earth half-space"
    )
    environment = problem.system.environment
    environment isa Union{Nothing, EarthProps.EarthModel} || _fem_error(
        :unsupported,
        problem.system.system_id,
        :environment,
        "the declared environment type $(typeof(environment)) has no FEM adaptation"
    )

    system = problem.system
    terminal_count = length(system.terminal_order)
    terminal_count > 0 || _fem_error(
        :adaptation, system.system_id, :terminal_order, "no terminals were resolved"
    )
    terminal_ids = [@sprintf("cable_%04d/%s/%s",
                        entry.cable,
                        system.designs[entry.cable].cable_id,
                        entry.terminal)
                    for entry in system.terminal_order]
    length(unique(terminal_ids)) == terminal_count || _fem_error(
        :adaptation,
        system.system_id,
        :terminal_order,
        "global cable/terminal identities must be unique"
    )
    terminal_names = [@sprintf("LCM/terminal/%04d/%s", index, terminal_ids[index])
                      for index in eachindex(terminal_ids)]

    region_plans = FEMRegionPlan[]
    material_plans = FEMMaterialPlan{T}[]
    global_region = 0
    for (cable_index, design) in enumerate(system.designs)
        _validate_material_partition(design)
        count = length(design.geometry.regions)
        first_global = global_region + 1
        last_global = global_region + count
        regions = @view system.geometry[first_global:last_global]
        terminals = @view system.terminal_map[first_global:last_global]
        cable_id = @sprintf("cable_%04d/%s", cable_index, design.cable_id)
        formations = _formations(regions, terminals, cable_id)
        collapsed = Dict(
            first(formation.members) => formation
            for formation in formations if formation.complete
        )
        skipped = Set(
            member
            for formation in formations if formation.complete
            for member in formation.members[2:end]
        )
        for local_region in eachindex(design.geometry.regions)
            global_region += 1
            local_region in skipped && continue
            placed = system.geometry[global_region]
            source = placed.source
            object_id = @sprintf("cable_%04d/%s/%s/region_%04d",
                cable_index,
                design.cable_id,
                source.tag,
                local_region)
            _validate_material(source.material, object_id)
            terminal_index = system.terminal_map[global_region]
            if source.material.kind === :conductor
                terminal_index > 0 || _fem_error(
                    :adaptation,
                    object_id,
                    :terminal,
                    "a conductor surface has no electrical terminal owner"
                )
            elseif terminal_index != 0
                _fem_error(
                    :adaptation,
                    object_id,
                    :terminal,
                    "a passive material surface cannot own an electrical terminal"
                )
            end
            rho = convert(T, _temperature_resistivity(
                source.material, problem, formulation
            ))
            formation = get(collapsed, local_region, nothing)
            shape = formation === nothing ?
                    _coalesce(placed.primitive, formations, object_id) :
                    formation.boundary
            mesh_size = formation === nothing ?
                        _fem_region_mesh_size(placed) : formation.mesh_size
            _validate_fem_shape(shape, object_id)
            material_index = length(material_plans) + 1
            physical_tag = 10_000 + material_index
            physical_name = @sprintf("LCM/material/%04d/%s/%s",
                material_index,
                source.material.kind,
                object_id)
            push!(material_plans,
                FEMMaterialPlan{T}(
                    object_id,
                    source.tag,
                    source.material.kind,
                    rho,
                    convert(T, source.material.eps_r),
                    convert(T, source.material.mu_r),
                    convert(T, source.material.tan_delta),
                    physical_tag,
                    physical_name
                ))
            push!(region_plans,
                FEMRegionPlan(
                    object_id,
                    cable_index,
                    local_region,
                    terminal_index,
                    material_index,
                    shape,
                    convert(T, mesh_size)
                ))
        end
    end
    global_region == length(system.geometry) || _fem_error(
        :adaptation,
        system.system_id,
        :geometry,
        "cable-design and system geometry orders are inconsistent"
    )

    cable_boundaries = Any[LineCableModels.resolve(position, design.geometry.outer)
                           for (design, position) in zip(system.designs, system.positions)]
    for (index, boundary) in enumerate(cable_boundaries)
        _validate_fem_shape(boundary, @sprintf("cable_%04d/%s", index,
            system.designs[index].cable_id))
    end
    cable_hosts = Symbol[]
    for (design, position) in zip(system.designs, system.positions)
        radius = LineCableModels.outer_radius(design)
        if position.y - radius > 0
            push!(cable_hosts, :air)
        elseif position.y + radius < 0
            push!(cable_hosts, :earth)
        else
            _fem_error(
                :unsupported,
                design.cable_id,
                :position,
                "a cable cannot cross or touch the air/earth interface in " *
                "the two-dimensional FEM adaptation"
            )
        end
    end
    centre_x = convert(T, sum(position.x for position in system.positions) /
                          length(system.positions))
    maximum_cable_distance = maximum(
        (
            hypot(
                system.positions[right].x - system.positions[left].x,
                system.positions[right].y - system.positions[left].y
            )
        for left in eachindex(system.positions)
        for right in (left + 1):length(system.positions)
        );
        init = zero(T))
    envelope_radius = maximum(
        hypot(position.x - centre_x, position.y) +
        LineCableModels.outer_radius(design)
    for (design, position) in zip(system.designs, system.positions)
    )
    layout_radius = max(
        convert(T, 5),
        convert(T, 2) * maximum_cable_distance,
        convert(T, envelope_radius)
    )
    cable_outer_mesh_sizes = convert(
        Vector{T},
        _fem_cable_outer_mesh_sizes(
            region_plans, cable_boundaries, length(system.designs)
        )
    )
    mesh_growth_factor = convert(T, 1.2)
    mesh_plans = _fem_mesh_plans(
        problem,
        centre_x,
        layout_radius,
        cable_outer_mesh_sizes,
        mesh_growth_factor
    )
    display_plan = last(mesh_plans)
    domain_radius = display_plan.domain_radius
    shell_outer_radius = display_plan.shell_outer_radius
    fine_mesh_size = minimum(getproperty.(region_plans, :mesh_size))
    coarse_mesh_size = display_plan.domain_mesh_size
    maximum_frequency = maximum(problem.frequencies)
    tags = (
        air = 1_001,
        earth = 1_002,
        outer_boundary = 2_001,
        interface = 2_002,
        inner_shell_boundary = 2_003,
        outer_air_boundary = 2_004,
        outer_earth_boundary = 2_005,
        air_infinite = 1_003,
        earth_infinite = 1_004,
        infinite_domain = 1_005,
        terminal_base = 3_000,
        terminal_contour_base = 4_000,
        cable_contour_base = 5_000
    )
    return FEMResolvedModel(
        problem,
        terminal_ids,
        terminal_names,
        region_plans,
        material_plans,
        cable_boundaries,
        cable_hosts,
        tags,
        (centre_x, zero(T)),
        domain_radius,
        shell_outer_radius,
        fine_mesh_size,
        coarse_mesh_size,
        maximum_frequency,
        cable_outer_mesh_sizes,
        mesh_growth_factor,
        mesh_plans
    )
end
