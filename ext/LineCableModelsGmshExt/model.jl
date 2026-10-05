@enum FEMRunState begin
    created
    geometry_ready
    mesh_ready
    running
    completed
    failed
    cancelled
end

mutable struct FEMRun
    path::String
    state::FEMRunState
    message::String
    mesh_source::Symbol
    mesh_fingerprint::String
    getdp_invocations::Int
    completed_columns::Int
    completed_frequencies::Int
end

function FEMRun(path, state, message, mesh_source, mesh_fingerprint, invocations = 0)
    FEMRun(
        path, state, message, mesh_source, mesh_fingerprint, invocations, 0, 0
    )
end

struct FEMMaterialPlan{T <: Real}
    object_id::String
    field::Symbol
    kind::Symbol
    mu_r::T
    admittivity::Vector{Complex{T}}
    physical_tag::Int
    physical_name::String
end

struct FEMRegionPlan{S}
    object_id::String
    cable_index::Int
    region_index::Int
    terminal_index::Int
    material_index::Int
    shape::S
end

struct FEMResolvedModel{T <: Real, P <: LineParametersProblem}
    problem::P
    terminal_ids::Vector{String}
    terminal_names::Vector{String}
    region_plans::Vector{FEMRegionPlan}
    material_plans::Vector{FEMMaterialPlan{T}}
    earth_materials::Vector{Earth.EarthMaterial{T}}
    cable_boundaries::Vector{Any}
    cable_hosts::Vector{Symbol}
    tags::NamedTuple
    centre::Tuple{T, T}
    cad_scale::T
    prescribed_gamma::Vector{Complex{T}}
end

# Native topology of the graded strips and small core inside one sector.
struct FEMSectorPartition
    core::Int
    curve_pairs::Vector{Tuple{Int, Int, Float64, Float64}} # outer, inner, length, turn
    spokes::Vector{Int} # oriented from the physical contour towards the core
    depth::Float64
    patches::Dict{Int, NTuple{4, Int}}
end

# Intrinsic conductor section geometry; numerical grading belongs to mesh.geo.
_conductor_section(shape) = nothing
_conductor_section(shape::DataModel.Disk) = (; kind=1, width=shape.r)
_conductor_section(shape::DataModel.Annulus) = (; kind=2, width=shape.ro-shape.ri)
_conductor_section(shape::DataModel.SectorShape) = (; kind=3, width=shape.primitive.r_back)

struct FEMScan{T <: Real}
    Z::Array{Complex{T}, 3}
    P::Array{Complex{T}, 3}
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

function _fem_float64_scalar(value::Real)
    scalar = LineCableModels.nominal(value)
    converted = Float64(scalar)
    isfinite(scalar) && !isfinite(converted) && throw(OverflowError(
        "finite value $scalar is outside the Float64 range"))
    return converted
end

function _fem_float64_scalar(value::Complex)
    return complex(_fem_float64_scalar(real(value)), _fem_float64_scalar(imag(value)))
end

function _fem_float64_scalar(value::AbstractDict)
    return ImportExport.serialize_value(
        _fem_float64_scalar(ImportExport.deserialize_value(value)))
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
        LineParametersProblem(
            system;
            temperature,
            earth_props,
            frequencies
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
    if material isa LineCableModels.RadialDielectric
        foreach(m -> _validate_material(m, object_id), material.materials)
    else
        isfinite(material.tan_delta) || _fem_error(
            :adaptation, object_id, :tan_delta, "loss tangent must be finite"
        )
    end
    return nothing
end

function _validate_fem_shape(
        shape::Union{
            DataModel.Disk,
            DataModel.Rectangle,
            DataModel.Ellipse,
            DataModel.EllipseOffset,
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
            index -> LineCableModels.area(regions[index].primitive),
            members
        )
        boundary_area = LineCableModels.area(boundary)
        complete = isapprox(
            occupied_area,
            boundary_area;
            rtol = 0,
            atol = DataModel.geometry_tolerance(boundary_area)
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
                if shape isa DataModel.Annulus && boundary isa DataModel.Disk
                    concentric = isapprox(shape.at.x, boundary.at.x) &&
                                 isapprox(shape.at.y, boundary.at.y)
                    concentric && isapprox(shape.ro, boundary.r) || return false
                    isapprox(occupied_area, π * shape.ri^2;
                        rtol = 5e-6, atol = 0) || return false
                    shift = DataModel.Pose2(-boundary.at.x, -boundary.at.y)
                    tolerance = 64eps(shape.ri)
                    return all(member_shapes) do member_shape
                        DataModel.support(DataModel.resolve(shift, member_shape)) <=
                        shape.ri + tolerance
                    end
                end
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
        push!(formations,
            (;
                members,
                member_shapes,
                boundary,
                complete
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
            "an enclosing material excludes only part of a complete bounded formation"
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
    supplied = formulation.options.data.Γ
    supplied isa AbstractVector && length(supplied) != length(problem.frequencies) &&
        throw(DimensionMismatch("FEM Γ must contain one value per frequency sample"))
    prescribed = supplied isa Number ? fill(Complex{T}(supplied), length(problem.frequencies)) :
        Complex{T}.(supplied)
    all(isfinite, prescribed) || throw(ArgumentError("FEM Γ must be finite in the solver precision"))
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
    environment isa Union{Nothing, Earth.EarthModel} || _fem_error(
        :unsupported,
        problem.system.system_id,
        :environment,
        "the declared environment type $(typeof(environment)) has no FEM adaptation"
    )

    soil = earth.layers[2]
    isinf(soil.thickness) || _fem_error(:unsupported, problem.system.system_id,
        :earth_props, "the FEM soil region must be a semi-infinite half-space")
    earth_materials = Earth.EarthMaterial{T}[]
    for frequency in problem.frequencies
        evaluated = try
            state = LineCableModels.constitutive(formulation.methods.earth_properties,
                Earth.EarthMaterial(soil), frequency)
            material = Earth.EarthMaterial(_fem_float64_scalar(state.rho),
                _fem_float64_scalar(state.eps_r), _fem_float64_scalar(state.mu_r))
            isfinite(inv(material.rho)) && inv(material.rho) >= 0 || throw(DomainError(
                material.rho, "soil conductivity must be nonnegative and finite"))
            material
        catch exception
            _fem_error(:adaptation, problem.system.system_id, :earth_properties,
                "soil constitutive evaluation at $frequency Hz: $(sprint(showerror, exception))")
        end
        push!(earth_materials, evaluated)
    end

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
        complete = Dict(
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
            material = source.material
            admittivity = try
                if material.kind === :conductor
                    rho = _fem_float64_scalar(LineCableModels.constitutive(
                        formulation.methods.temperature_dependence, material, problem.temperature))
                    isfinite(inv(rho)) && inv(rho) > 0 || throw(DomainError(rho,
                        "conductor conductivity must be positive and finite"))
                    epsilon = convert(T, material.eps_r * 8.8541878128e-12)
                    Complex{T}[complex(inv(rho) + 2π * f * epsilon * material.tan_delta,
                        2π * f * epsilon) for f in problem.frequencies]
                else
                    selected = material isa LineCableModels.RadialDielectric ?
                        (formulation.methods.insulation_admittance,
                         formulation.methods.semicon_admittance) :
                        material.kind === :semicon ? formulation.methods.semicon_admittance :
                        formulation.methods.insulation_admittance
                    Complex{T}[_fem_float64_scalar(LineCableModels.constitutive(
                        selected, material, frequency, problem.temperature;
                        temperature_dependence=formulation.methods.temperature_dependence))
                        for frequency in problem.frequencies]
                end
            catch exception
                _fem_error(:adaptation, object_id, :constitutive,
                    "material evaluation at $(problem.temperature) °C: $(sprint(showerror, exception))")
            end
            all(isfinite, admittivity) || _fem_error(:adaptation, object_id,
                :constitutive, "evaluated material admittivity must be finite")
            formation = get(complete, local_region, nothing)
            shape = formation === nothing ?
                    _coalesce(placed.primitive, formations, object_id) :
                    formation.boundary
            _validate_fem_shape(shape, object_id)
            # Geometric regions and terminal ownership remain independent of
            # constitutive identity. Equal evaluated laws share one physical
            # material group, even across disconnected strand surfaces.
            mu_r = convert(T, material.mu_r)
            material_index = something(findfirst(material_plans) do plan
                plan.kind === material.kind && plan.mu_r == mu_r &&
                    plan.admittivity == admittivity
            end, 0)
            if iszero(material_index)
                material_index = length(material_plans) + 1
                push!(material_plans,
                    FEMMaterialPlan{T}(
                        object_id, source.tag, material.kind, mu_r, admittivity,
                        10_000 + material_index,
                        @sprintf("LCM/material/%04d/%s", material_index, material.kind)
                    ))
            end
            push!(region_plans,
                FEMRegionPlan(
                    object_id,
                    cable_index,
                    local_region,
                    terminal_index,
                    material_index,
                    shape
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
    # A frequency-independent scale for CAD vertex lookup only. Native point
    # sizes replace all provisional CAD sizes before meshing.
    cad_scale = convert(T,minimum(LineCableModels.outer_radius.(system.designs)))
    tags = (
        air = 1_001,
        earth = 1_002,
        outer_boundary = 2_001,
        interface = 2_002,
        pml_inner_boundary = 2_003,
        outer_air_boundary = 2_004,
        outer_earth_boundary = 2_005,
        air_pml = 1_003,
        earth_pml = 1_004,
        pml = 1_005,
        terminal_base = 3_000,
        terminal_contour_base = 4_000,
        measurement_line_base = 7_000,
        cable_contour_base = 5_000
    )
    return FEMResolvedModel(
        problem,
        terminal_ids,
        terminal_names,
        region_plans,
        material_plans,
        earth_materials,
        cable_boundaries,
        cable_hosts,
        tags,
        (centre_x, zero(T)),
        cad_scale,
        prescribed
    )
end
