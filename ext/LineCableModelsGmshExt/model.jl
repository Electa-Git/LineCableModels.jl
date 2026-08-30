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
    cable_boundaries::Vector{Any}
    cable_hosts::Vector{Symbol}
    tags::NamedTuple
    centre::Tuple{T, T}
    domain_radius::T
    shell_outer_radius::T
    fine_mesh_size::T
    coarse_mesh_size::T
    maximum_frequency::T
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

function _validate_fem_shape(shape, object_id::String)
    if shape isa Union{
        DataModel.Disk,
        DataModel.Rectangle,
        DataModel.Ellipse,
        DataModel.Sector,
        DataModel.Annulus,
        DataModel.Polygon,
        DataModel.RoundedSectorShape
    }
        return nothing
    end
    if hasproperty(shape, :inner) && hasproperty(shape, :outer)
        _validate_fem_shape(getproperty(shape, :inner), object_id)
        _validate_fem_shape(getproperty(shape, :outer), object_id)
        return nothing
    end
    if hasproperty(shape, :holes) && hasproperty(shape, :outer)
        _validate_fem_shape(getproperty(shape, :outer), object_id)
        for hole in getproperty(shape, :holes)
            _validate_fem_shape(hole, object_id)
        end
        return nothing
    end
    if hasproperty(shape, :members)
        for member in getproperty(shape, :members)
            _validate_fem_shape(member, object_id)
        end
        return nothing
    end
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
        for local_region in eachindex(design.geometry.regions)
            global_region += 1
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
            _validate_fem_shape(placed.primitive, object_id)
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
                    placed.primitive
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
    maximum_extent = maximum(
        LineCableModels.support(boundary)
    for boundary in cable_boundaries
    )
    x_min = minimum(
        position.x - LineCableModels.outer_radius(design)
    for (design, position) in zip(system.designs, system.positions)
    )
    x_max = maximum(
        position.x + LineCableModels.outer_radius(design)
    for (design, position) in zip(system.designs, system.positions)
    )
    centre_x = convert(T, (x_min + x_max) / 2)
    layout_radius = convert(T, max(
        4maximum_extent,
        4max(abs(x_min - centre_x), abs(x_max - centre_x)),
        one(T)
    ))

    maximum_frequency = maximum(problem.frequencies)
    minimum_frequency = minimum(problem.frequencies)
    earth_skin_depth = sqrt(
        earth.layers[2].rho /
        (convert(T, π) * minimum_frequency * earth.layers[2].mu_r *
         convert(T, 4π * 1e-7))
    )
    domain_radius = max(layout_radius, convert(T, earth_skin_depth))
    shell_outer_radius = convert(T, 1.25) * domain_radius
    conductor_scales = T[]
    for (region, material) in zip(region_plans, material_plans)
        material.kind === :conductor || continue
        geometric_scale = sqrt(convert(T, LineCableModels.area(region.shape))) / 4
        skin_depth = sqrt(
            material.rho /
            (convert(T, π) * maximum_frequency * material.mu_r * convert(T, 4π * 1e-7))
        )
        push!(conductor_scales, min(geometric_scale, skin_depth / 3))
    end
    isempty(conductor_scales) && _fem_error(
        :adaptation,
        system.system_id,
        :geometry,
        "at least one conductor region is required"
    )
    fine_mesh_size = max(minimum(conductor_scales), domain_radius * convert(T, 1e-6))
    coarse_mesh_size = max(domain_radius / 14, 8fine_mesh_size)
    tags = (
        air = 1_001,
        earth = 1_002,
        outer_boundary = 2_001,
        interface = 2_002,
        inner_shell_boundary = 2_003,
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
        maximum_frequency
    )
end
