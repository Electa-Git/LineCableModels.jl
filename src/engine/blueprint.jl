"""One physical circular material layer retained before conductor homogenization."""
struct BlueprintConductorLayer{T <: Real}
    r_in::T
    r_ex::T
    material::Material{T}
end

"""Physical input of the linear one-plus-six-wire ACSR reduction."""
struct BlueprintSingleLayerACSR{T <: Real}
    radius::T
    pitch::T
    core::Material{T}
    strands::Material{T}
end

"""Physical area, exposed perimeter, and material of one homogeneous section."""
struct BlueprintHomogeneousSection{T <: Real}
    area::T
    perimeter::T
    material::Material{T}
end

"""
$(TYPEDEF)

Store one frequency-independent equivalent coaxial conductor row.

$(TYPEDFIELDS)
"""
struct BlueprintConductor{T <: Real}
    "Retained terminal name."
    terminal::Symbol
    "Concentric assembly containing the terminal."
    assembly::Int
    "Equivalent inner radius [m]."
    r_in::T
    "Equivalent outer radius [m]."
    r_ex::T
    "Physical conductor cross-section [m²]."
    cross_section::T
    "Number of explicitly represented wires."
    num_wires::Int
    "Equivalent helical turns per unit length [1/m]."
    num_turns::T
    "Equivalent resistance at the material reference temperature [Ω/m]."
    resistance::T
    "Equivalent temperature coefficient [1/°C]."
    alpha::T
    "Equivalent geometric-mean radius [m]."
    gmr::T
    "Assembly-local centre in the design frame [m]."
    position::Tuple{T, T}
    "Artificial homogeneous conductor material."
    material::Material{T}
    "Physical concentric metal layers; empty for a nonradial or stranded reduction."
    layers::Vector{BlueprintConductorLayer{T}}
    "Physical one-plus-six-wire ACSR data, when the exact source geometry matches."
    acsr::Union{Nothing,BlueprintSingleLayerACSR{T}}
    "Homogeneous physical section, without strand or material homogenization."
    section::Union{Nothing,BlueprintHomogeneousSection{T}}
end

# Preserve the manual homogeneous-blueprint constructor.
function BlueprintConductor{T}(terminal,assembly,r_in,r_ex,cross_section,
        num_wires,num_turns,resistance,alpha,gmr,position,material) where {T <: Real}
    layers=BlueprintConductorLayer{T}[
        BlueprintConductorLayer{T}(r_in,r_ex,material)]
    return BlueprintConductor{T}(terminal,assembly,r_in,r_ex,cross_section,
        num_wires,num_turns,resistance,alpha,gmr,position,material,layers,nothing,
        BlueprintHomogeneousSection{T}(cross_section,2π*r_ex,material))
end

function _homogeneous_section(design,terminal,::Type{T}) where {T}
    sources=filter(region->region.terminal===terminal,design.geometry.regions)
    length(sources)==1 || return nothing
    source=only(sources)
    isempty(source.paths) && isempty(source.placement.patterns) || return nothing
    shape=source.primitive
    # A complete annulus carries isolated-conductor current at its outer
    # boundary. An open sector uses its entire connected contour.
    contour=shape isa DataModel.Annulus ? 2π*shape.ro : DataModel.perimeter(shape)
    return BlueprintHomogeneousSection{T}(T(DataModel.area(shape)),T(contour),
        convert(Material{T},source.source.material))
end

function _single_layer_acsr(design,terminal,conductor,::Type{T}) where {T}
    sources=filter(region->region.terminal===terminal,design.geometry.regions)
    length(sources)==7 && all(r->r.primitive isa DataModel.Disk,sources) ||
        return nothing
    points=[DataModel.centroid(r.primitive) for r in sources]
    offsets=[(p[1]-conductor.position[1],p[2]-conductor.position[2]) for p in points]
    distances=[hypot(p...) for p in offsets]
    centre=argmin(distances); R=T(sources[centre].primitive.r)
    tolerance=128eps(T)*max(R,one(T))
    distances[centre]<=tolerance && isempty(sources[centre].paths) || return nothing
    ring=setdiff(eachindex(sources),[centre])
    all(k->isapprox(T(sources[k].primitive.r),R) &&
        isapprox(distances[k],2R),ring) || return nothing
    angles=sort!([atan(offsets[k][2],offsets[k][1]) for k in ring])
    gaps=diff(vcat(angles,first(angles)+2*(one(T)*π)))
    all(x->isapprox(x,(one(T)*π)/3;rtol=128eps(T)),gaps) || return nothing
    metal=sources[first(ring)].source.material
    isapprox(metal.mu_r,one(T)) &&
        all(k->sources[k].source.material==metal,ring) || return nothing
    pitches=T[]
    directions=Int[]
    for k in ring
        paths=sources[k].paths
        if isempty(paths)
            push!(pitches,T(Inf));push!(directions,1)
        elseif length(paths)==1 && isapprox(only(paths).radius,2R)
            entry=only(paths)
            push!(pitches,T(DataModel.pitch(entry.path,entry.radius)))
            push!(directions,entry.path.dir)
        else
            return nothing
        end
    end
    all(==(first(pitches)),pitches) && all(==(first(directions)),directions) ||
        return nothing
    iszero(conductor.r_in) && isapprox(conductor.r_ex,3R) || return nothing
    return BlueprintSingleLayerACSR{T}(R,first(pitches),
        convert(Material{T},sources[centre].source.material),convert(Material{T},metal))
end

function _physical_conductor_layers(design,terminal,conductor,::Type{T}) where {T}
    layers=BlueprintConductorLayer{T}[]
    iszero(conductor.num_turns) || return layers
    for region in design.geometry.regions
        region.terminal===terminal || continue
        region.primitive isa Union{DataModel.Disk,DataModel.Annulus} ||
            return BlueprintConductorLayer{T}[]
        isempty(region.placement.patterns) && isempty(region.paths) ||
            return BlueprintConductorLayer{T}[]
        DataModel.same_radial_position(
            DataModel.conductor_zone_position([region]),conductor.position) ||
            return BlueprintConductorLayer{T}[]
        push!(layers,BlueprintConductorLayer{T}(
            T(DataModel.r_in(region.primitive)),T(DataModel.r_ex(region.primitive)),
            convert(Material{T},region.source.material)))
    end
    sort!(layers;by=layer->layer.r_in)
    return layers
end

"""
$(TYPEDEF)

Store one physical dielectric layer owned by a coaxial conductor interval.

$(TYPEDFIELDS)
"""
struct BlueprintDielectric{T <: Real}
    "Index of the conductor immediately inside this radial interval."
    conductor::Int
    "Layer inner radius [m]."
    r_in::T
    "Layer outer radius [m]."
    r_ex::T
    "Unmodified physical material."
    material::Material{T}
end

"""
$(TYPEDEF)

Store the frequency-independent, unreduced numerical description of one cable
design consumed by the coaxial backend.

The blueprint is the result of computational flattening. It retains equivalent
conductor annuli and every physical dielectric layer, but contains no evaluated
constitutive relation, frequency, temperature correction, earth property, or
matrix result.

$(TYPEDFIELDS)
"""
struct CableBlueprint{T <: Real}
    "Source cable identifier."
    cable_id::String
    "Equivalent conductors in DataModel terminal order."
    conductors::Vector{BlueprintConductor{T}}
    "Physical dielectric layers in radial order."
    dielectrics::Vector{BlueprintDielectric{T}}
    "Dielectric-layer range owned by every conductor row."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Contiguous conductor ranges for independent concentric assemblies."
    assembly_ranges::Vector{UnitRange{Int}}
    "Explicit nonconcentric common-pipe enclosure hierarchy."
    pipes::Vector{PipeAssembly{T}}

    function CableBlueprint{T}(
            cable_id::String,
            conductors::Vector{BlueprintConductor{T}},
            dielectrics::Vector{BlueprintDielectric{T}},
            dielectric_ranges::Vector{UnitRange{Int}},
            assembly_ranges::Vector{UnitRange{Int}},
            pipes::Vector{PipeAssembly{T}}=PipeAssembly{T}[]
    ) where {T <: Real}
        return validate(new{T}(
            cable_id,
            conductors,
            dielectrics,
            dielectric_ranges,
            assembly_ranges,
            pipes
        ))
    end
end

Base.eltype(::CableBlueprint{T}) where {T} = T
Base.eltype(::Type{<:CableBlueprint{T}}) where {T} = T
Base.length(blueprint::CableBlueprint) = length(blueprint.conductors)

function _assembly_ranges(components)
    isempty(components) && throw(ArgumentError(
        "a coaxial blueprint requires at least one retained terminal",
    ))
    starts = Int[1]
    positions = [first(components).conductor.position]
    @inbounds for index in 2:length(components)
        position = components[index].conductor.position
        if !DataModel.same_radial_position(position, last(positions))
            any(reference -> DataModel.same_radial_position(position, reference), positions) &&
                throw(ArgumentError(
                    "a concentric assembly cannot reappear after another assembly",
                ))
            push!(starts, index)
            push!(positions, position)
        end
    end
    return UnitRange{Int}[start:(index == length(starts) ? length(components) :
                                 starts[index + 1] - 1)
                          for (index, start) in pairs(starts)]
end

function validate(blueprint::CableBlueprint)
    isempty(blueprint.cable_id) && throw(ArgumentError(
        "CableBlueprint.cable_id cannot be empty"
    ))
    count = length(blueprint.conductors)
    count > 0 || throw(ArgumentError(
        "CableBlueprint.conductors must contain at least one conductor"
    ))
    length(blueprint.dielectric_ranges) == count || throw(DimensionMismatch(
        "CableBlueprint.dielectric_ranges must contain one range per conductor; " *
        "received $(length(blueprint.dielectric_ranges)) ranges for $count conductors",
    ))
    isempty(blueprint.assembly_ranges) && throw(ArgumentError(
        "CableBlueprint.assembly_ranges must contain at least one assembly",
    ))
    collect(Iterators.flatten(blueprint.assembly_ranges)) == collect(1:count) ||
        throw(DimensionMismatch(
            "CableBlueprint.assembly_ranges must partition conductor indices 1:$count in order",
        ))
    layer_count = length(blueprint.dielectrics)
    collected_layers = collect(Iterators.flatten(blueprint.dielectric_ranges))
    collected_layers == collect(1:layer_count) || throw(DimensionMismatch(
        "CableBlueprint.dielectric_ranges must partition dielectric indices " *
        "1:$layer_count in order",
    ))
    @inbounds for (index, conductor) in pairs(blueprint.conductors)
        isempty(String(conductor.terminal)) && throw(ArgumentError(
            "CableBlueprint.conductors[$index].terminal cannot be empty"
        ))
        conductor.assembly in eachindex(blueprint.assembly_ranges) || throw(DomainError(
            conductor.assembly,
            "CableBlueprint.conductors[$index].assembly must index assembly_ranges"
        ))
        index in blueprint.assembly_ranges[conductor.assembly] || throw(DimensionMismatch(
            "CableBlueprint.conductors[$index].assembly does not own conductor $index",
        ))
        isfinite(conductor.r_in) && conductor.r_in >= zero(conductor.r_in) ||
            throw(DomainError(
                conductor.r_in,
                "CableBlueprint.conductors[$index].r_in must be nonnegative and finite"
            ))
        isfinite(conductor.r_ex) && conductor.r_ex > conductor.r_in ||
            throw(DomainError(
                conductor.r_ex,
                "CableBlueprint.conductors[$index].r_ex must be finite and greater than r_in"
            ))
        isfinite(conductor.cross_section) &&
        conductor.cross_section > zero(conductor.cross_section) || throw(DomainError(
            conductor.cross_section,
            "CableBlueprint.conductors[$index].cross_section must be positive and finite"
        ))
        conductor.num_wires >= 0 || throw(DomainError(
            conductor.num_wires,
            "CableBlueprint.conductors[$index].num_wires must be nonnegative"
        ))
        isfinite(conductor.num_turns) || throw(DomainError(
            conductor.num_turns,
            "CableBlueprint.conductors[$index].num_turns must be finite"
        ))
        isfinite(conductor.resistance) &&
        conductor.resistance > zero(conductor.resistance) || throw(DomainError(
            conductor.resistance,
            "CableBlueprint.conductors[$index].resistance must be positive and finite"
        ))
        isfinite(conductor.alpha) || throw(DomainError(
            conductor.alpha,
            "CableBlueprint.conductors[$index].alpha must be finite"
        ))
        isfinite(conductor.gmr) && conductor.gmr > zero(conductor.gmr) ||
            throw(DomainError(
                conductor.gmr,
                "CableBlueprint.conductors[$index].gmr must be positive and finite"
            ))
        all(isfinite, conductor.position) || throw(DomainError(
            conductor.position,
            "CableBlueprint.conductors[$index].position must be finite"
        ))
        validate(conductor.material)
        conductor.material.kind === :conductor || throw(ArgumentError(
            "CableBlueprint.conductors[$index].material.kind must be :conductor; " *
            "received $(repr(conductor.material.kind))",
        ))
        if !isempty(conductor.layers)
            first(conductor.layers).r_in==conductor.r_in &&
                last(conductor.layers).r_ex==conductor.r_ex ||
                throw(ArgumentError("physical conductor layers must span the conductor annulus"))
            for (n,layer) in pairs(conductor.layers)
                isfinite(layer.r_in) && isfinite(layer.r_ex) &&
                    0<=layer.r_in<layer.r_ex ||
                    throw(DomainError(layer,"invalid physical conductor-layer radii"))
                layer.material.kind===:conductor ||
                    throw(ArgumentError("physical metal profiles require conductor materials"))
                validate(layer.material)
                n==1 || conductor.layers[n-1].r_ex==layer.r_in ||
                    throw(ArgumentError("bonded metal layers must be contiguous"))
            end
        end
        if conductor.acsr !== nothing
            profile=conductor.acsr
            isfinite(profile.radius) && profile.radius>0 &&
                profile.pitch>0 && !isnan(profile.pitch) ||
                throw(DomainError(profile,"invalid physical ACSR radius or pitch"))
            iszero(conductor.r_in) && isapprox(conductor.r_ex,3profile.radius) &&
                isapprox(conductor.cross_section,7π*profile.radius^2) &&
                isempty(conductor.layers) ||
                throw(ArgumentError("ACSR data must describe the retained seven-wire conductor"))
            for material in (profile.core,profile.strands)
                validate(material)
                material.kind===:conductor ||
                    throw(ArgumentError("ACSR wire materials must be conductors"))
            end
            isapprox(profile.strands.mu_r,one(profile.strands.mu_r)) ||
                throw(ArgumentError("Merkushev's aluminium strands must be nonmagnetic"))
        end
        if conductor.section !== nothing
            section=conductor.section
            all(isfinite,(section.area,section.perimeter)) &&
                section.area>0 && section.perimeter>0 ||
                throw(DomainError(section,"invalid physical section area or perimeter"))
            isapprox(section.area,conductor.cross_section) ||
                throw(ArgumentError("physical section area must match the retained conductor"))
            validate(section.material)
            section.material.kind===:conductor ||
                throw(ArgumentError("a physical homogeneous section needs a conductor material"))
        end
        for layer_index in blueprint.dielectric_ranges[index]
            layer = blueprint.dielectrics[layer_index]
            layer.conductor == index || throw(DimensionMismatch(
                "CableBlueprint.dielectrics[$layer_index].conductor must be $index; " *
                "received $(layer.conductor)",
            ))
            isfinite(layer.r_in) && layer.r_in >= zero(layer.r_in) ||
                throw(DomainError(
                    layer.r_in,
                    "CableBlueprint.dielectrics[$layer_index].r_in must be nonnegative and finite"
                ))
            isfinite(layer.r_ex) && layer.r_ex > layer.r_in ||
                throw(DomainError(
                    layer.r_ex,
                    "CableBlueprint.dielectrics[$layer_index].r_ex must be finite and greater than r_in"
                ))
            validate(layer.material)
            layer.material.kind in (:insulator, :semicon) || throw(ArgumentError(
                "CableBlueprint.dielectrics[$layer_index].material.kind must be " *
                ":insulator or :semicon; received $(repr(layer.material.kind))",
            ))
        end
    end
    if !isempty(blueprint.pipes)
        conductors = blueprint.conductors
        ranges = blueprint.assembly_ranges
        all(pipe -> pipe.conductor in eachindex(conductors), blueprint.pipes) ||
            throw(ArgumentError("common pipe conductor index is outside the blueprint"))
        allunique(pipe.conductor for pipe in blueprint.pipes) ||
            throw(ArgumentError("common pipe conductors must be distinct"))
        issorted([conductors[pipe.conductor].r_in for pipe in blueprint.pipes]) ||
            throw(ArgumentError("common pipes must be ordered from inner to outer"))
        allunique(Iterators.flatten(pipe.children for pipe in blueprint.pipes)) ||
            throw(ArgumentError("common pipe children must be distinct and have one direct parent"))
        radii = [
            isempty(blueprint.dielectric_ranges[last(range)]) ?
            conductors[last(range)].r_ex :
            blueprint.dielectrics[last(blueprint.dielectric_ranges[last(range)])].r_ex
            for range in ranges
        ]
        centres = [conductors[first(range)].position for range in ranges]
        for (index, range) in pairs(ranges), conductor in conductors[range]
            left, right = conductor.position, centres[index]
            scale = max(one(eltype(left)), maximum(abs, (left..., right...)))
            tolerance = sqrt(eps(eltype(left))) * scale
            hypot((left .- right)...) <= tolerance ||
                throw(ArgumentError("a radial assembly must have one axis"))
        end
        for pipe in blueprint.pipes
            wall = conductors[pipe.conductor]
            first(ranges[wall.assembly]) == pipe.conductor && wall.r_in > 0 ||
                throw(ArgumentError("a common pipe must begin its own hollow radial assembly"))
            isempty(pipe.children) && throw(ArgumentError("a common pipe needs direct children"))
            validate(pipe.material)
            pipe.material.kind === :insulator &&
                isapprox(pipe.material.mu_r, one(pipe.material.mu_r)) ||
                throw(ArgumentError("common pipe cavities require nonmagnetic insulation"))
            for child in pipe.children
                child in eachindex(ranges) && child != wall.assembly ||
                    throw(ArgumentError("invalid common pipe child index"))
                distance = hypot((centres[child] .- wall.position)...)
                distance > 0 && distance + radii[child] < wall.r_in ||
                    throw(ArgumentError("a common pipe child must fit inside on a distinct axis"))
            end
            for (k, left) in pairs(pipe.children), right in pipe.children[k+1:end]
                hypot((centres[left] .- centres[right])...) > radii[left] + radii[right] ||
                    throw(ArgumentError("common pipe child assemblies must not overlap"))
            end
        end
        for assembly in eachindex(ranges)
            candidates = [
                index for (index, pipe) in pairs(blueprint.pipes)
                if assembly != conductors[pipe.conductor].assembly &&
                   hypot((centres[assembly] .- conductors[pipe.conductor].position)...) +
                   radii[assembly] < conductors[pipe.conductor].r_in
            ]
            expected = isempty(candidates) ? nothing :
                argmin(i -> conductors[blueprint.pipes[i].conductor].r_in, candidates)
            actual = findfirst(pipe -> assembly in pipe.children, blueprint.pipes)
            actual == expected ||
                throw(ArgumentError("common pipe hierarchy must retain the nearest enclosing wall"))
        end
    end
    return blueprint
end

"""
$(TYPEDSIGNATURES)

Flatten one completed cable design into the frequency-independent numerical
description consumed by the coaxial backend.

# Arguments

- `engine`: Coaxial backend identity.
- `design`: Completed physical cable design.
- `T`: Scalar type used by the numerical payload.

# Returns

- A validated [`CableBlueprint`](@ref) with conductor annuli, physical
  dielectric layers, and assembly partitions.
"""
function flatten(
        ::LineCableModelsCoaxial,
        design::CableDesign,
        ::Type{T}
) where {T <: Real}
    components = DataModel.radial_components(design, T)
    ranges = _assembly_ranges(components)
    assembly_by_conductor = Vector{Int}(undef, length(components))
    @inbounds for (assembly, indices) in pairs(ranges), index in indices

        assembly_by_conductor[index] = assembly
    end

    conductors = Vector{BlueprintConductor{T}}(undef, length(components))
    layer_count = sum(component -> length(component.dielectric.layers), components)
    dielectrics = Vector{BlueprintDielectric{T}}(undef, layer_count)
    dielectric_ranges = Vector{UnitRange{Int}}(undef, length(components))
    layer_index = 0
    @inbounds for (index, component) in pairs(components)
        conductor = component.conductor
        conductors[index] = BlueprintConductor{T}(
            component.name,
            assembly_by_conductor[index],
            conductor.r_in,
            conductor.r_ex,
            conductor.cross_section,
            conductor.num_wires,
            conductor.num_turns,
            conductor.resistance,
            conductor.alpha,
            conductor.gmr,
            conductor.position,
            conductor.material,
            _physical_conductor_layers(design,component.name,conductor,T),
            _single_layer_acsr(design,component.name,conductor,T),
            _homogeneous_section(design,component.name,T)
        )
        first_layer = layer_index + 1
        for layer in component.dielectric.layers
            layer_index += 1
            dielectrics[layer_index] = BlueprintDielectric{T}(
                index,
                layer.r_in,
                layer.r_ex,
                layer.material
            )
        end
        dielectric_ranges[index] = first_layer:layer_index
    end
    return CableBlueprint{T}(
        design.cable_id,
        conductors,
        dielectrics,
        dielectric_ranges,
        ranges,
        _pipe_assemblies(design,conductors,dielectrics,dielectric_ranges,ranges)
    )
end

function flatten(
        engine::LineCableModelsCoaxial,
        design::CableDesign
)
    T = eltype(design)
    return flatten(engine, design, T)
end

"""
$(TYPEDEF)

Store the concrete array representation shared by the coaxial local primitive
impedance and potential-coefficient assemblers.

$(TYPEDFIELDS)
"""
struct LocalCableData{T <: Real}
    "Terminal names in DataModel order."
    terminals::Vector{Symbol}
    "Contiguous conductor-index ranges for concentric assemblies."
    assemblies::Vector{UnitRange{Int}}
    "Source-design index for each concentric assembly."
    assembly_designs::Vector{Int}
    "Assembly-local conductor centres [m]."
    positions::Vector{Tuple{T, T}}
    "Equivalent conductor inner radii [m]."
    r_in::Vector{T}
    "Equivalent conductor outer radii [m]."
    r_ext::Vector{T}
    "Inner radii of conductor-owned dielectric intervals [m]."
    r_ins_in::Vector{T}
    "Outer radii of conductor-owned dielectric intervals [m]."
    r_ins_ext::Vector{T}
    "Conductor resistivities at their material reference temperatures [Ω·m]."
    rho0_cond::Vector{T}
    "Conductor material reference temperatures [°C]."
    T0_cond::Vector{T}
    "Conductor temperature coefficients [1/°C]."
    alpha_cond::Vector{T}
    "Conductor relative permeabilities."
    mu_cond::Vector{T}
    "Equivalent relative permeabilities of conductor-owned dielectric intervals."
    mu_ins::Vector{T}
    "Physical dielectric-layer range owned by each conductor."
    dielectric_ranges::Vector{UnitRange{Int}}
    "Physical dielectric-layer inner radii [m]."
    r_layer_in::Vector{T}
    "Physical dielectric-layer outer radii [m]."
    r_layer_ext::Vector{T}
    "Physical dielectric materials in radial order."
    dielectric_materials::Vector{Material{T}}
    "Indices of dielectric layers classified as insulation."
    insulation_indices::Vector{Int}
    "Indices of dielectric layers classified as semiconducting material."
    semicon_indices::Vector{Int}
    "Common pipes with global conductor and assembly indices."
    pipes::Vector{PipeAssembly{T}}
    "Uncombined physical conductor layers for each retained terminal."
    conductor_layers::Vector{Vector{BlueprintConductorLayer{T}}}
    "Physical one-plus-six-wire data for each matching terminal."
    acsr::Vector{Union{Nothing,BlueprintSingleLayerACSR{T}}}
    "Physical homogeneous sections for shape-dependent scalar formulas."
    sections::Vector{Union{Nothing,BlueprintHomogeneousSection{T}}}
end

function LocalCableData(blueprints::AbstractVector{<:CableBlueprint{T}}) where {T <: Real}
    isempty(blueprints) && throw(ArgumentError(
        "local cable data require at least one blueprint",
    ))
    conductor_count = sum(length, blueprints)
    layer_count = sum(blueprint -> length(blueprint.dielectrics), blueprints)
    assembly_count = sum(blueprint -> length(blueprint.assembly_ranges), blueprints)

    terminals = Vector{Symbol}(undef, conductor_count)
    assemblies = Vector{UnitRange{Int}}(undef, assembly_count)
    assembly_designs = Vector{Int}(undef, assembly_count)
    positions = Vector{Tuple{T, T}}(undef, conductor_count)
    r_in_values = Vector{T}(undef, conductor_count)
    r_ext_values = Vector{T}(undef, conductor_count)
    r_ins_in = Vector{T}(undef, conductor_count)
    r_ins_ext = Vector{T}(undef, conductor_count)
    rho0_cond = Vector{T}(undef, conductor_count)
    T0_cond = Vector{T}(undef, conductor_count)
    alpha_cond = Vector{T}(undef, conductor_count)
    mu_cond = Vector{T}(undef, conductor_count)
    mu_ins = Vector{T}(undef, conductor_count)
    dielectric_ranges = Vector{UnitRange{Int}}(undef, conductor_count)
    r_layer_in = Vector{T}(undef, layer_count)
    r_layer_ext = Vector{T}(undef, layer_count)
    dielectric_materials = Vector{Material{T}}(undef, layer_count)
    insulation_indices = Int[]
    semicon_indices = Int[]
    pipes = PipeAssembly{T}[]
    conductor_layers=Vector{Vector{BlueprintConductorLayer{T}}}(undef,conductor_count)
    acsr=Vector{Union{Nothing,BlueprintSingleLayerACSR{T}}}(undef,conductor_count)
    sections=Vector{Union{Nothing,BlueprintHomogeneousSection{T}}}(undef,conductor_count)
    sizehint!(insulation_indices, layer_count)
    sizehint!(semicon_indices, layer_count)

    conductor_offset = 0
    layer_offset = 0
    assembly_offset = 0
    @inbounds for (design_index, blueprint) in pairs(blueprints)
        for pipe in blueprint.pipes
            push!(pipes,PipeAssembly{T}(
                pipe.conductor+conductor_offset,
                pipe.children .+ assembly_offset,
                pipe.material
            ))
        end
        for local_range in blueprint.assembly_ranges
            assembly_offset += 1
            assemblies[assembly_offset] = (
                first(local_range) + conductor_offset
            ):(last(local_range) + conductor_offset)
            assembly_designs[assembly_offset] = design_index
        end
        for local_index in eachindex(blueprint.conductors)
            index = conductor_offset + local_index
            conductor = blueprint.conductors[local_index]
            terminals[index] = conductor.terminal
            positions[index] = conductor.position
            r_in_values[index] = conductor.r_in
            r_ext_values[index] = conductor.r_ex
            rho0_cond[index] = conductor.material.rho
            T0_cond[index] = conductor.material.T0
            alpha_cond[index] = conductor.material.alpha
            mu_cond[index] = conductor.material.mu_r
            conductor_layers[index]=copy(conductor.layers)
            acsr[index]=conductor.acsr
            sections[index]=conductor.section

            local_layers = blueprint.dielectric_ranges[local_index]
            first_layer = layer_offset + 1
            for local_layer in local_layers
                layer_offset += 1
                layer = blueprint.dielectrics[local_layer]
                r_layer_in[layer_offset] = layer.r_in
                r_layer_ext[layer_offset] = layer.r_ex
                dielectric_materials[layer_offset] = layer.material
                if layer.material.kind === :insulator
                    push!(insulation_indices, layer_offset)
                elseif layer.material.kind === :semicon
                    push!(semicon_indices, layer_offset)
                else
                    throw(ArgumentError(
                        "unsupported coaxial dielectric kind :$(layer.material.kind)",
                    ))
                end
            end
            dielectric_ranges[index] = first_layer:layer_offset
            if isempty(local_layers)
                r_ins_in[index] = conductor.r_ex
                r_ins_ext[index] = conductor.r_ex
                mu_ins[index] = one(T)
            else
                r_ins_in[index] = blueprint.dielectrics[first(local_layers)].r_in
                r_ins_ext[index] = blueprint.dielectrics[last(local_layers)].r_ex
                layers = @view blueprint.dielectrics[local_layers]
                mu_ins[index] = DataModel.equivalent_dielectric_permeability(
                    layers,
                    conductor.num_turns,
                    conductor.r_ex,
                    r_ins_ext[index]
                )
            end
        end
        conductor_offset += length(blueprint.conductors)
    end

    return LocalCableData{T}(
        terminals,
        assemblies,
        assembly_designs,
        positions,
        r_in_values,
        r_ext_values,
        r_ins_in,
        r_ins_ext,
        rho0_cond,
        T0_cond,
        alpha_cond,
        mu_cond,
        mu_ins,
        dielectric_ranges,
        r_layer_in,
        r_layer_ext,
        dielectric_materials,
        insulation_indices,
        semicon_indices,
        pipes,
        conductor_layers,
        acsr,
        sections
    )
end

function LocalCableData(blueprint::CableBlueprint{T}) where {T <: Real}
    LocalCableData(CableBlueprint{T}[blueprint])
end
