"Tags of the physical groups in the finite-element model."
const TAGS = (air = 1_001,
    earth = 1_002,
    air_pml = 1_003,
    earth_pml = 1_004,
    pml = 1_005,
    outer_boundary = 2_001,
    interface = 2_002,
    pml_inner_boundary = 2_003,
    outer_air_boundary = 2_004,
    outer_earth_boundary = 2_005,
    terminal_base = 3_000,
    terminal_contour_base = 4_000,
    cable_contour_base = 5_000,
    conductors = 6_001,
    passive = 6_002,
    field_maps = 6_003,
    measurement_line_base = 7_000,
    gauge_tree = 9_001,
    sampling_nodes = 9_002,
    material_base = 10_000)

"""
$(TYPEDEF)

Store the terminal labels, material groups, repeated regions, cable media, and geometric
scales of a finite-element computation.

Resistivity `rho` is in \\[Ω·m\\]. Relative permittivity `eps_r` and permeability `mu_r`
are in \\[dimensionless\\].

$(TYPEDFIELDS)
"""
struct FieldModel{T <: Real, P <: LineParametersProblem{T}}
    "Problem that supplies the cable system, earth, temperature, and frequencies."
    problem::P
    "Label of each terminal row, in system order."
    label::Vector{String}
    """
    The groups of regions that share a material. `conductor` states whether the regions are
    conductive. `mu_r` is the relative permeability \\[dimensionless\\]. `admittivity` holds
    a value at each frequency \\[S/m\\]. `tag` is the physical tag, and `name` is the
    physical name, `LCM/material/<index>/<class>`.
    """
    materials::Vector{@NamedTuple{
        conductor::Bool, mu_r::T, admittivity::Vector{Complex{T}}, tag::Int, name::String}}
    "Material group of each region, in system order."
    material::Vector{Int}
    "Whether each region belongs to a copy of a repeat."
    strand::Vector{Bool}
    "Whether each region is a passive annulus, partitioned into thin-layer strips."
    layer::Vector{Bool}
    "Whether each cable lies in the air."
    air::Vector{Bool}
    "Longitudinal wavenumber Γ at each frequency \\[1/m\\]."
    gamma::Vector{Complex{T}}
    "Properties of the earth at each frequency."
    earth::Vector{EarthLayer{T}}
    "Mean horizontal placement of the cables \\[m\\]."
    center::T
    "Smallest radius of the cable outlines, used to size CAD points \\[m\\]."
    scale::T

    @doc """
    $(TYPEDSIGNATURES)

    Evaluate the material laws at every frequency and group regions with the same material
    type, relative permeability, and admittivity. For a conductive region, evaluate

    ```math
    \\kappa(f) = \\frac{1}{\\rho(\\theta)} + \\mathrm{j} 2\\pi f \\varepsilon_0 \\varepsilon_r,
    ```

    where ``\\kappa`` is admittivity \\[S/m\\], ``\\rho`` is resistivity \\[Ω·m\\],
    ``\\theta`` is temperature \\[°C\\], ``f`` is frequency \\[Hz\\],
    ``\\varepsilon_0`` is vacuum permittivity \\[F/m\\], and ``\\varepsilon_r`` is
    relative permittivity \\[dimensionless\\]. Dielectric regions use the formulation's
    admittivity law for their class.
    Evaluate the formulation's earth law at every frequency.

    # Arguments

    - `problem`: cable system, earth, temperature \\[°C\\], and frequencies \\[Hz\\].
    - `formulation`: material laws and formulation options.

    # Returns

    - A validated [`FieldModel`](@ref).

    # Errors

    - Throws `DomainError` when a conductive region's conductivity is not finite and
      positive.
    - Throws `DomainError` when a cable's outer circle touches the ground surface.
    - Throws `DomainError` when earth conductivity is not finite and nonnegative.
    """
    function FieldModel(problem::LineParametersProblem{T}, formulation::FEM) where {T <:
                                                                                    Real}
        label = reductions(problem.system, formulation.options).label
        sources = [convert(Material{T}, region.source.material)
                   for region in problem.system.regions]
        κ = Vector{Vector{Complex{T}}}(undef, length(sources))
        conductor = fill(false, length(sources))
        for r in findall(source -> source isa Material{T, :conductor}, sources)
            source = sources[r]
            ρ = formulation.resistivity(source, problem.temperature)
            κ[r] = Complex{T}[complex(inv(ρ), 2π * f * vacuum_permittivity * source.eps_r)
                              for f in problem.frequencies]
            conductor[r] = true
        end
        for r in findall(source -> !(source isa Material{T, :conductor}), sources)
            κ[r] = Complex{T}[admittivity(formulation, sources[r], f)
                              for f in problem.frequencies]
        end
        class = map(function (material::Material{<:Real, K}) where {K}
                return string(K)
            end, sources)
        materials = @NamedTuple{
            conductor::Bool, mu_r::T, admittivity::Vector{Complex{T}}, tag::Int,
            name::String}[]
        material = Vector{Int}(undef, length(sources))
        first_region = Int[]
        for r in eachindex(sources)
            mu_r = convert(T, sources[r].mu_r)
            found = findfirst(
                i -> typeof(sources[r]) == typeof(sources[first_region[i]]) &&
                     materials[i].mu_r == mu_r && materials[i].admittivity == κ[r],
                eachindex(materials))
            if found === nothing
                g = length(materials) + 1
                push!(materials,
                    (conductor = conductor[r], mu_r = mu_r,
                        admittivity = κ[r], tag = TAGS.material_base + g,
                        name = "LCM/material/$(lpad(g, 4, '0'))/$(class[r])"))
                push!(first_region, r)
                material[r] = g
            else
                material[r] = found
            end
        end
        strand = fill(false, length(sources))
        offset = 0
        for design in problem.system.designs
            for row in design.repeats, indices in row.copies, r in indices
                strand[offset + r] = true
            end
            offset += length(design.regions)
        end
        layer = fill(false, length(sources))
        for r in findall(
            i -> problem.system.regions[i].shape isa Annulus &&
                 !materials[material[i]].conductor,
            eachindex(problem.system.regions))
            layer[r] = true
        end
        air = [at.y - support(design.outer) > 0
               for (design, at) in zip(problem.system.designs, problem.system.placements)]
        gamma = Vector{Complex{T}}(undef, length(problem.frequencies))
        gamma .= formulation.options.Γ
        earth = [formulation.earth_properties(problem.earth.layers[2], f)
                 for f in problem.frequencies]
        center = convert(T, sum(at.x for at in problem.system.placements) /
                            length(problem.system.placements))
        scale = convert(T, minimum(design -> support(design.outer), problem.system.designs))
        return validate(new{T, typeof(problem)}(
            problem, label, materials, material, strand, layer, air, gamma, earth, center, scale))
    end
end

"""
$(TYPEDSIGNATURES)

Check the regions' and earth's conductivity and the cables' contact with the ground surface.

# Arguments

- `model`: finite-element model with evaluated material groups.

# Returns

- `model`, unchanged.

# Errors

- Throws `DomainError` when a conductive region's conductivity \\[S/m\\] is not finite
  and positive at every frequency.
- Throws `DomainError` when a cable's outer circle touches the ground surface.
- Throws `DomainError` when earth conductivity is not finite and nonnegative.

An error about a cable identifies its index. A conductive region's error also identifies its
tag. An error about earth conductivity identifies the frequency.
"""
function validate(model::FieldModel)
    for (k, state) in enumerate(model.earth)
        σ = inv(state.rho)
        isfinite(σ) && σ >= 0 || throw(DomainError(σ,
            "earth conductivity must be finite and nonnegative at frequency " *
            "$(model.problem.frequencies[k]) Hz"))
    end
    offset = 0
    for (k, design) in enumerate(model.problem.system.designs)
        for r in eachindex(design.regions)
            group = model.materials[model.material[offset + r]]
            if group.conductor &&
               any(κ -> !(isfinite(real(κ)) && real(κ) > 0), group.admittivity)
                throw(DomainError(group.admittivity,
                    "conductive region $(design.regions[r].source.tag) of cable $k " *
                    "($(design.id)) must have finite positive conductivity"))
            end
        end
        at = model.problem.system.placements[k]
        R = support(design.outer)
        abs(at.y) == R && throw(DomainError(at.y,
            "cable $k ($(design.id)) touches the ground surface at outer radius $R m"))
        offset += length(design.regions)
    end
    return model
end
