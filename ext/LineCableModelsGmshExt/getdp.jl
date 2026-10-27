"""
$(TYPEDSIGNATURES)

Write the native model data of `model` to `path`, with the mesh controls in `options`.
Use 17 significant digits for numbers and braces for numerical lists. Region dimensions
come from [`section`](@ref), repeated regions from `model.strand`, passive annuli from
`model.layer`, and terminal names from the model's labels.

The material data evaluates

```math
\\sigma = \\operatorname{Re}\\kappa, \\qquad
\\varepsilon = \\frac{\\operatorname{Im}\\kappa}{2\\pi f}, \\qquad
\\mu = \\mu_0\\mu_r,
```

Here ``\\kappa`` is admittivity \\[S/m\\], ``f`` is frequency \\[Hz\\],
``\\sigma`` is conductivity \\[S/m\\], and ``\\varepsilon`` is permittivity \\[F/m\\].
Permeability is ``\\mu`` \\[H/m\\], and vacuum permeability is ``\\mu_0`` \\[H/m\\].
The factor ``\\mu_r`` is relative permeability \\[dimensionless\\].

Earth uses the formulation's properties from `model.earth` at each frequency. Air uses the
declared properties of the first earth layer, written once. Both use resistivity ``\\rho``
\\[Ω·m\\] with ``\\sigma = 1/\\rho``, and relative permittivity ``\\varepsilon_r``
\\[dimensionless\\] with ``\\varepsilon = \\varepsilon_0\\varepsilon_r``, where
``\\varepsilon_0`` is vacuum permittivity \\[F/m\\]. Coordinates and lengths are in \\[m\\],
and longitudinal wavenumbers are in \\[1/m\\]. Return `path`.
"""
function Base.write(path::AbstractString, model::FieldModel,
        options::ComputationOptions{GmshBackend})
    system = model.problem.system
    frequencies = model.problem.frequencies
    air = model.problem.earth.layers[1]
    controls = options.overrides
    open(path, "w") do io
        println(io, "// Resolved LineCableModels FEM inputs; solver logic lives in getdp/*.pro")
        for (name, number) in (("NumTerminals", length(model.label)),
            ("NumCables", length(system.designs)),
            ("NumMaterialRegions", length(model.materials)),
            ("CADGeometryTolerance", 1e-8), ("FrequencyCount", length(frequencies)))
            println(io, name, " = ", @sprintf("%.17g", Float64(number)), ";")
        end
        for (name, entries) in (("Frequencies", frequencies),
            ("GammaReValues", real.(model.gamma)), ("GammaImValues", imag.(model.gamma)))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (name, number) in (("PmlSideThicknessFactor", controls.PmlSideThicknessFactor),
            ("PmlSideLayers", controls.PmlSideLayers),
            ("PmlSideGrading", controls.PmlSideGrading),
            ("PmlTopThicknessFactor", controls.PmlTopThicknessFactor),
            ("PmlTopLayers", controls.PmlTopLayers),
            ("PmlTopGrading", controls.PmlTopGrading),
            ("PmlBottomThicknessFactor", controls.PmlBottomThicknessFactor),
            ("PmlBottomLayers", controls.PmlBottomLayers),
            ("PmlBottomGrading", controls.PmlBottomGrading),
            ("DomainSizeFactor", options.domain_size_factor),
            ("MeshSizeFactor", options.mesh_size_factor),
            ("ExteriorMeshSizeFactor", controls.ExteriorMeshSizeFactor),
            ("InterfaceRefinementFactor", controls.InterfaceRefinementFactor),
            ("PmlQuadrangles", controls.PmlQuadrangles),
            ("ThinLayerStrips", controls.ThinLayerStrips),
            ("VolumeQuadrature", controls.VolumeQuadrature),
            ("PhysicalVolumeQuadrature", controls.PhysicalVolumeQuadrature),
            ("FEMMumpsForwardErrorBudget", 0.01),
            ("PmlQuadrature", controls.PmlQuadrature),
            ("ConductorGeometryTolerance", controls.ConductorGeometryTolerance),
            ("ConductorSkinDepthElements", controls.ConductorSkinDepthElements),
            ("ConductorMeshGrowth", controls.ConductorMeshGrowth),
            ("ConductorSkinDepths", controls.ConductorSkinDepths),
            ("ConductorThicknessElements", controls.ConductorThicknessElements))
            println(io, name, " = ", @sprintf("%.17g", Float64(number)), ";")
        end
        radii = [support(design.outer) for design in system.designs]
        for (name, entries) in (("CableX", [at.x for at in system.placements]),
            ("CableY", [at.y for at in system.placements]), ("CableRadius", radii))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        println(io, "NumPhysicalRegions = ", @sprintf("%.17g", length(system.regions)), ";")
        for (i, region) in enumerate(system.regions)
            s = section(region)
            println(io, "RegionAnnularPassive~{", i - 1, "} = ", Int(model.layer[i]), ";")
            for (name, number) in (("RegionDimensionKinds", s.kind),
                ("RegionDimension1", s.a), ("RegionDimension2", s.b),
                ("RegionRepeated", model.strand[i]))
                println(io, name, "~{", i - 1, "}() = {",
                    @sprintf("%.17g", Float64(number)), "};")
            end
        end
        offset = 0
        for (k, design) in enumerate(system.designs)
            offset += length(design.regions)
            println(io, "CableSizingRegions~{", k - 1, "}() = {", offset - 1, "};")
        end
        println(io, "PmlReflection = ", @sprintf("%.17g", options.pml_reflection), ";")
        sampling_tag = ("VOLTAGE_SAMPLING_NODES", TAGS.sampling_nodes)
        for (name, tag) in (("AIR_EM", TAGS.air), ("EARTH_EM", TAGS.earth),
            ("AIR_PML", TAGS.air_pml), ("EARTH_PML", TAGS.earth_pml),
            ("DOMAIN_PML", TAGS.pml), ("OUTBND_EM", TAGS.outer_boundary),
            ("OUTBND_ELE_AIR", TAGS.outer_air_boundary),
            ("OUTBND_ELE_REF", TAGS.outer_earth_boundary),
            ("INTERFACE_AIR_SOIL", TAGS.interface),
            ("INNER_PML_BND", TAGS.pml_inner_boundary),
            ("TERMINAL", TAGS.terminal_base + 1),
            ("TERMINAL_CONTOUR", TAGS.terminal_contour_base + 1),
            ("GAUGE_TREE", TAGS.gauge_tree), sampling_tag,
            ("MEASUREMENT_LINE", TAGS.measurement_line_base + 1))
            println(io, name, " = ", @sprintf("%.17g", tag), ";")
        end
        receivers = [Int(model.air[terminal.cable]) for terminal in system.terminals]
        println(io, "ReceiverInAir() = {", join(receivers, ", "), "};")
        println(io, "CABLE_CONTOUR = ", TAGS.cable_contour_base + 1, ";")
        labels = ["\"" * replace(label, '\\' => "\\\\", '"' => "\\\"") * "\""
                  for label in model.label]
        println(io, "TerminalNames() = Str[", join(labels, ", "), "];")
        tags = [material.tag for material in model.materials]
        conductive = [material.conductor for material in model.materials]
        losses = [any(value -> !iszero(real(value)), material.admittivity)
                  for material in model.materials]
        permeabilities = [vacuum_permeability * material.mu_r
                          for material in model.materials]
        for (name, entries) in (("MaterialRegionTags", tags),
            ("MaterialIsConductor", conductive), ("MaterialHasLoss", losses),
            ("MaterialMu", permeabilities))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (i, material) in enumerate(model.materials)
            epsilon = imag.(material.admittivity) ./ (2π .* frequencies)
            for (name, entries) in
                (("Sigma", real.(material.admittivity)), ("Epsilon", epsilon))
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "Material", name, "_", i, "() = {", join(numbers, ", "), "};")
            end
        end
        earth_sigma = [inv(state.rho) for state in model.earth]
        earth_epsilon = [vacuum_permittivity * state.eps_r for state in model.earth]
        earth_mu = [vacuum_permeability * state.mu_r for state in model.earth]
        for (name, entries) in (("EarthSigma", earth_sigma),
            ("EarthEpsilon", earth_epsilon), ("EarthMu", earth_mu))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (name, number) in (("AirSigma", inv(air.rho)),
            ("AirEpsilon", vacuum_permittivity * air.eps_r),
            ("AirMu", vacuum_permeability * air.mu_r),
            ("Xcenter", model.center), ("Ycenter", 0), ("Zcenter", 0))
            println(io, name, " = ", @sprintf("%.17g", Float64(number)), ";")
        end
        return
    end
    return path
end
