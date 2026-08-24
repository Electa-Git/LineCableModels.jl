"""
$(TYPEDSIGNATURES)

Write an ATPDraw LCC project for a materialised line-cable system.

The XML contains cable positions, conductor and insulation properties, line
length, system frequency, and the resistivity of the last earth layer.

# Arguments

- `::Val{:atp}`: ATP format selector.
- `cable_system`: Positioned cable designs and line length.
- `earth_props`: Earth model. The last layer supplies `Grnd resis`.

# Keywords

- `base_freq`: ATP system frequency \\[Hz\\]. Default: `50.0`.
- `file_name`: Output path. Relative paths are resolved beside this exporter.
  a supplied basename is prefixed with `cable_system.system_id`. Default:
  `nothing`, which writes `<system_id>_export.xml`.

# Returns

- The absolute output path.

# Units

Units are printed in the XML file according to the ATPDraw specifications:

- Radii (`Rin`, `Rout`, `Rout` of cable): \\[m\\]
- Coordinates (`PosX`, `PosY`): \\[m\\]
- Length (`Length` tag): \\[m\\]
- Frequency (`SysFreq`/`Freq`): \\[Hz\\]
- Resistivity (`rho`, *Grnd resis*): \\[Ω·m\\]
- Relative permittivity (`epsI`) / permeability (`muC`, `muI`): \\[dimensionless\\]
- Shunt capacitance (`Cext`): \\[F/m\\]
- Shunt conductance (`Gext`): \\[S/m\\]

# Notes

- The exporter writes eager equivalent properties at the common material
  reference state. The exporter does not apply operating-temperature correction.
- [`nominal`](@ref) removes uncertainty before numeric values are written.
- LineCableSystem construction, rather than export, checks cable overlap.

  """
function export_data(::Val{:atp},
        cable_system::LineCableSystem,
        earth_props::EarthModel;
        base_freq = 50.0,
        file_name::Union{String, Nothing} = nothing
)::String
    function _set_attributes!(element::EzXML.Node, attrs::Dict)
        for (k, v) in attrs
            element[k] = string(v)
        end
    end
    # --- 1. Setup Constants and Variables ---
    if isnothing(file_name)
        # Derive the name from cable_system when the caller omits it.
        file_name = joinpath(@__DIR__, "$(cable_system.system_id)_export.xml")
    else
        # caller supplied a path/name -> respect directory, but prepend system_id to basename
        requested = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
        if isnothing(cable_system)
            file_name = requested
        else
            dir = dirname(requested)
            base = basename(requested)
            file_name = joinpath(dir, "$(cable_system.system_id)_$base")
        end
    end

    num_phases = length(cable_system.cables)

    # Create XML Structure and LCC Component
    doc = XMLDocument()
    project = ElementNode("project")
    setroot!(doc, project)
    _set_attributes!(
        project,
        Dict("Application" => "ATPDraw", "Version" => "7.3", "VersionXML" => "1")
    )
    header = addelement!(project, "header")
    _set_attributes!(
        header,
        Dict(
            "Timestep" => 1e-6,
            "Tmax" => 0.1,
            "XOPT" => 0,
            "COPT" => 0,
            "SysFreq" => base_freq,
            "TopLeftX" => 200,
            "TopLeftY" => 0
        )
    )
    objects = addelement!(project, "objects")
    variables = addelement!(project, "variables")
    comp = addelement!(objects, "comp")
    _set_attributes!(
        comp,
        Dict(
            "Name" => "LCC",
            "Id" => "$(cable_system.system_id)_1",
            "Capangl" => 90,
            "CapPosX" => -10,
            "CapPosY" => -25,
            "Caption" => ""
        )
    )
    comp_content = addelement!(comp, "comp_content")
    _set_attributes!(
        comp_content,
        Dict(
            "PosX" => 280,
            "PosY" => 360,
            "NumPhases" => num_phases,
            "Icon" => "default",
            "SinglePhaseIcon" => "true"
        )
    )
    for side in ["IN", "OUT"]
        y0 = -20
        for k in 1:num_phases
            y0 += 10
            node = addelement!(comp_content, "node")
            _set_attributes!(
                node,
                Dict(
                    "Name" => "$side$k",
                    "Value" => "C$(k)$(side=="IN" ? "SND" : "RCV")",
                    "UserNamed" => "true",
                    "Kind" => k,
                    "PosX" => side == "IN" ? -20 : 20,
                    "PosY" => y0,
                    "NamePosX" => 0,
                    "NamePosY" => 0
                )
            )
        end
    end

    line_length = nominal(cable_system.line_length)
    soil_rho = nominal(earth_props.layers[end].rho)
    for (name, value) in [
        ("Length", line_length), ("Freq", base_freq), ("Grnd resis", soil_rho)]
        data_node = addelement!(comp_content, "data")
        _set_attributes!(data_node, Dict("Name" => name, "Value" => value))
    end

    # Populate the LCC Sub-structure with CORRECTLY Structured Cable Data
    lcc_node = addelement!(comp, "LCC")
    _set_attributes!(
        lcc_node,
        Dict(
            "NumPhases" => num_phases,
            "IconLength" => "true",
            "LineCablePipe" => 2,
            "ModelType" => 1
        )
    )
    cable_header = addelement!(lcc_node, "cable_header")
    _set_attributes!(
        cable_header,
        Dict("InAirGrnd" => 1, "MatrixOutput" => "true", "ExtraCG" => "$(num_phases)")
    )

    for (k, cable) in enumerate(cable_system.cables)
        cable_node = addelement!(cable_header, "cable")

        num_components = length(cable.design_data.components)
        outermost_radius = nominal(cable.design_data.components[end].insulator_group.r_ex)

        _set_attributes!(
            cable_node,
            Dict(
                "NumCond" => num_components,
                "Rout" => outermost_radius,
                "PosX" => nominal(cable.horz),
                "PosY" => nominal(cable.vert)
            )
        )

        for component in cable.design_data.components
            conductor_node = addelement!(cable_node, "conductor")

            cond_group = component.conductor_group
            cond_props = component.conductor_props
            ins_group = component.insulator_group
            ins_props = component.insulator_props

            rho_eq = (cond_props.rho)
            mu_r_cond = (cond_props.mu_r)
            mu_r_ins = (ins_props.mu_r)
            eps_eq = (ins_props.eps_r)

            _set_attributes!(
                conductor_node,
                Dict(
                    "Rin" => nominal(cond_group.r_in),
                    "Rout" => nominal(cond_group.r_ex),
                    "rho" => nominal(rho_eq),
                    "muC" => nominal(mu_r_cond),
                    "muI" => nominal(mu_r_ins),
                    "epsI" => nominal(eps_eq),
                    "Cext" => nominal(ins_group.shunt_capacitance),
                    "Gext" => nominal(ins_group.shunt_conductance)
                )
            )
        end
    end

    # Write the completed XML document.
    _set_attributes!(variables, Dict("NumSim" => 1, "IOPCVP" => 0, "UseParser" => "false"))

    open(file_name, "w") do fid
        prettyprint(fid, doc)
    end
    @info "XML file saved to: $(_display_path(file_name))"
    return file_name
end

"""
$(TYPEDSIGNATURES)

Write frequency-indexed series-impedance and shunt-admittance matrices to an
ATPDraw `ZY` XML file.

Each frequency produces one `<Z>` and one `<Y>` block. Matrix rows use
comma-separated `R+Xi` and `G+Bi` entries.

# Arguments

- `::Val{:atp}`: ATP format selector.
- `line_params`: Frequency-dependent matrices and frequency vector.

# Keywords

- `file_name`: Output path. Default: `ZY_export.xml` beside this exporter.
- `cable_system`: Optional system supplying the line length and output-name
  prefix. Default: `nothing`.

# Returns

- The absolute output path.

# Units

Units are printed in the XML file according to the ATPDraw specifications:

- `freq` (XML `Freq` attribute): \\[Hz\\]
- `Z` and `Y` entries retain the numerical basis stored by `line_params`.
- XML `Length`: cable-system length \\[m\\], or `1.0` when `cable_system` is
  absent.

# Notes

- The exporter does not scale or recompute the matrices.
- [`nominal`](@ref) removes uncertainty before numeric values are written.

  """
function export_data(::Val{:atp},
        line_params::LineParameters;
        file_name::Union{String, Nothing} = nothing,
        cable_system::Union{LineCableSystem, Nothing} = nothing
)::String

    # Resolve final file_name while preserving any user-supplied path.
    if isnothing(file_name)
        # Derive the name from cable_system when the caller omits it.
        if isnothing(cable_system)
            file_name = joinpath(@__DIR__, "ZY_export.xml")
        else
            file_name = joinpath(@__DIR__, "$(cable_system.system_id)_ZY_export.xml")
        end
    else
        # caller supplied a path/name -> respect directory, but prepend system_id to basename if cable_system provided
        requested = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
        if isnothing(cable_system)
            file_name = requested
        else
            dir = dirname(requested)
            base = basename(requested)
            file_name = joinpath(dir, "$(cable_system.system_id)_$base")
        end
    end

    freq = observe(line_params, frequencies)

    @debug ("ZY export called",
        :method => "ZY",
        :cable_system_isnothing => isnothing(cable_system),
        :cable_system_type => (isnothing(cable_system) ? :nothing : typeof(cable_system)),
        :file_name_in => file_name)

    cable_length = isnothing(cable_system) ? 1.0 : nominal(cable_system.line_length)
    atp_format = "G+Bi"
    # file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)

    open(file_name, "w") do fid
        num_phases = size(observe(line_params, Z), 1)
        y_fmt = (atp_format == "C") ? "C" : "G+Bi"

        @printf(fid,
            "<ZY NumPhases=\"%d\" Length=\"%.4f\" ZFmt=\"R+Xi\" YFmt=\"%s\">\n",
            num_phases,
            cable_length,
            y_fmt)

        # --- Z Matrix Printing ---
        for (k, freq_val) in enumerate(freq)
            @printf(fid, "  <Z Freq=\"%.16E\">\n", nominal(freq_val))
            for i in 1:num_phases
                row_str = join(
                    [@sprintf("%.16E%+.16Ei",
                         nominal(real(observe(line_params, Z, i, j, k))),
                         nominal(imag(observe(line_params, Z, i, j, k))))
                     for j in 1:num_phases],
                    ","
                )
                println(fid, row_str)
            end
            @printf(fid, "  </Z>\n")
        end

        # --- Y Matrix Printing ---
        if atp_format == "C"
            freq1 = nominal(freq[1])
            @printf(fid, "  <Y Freq=\"%.16E\">\n", freq1)
            for i in 1:num_phases
                row_str = join(
                    [@sprintf("%.16E",
                         nominal(observe(line_params, C, i, j, 1)))
                     for j in 1:num_phases],
                    ","
                )
                println(fid, row_str)
            end
            @printf(fid, "  </Y>\n")
        else # Case for "G+Bi"
            for (k, freq_val) in enumerate(freq)
                @printf(fid, "  <Y Freq=\"%.16E\">\n", nominal(freq_val))
                for i in 1:num_phases
                    row_str = join(
                        [@sprintf("%.16E%+.16Ei",
                             nominal(real(observe(line_params, Y, i, j, k))),
                             nominal(imag(observe(line_params, Y, i, j, k))))
                         for j in 1:num_phases],
                        ","
                    )
                    println(fid, row_str)
                end
                @printf(fid, "  </Y>\n")
            end
        end

        # --- Footer ---
        println(fid, "</ZY>")
    end
    @info "XML file saved to: $(_display_path(file_name))"
    return file_name
end
