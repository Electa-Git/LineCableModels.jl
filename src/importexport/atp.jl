
"""
$(TYPEDSIGNATURES)

Exports calculated [`LineParameters`](@ref) to an ATP-style XML file.

This function takes all the system information, cables, ground parameters and frequency and assembles the XML to be used in ATPDraw

# Arguments

- `problem`: A [`LineParametersProblem`](@ref) object used to retrieve the frequency vector for the export.
- `file_name`: The path to the output XML file (default: "*_export.xml").

# Returns

- The absolute path of the saved file.
"""
function export_data(::Val{:atp},
    cable_system::LineCableSystem,
    earth_props::EarthModel;
    base_freq=f₀,
    file_name::String="$(cable_system.system_id)_export.xml"
)::Union{String,Nothing}

    function _set_attributes!(element::EzXML.Node, attrs::Dict)
        for (k, v) in attrs
            element[k] = string(v)
        end
    end
    # --- 1. Setup Constants and Variables ---
    file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
    num_phases = length(cable_system.cables)
    
    # Create XML Structure and LCC Component
    doc = XMLDocument()
    project = ElementNode("project")
    setroot!(doc, project)
    _set_attributes!(project, Dict("Application" => "ATPDraw", "Version" => "7.3", "VersionXML" => "1"))
    header = addelement!(project, "header")
    _set_attributes!(header, Dict("Timestep" => 1e-6, "Tmax" => 0.1, "XOPT" => 0, "COPT" => 0, "SysFreq" => base_freq, "TopLeftX" => 200, "TopLeftY" => 0))
    objects = addelement!(project, "objects")
    variables = addelement!(project, "variables")
    comp = addelement!(objects, "comp")
    _set_attributes!(comp, Dict("Name" => "LCC", "Id" => "$(cable_system.system_id)_1", "Capangl" => 90, "CapPosX" => -10, "CapPosY" => -25, "Caption" => ""))
    comp_content = addelement!(comp, "comp_content")
    _set_attributes!(comp_content, Dict("PosX" => 280, "PosY" => 360, "NumPhases" => num_phases, "Icon" => "default", "SinglePhaseIcon" => "true"))
    for side in ["IN", "OUT"]; y0 = -20; for k in 1:num_phases; y0 += 10; node = addelement!(comp_content, "node"); _set_attributes!(node, Dict("Name" => "$side$k", "Value" => "C$(k)$(side=="IN" ? "SND" : "RCV")", "UserNamed" => "true", "Kind" => k, "PosX" => side == "IN" ? -20 : 20, "PosY" => y0, "NamePosX" => 0, "NamePosY" => 0)); end; end
    soil_rho = earth_props.layers[end].base_rho_g
    for (name, value) in [("Length", cable_system.line_length), ("Freq", base_freq), ("Grnd resis", soil_rho)]; data_node = addelement!(comp_content, "data"); _set_attributes!(data_node, Dict("Name" => name, "Value" => value)); end
    
    # Populate the LCC Sub-structure with CORRECTLY Structured Cable Data
    lcc_node = addelement!(comp, "LCC")
    _set_attributes!(lcc_node, Dict("NumPhases" => num_phases, "IconLength" => "true", "LineCablePipe" => 2, "ModelType" => 1))
    cable_header = addelement!(lcc_node, "cable_header")
    _set_attributes!(cable_header, Dict("InAirGrnd" => 1, "MatrixOutput" => "true", "ExtraCG"=>"$(num_phases)"))
    
    for (k, cable) in enumerate(cable_system.cables)
        cable_node = addelement!(cable_header, "cable")
        
        num_components = length(cable.design_data.components)
        outermost_radius = cable.design_data.components[end].insulator_group.radius_ext

        _set_attributes!(cable_node, Dict(
            "NumCond" => num_components,
            "Rout" => outermost_radius,
            "PosX" => cable.horz,
            "PosY" => cable.vert
        ))
        
        for component in cable.design_data.components
            conductor_node = addelement!(cable_node, "conductor")

            cond_group = component.conductor_group
            ins_group = component.insulator_group

            rho_eq = calc_equivalent_rho(cond_group.resistance, cond_group.radius_ext, cond_group.radius_in)
            mu_r_cond = calc_equivalent_mu(cond_group.gmr, cond_group.radius_ext, cond_group.radius_in)
            mu_r_ins = ins_group.layers[1].material_props.mu_r
            eps_eq = calc_equivalent_eps(ins_group.shunt_capacitance, ins_group.radius_in, ins_group.radius_ext)

            _set_attributes!(conductor_node, Dict(
                "Rin" => cond_group.radius_in,
                "Rout" => cond_group.radius_ext,
                "rho" => rho_eq,
                "muC" => mu_r_cond,
                "muI" => mu_r_ins,
                "epsI" => eps_eq,
                "Cext" => ins_group.shunt_capacitance,
                "Gext" => ins_group.shunt_conductance
            ))
        end
    end

    # Finalize and Write to File
    _set_attributes!(variables, Dict("NumSim" => 1, "IOPCVP" => 0, "UseParser" => "false"))

    try
        open(file_name, "w") do fid
            prettyprint(fid, doc)
        end
        @info "XML file saved to: $(_display_path(file_name))"
        return file_name
    catch e
        @error "Failed to write XML file '$file_name'" exception=(e, catch_backtrace())
        return nothing
    end
end


"""
    read_atp_data(file_name::String, cable_system::LineCableSystem)

Reads an ATP `.lis` output file, extracts the Ze and Zi matrices, and dynamically
reorders them to a grouped-by-phase format based on the provided `cable_system`
structure. It correctly handles systems with a variable number of components per cable.

# Arguments
- `file_name`: The path to the `.lis` file.
- `cable_system`: The `LineCableSystem` object corresponding to the data in the file.

# Returns
- `Array{T, 2}`: A 2D complex matrix representing the total reordered series
  impedance `Z = Ze + Zi` for a single frequency.
- `nothing`: If the file cannot be found, parsed, or if the matrix dimensions in the
  file do not match the provided `cable_system` structure.
"""
function read_data(::Val{:atp},
    cable_system::LineCableSystem,
    freq::AbstractFloat;
    file_name::String="$(cable_system.system_id)_1.lis"
    )::Union{Array{COMPLEXSCALAR, 2}, Nothing}
    # --- Inner helper function to parse a matrix block from text lines ---
    function parse_block(block_lines::Vector{String})
        data_lines = filter(line -> !isempty(strip(line)), block_lines)
        if isempty(data_lines) return Matrix{ComplexF64}(undef, 0, 0) end
        matrix_size = length(split(data_lines[1]))
        real_parts = zeros(Float64, matrix_size, matrix_size)
        imag_parts = zeros(Float64, matrix_size, matrix_size)
        row_counter = 1
        for i in 1:2:length(data_lines)
            if i + 1 > length(data_lines) break end
            real_line, imag_line = data_lines[i], data_lines[i+1]
            try
                real_parts[row_counter, :] = [parse(Float64, s) for s in split(real_line)[1:matrix_size]]
                imag_parts[row_counter, :] = [parse(Float64, s) for s in split(imag_line)[1:matrix_size]]
            catch e; @error "Parsing failed" exception=(e, catch_backtrace()); return nothing end
            row_counter += 1
            if row_counter > matrix_size break end
        end
        return real_parts + im * imag_parts
    end

    # --- Main Function Logic ---
    if !isfile(file_name) @error "File not found: $file_name"; return nothing end
    lines = readlines(file_name)
    ze_start_idx = findfirst(occursin.("Earth impedance [Ze]", lines))
    zi_start_idx = findfirst(occursin.("Conductor internal impedance [Zi]", lines))
    if isnothing(ze_start_idx) || isnothing(zi_start_idx) @error "Could not find Ze/Zi headers."; return nothing end
    
    Ze = parse_block(lines[ze_start_idx + 1 : zi_start_idx - 1])
    Zi = parse_block(lines[zi_start_idx + 1 : end])
    if isnothing(Ze) || isnothing(Zi) return nothing end

    # --- DYNAMICALLY GENERATE PERMUTATION INDICES (Numerical Method) ---
    component_counts = [length(c.design_data.components) for c in cable_system.cables]
    total_conductors = sum(component_counts)
    num_phases = length(component_counts)
    max_components = isempty(component_counts) ? 0 : maximum(component_counts)

    if size(Ze, 1) != total_conductors
        @error "Matrix size from file ($(size(Ze,1))x$(size(Ze,1))) does not match total components in cable_system ($total_conductors)."
        return nothing
    end

    num_conductors_per_type = [sum(c >= i for c in component_counts) for i in 1:max_components]
    type_offsets = cumsum([0; num_conductors_per_type[1:end-1]])

    permutation_indices = Int[]
    sizehint!(permutation_indices, total_conductors)
    instance_counters = ones(Int, max_components)
    for phase_idx in 1:num_phases
        for comp_type_idx in 1:component_counts[phase_idx]
            instance = instance_counters[comp_type_idx]
            original_idx = type_offsets[comp_type_idx] + instance
            push!(permutation_indices, original_idx)
            instance_counters[comp_type_idx] += 1
        end
    end
    
    Ze_reordered = Ze[permutation_indices, permutation_indices]
    Zi_reordered = Zi[permutation_indices, permutation_indices]
    
    return Ze_reordered+Zi_reordered
end


"""
$(TYPEDSIGNATURES)

Exports calculated [`LineParameters`](@ref) to an ATP-style XML file.

This function takes the results of a simulation (Z and Y matrices) and writes them
into a structured XML format for use in other programs.

# Arguments

- `line_params`: A [`LineParameters`](@ref) object containing the calculated Z and Y matrices to be exported.
- `cable_system`: A [`LineCableSystem`](@ref) object used for metadata (e.g., the default filename).
- `problem`: A [`LineParametersProblem`](@ref) object used to retrieve the frequency vector for the export.
- `file_name`: The path to the output XML file (default: "*_export.xml").

# Returns

- The absolute path of the saved file.
"""
function export_data(::Val{:atp},
    line_params::LineParameters,
    freq::Vector{BASE_FLOAT};
    file_name::String="ZY_export.xml",
    cable_system::Union{LineCableSystem,Nothing}=nothing
    )::Union{String,Nothing}
    # Construct the file name with system_id if cable_system is provided
    file_name = isnothing(cable_system) ? file_name : "$(cable_system.system_id)_ZY_export.xml"
   
    cable_length = 1.0
    atp_format = "G+Bi"
    file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
    
    open(file_name, "w") do fid
        num_phases = size(line_params.Z, 1)
        y_fmt = (atp_format == "C") ? "C" : "G+Bi"

        @printf(fid, "<ZY NumPhases=\"%d\" Length=\"%.4f\" ZFmt=\"R+Xi\" YFmt=\"%s\">\n", num_phases, cable_length, y_fmt)

        # --- Z Matrix Printing ---
        for (k, freq_val) in enumerate(freq)
            @printf(fid, "  <Z Freq=\"%.16E\">\n", freq_val)
            for i in 1:num_phases
                row_str = join([@sprintf("%.16E%+.16Ei", real(line_params.Z[i, j, k]), imag(line_params.Z[i, j, k])) for j in 1:num_phases], ",")
                println(fid, row_str)
            end
            @printf(fid, "  </Z>\n")
        end

        # --- Y Matrix Printing ---
        if atp_format == "C"
            freq1 = f[1]
            @printf(fid, "  <Y Freq=\"%.16E\">\n", freq1)
            for i in 1:num_phases
                row_str = join([@sprintf("%.16E", imag(line_params.Y[i, j, 1]) / (2 * pi * freq1)) for j in 1:num_phases], ",")
                println(fid, row_str)
            end
            @printf(fid, "  </Y>\n")
        else # Case for "G+Bi"
            for (k, freq_val) in enumerate(freq)
                @printf(fid, "  <Y Freq=\"%.16E\">\n", freq_val)
                for i in 1:num_phases
                    row_str = join([@sprintf("%.16E%+.16Ei", real(line_params.Y[i, j, k]), imag(line_params.Y[i, j, k])) for j in 1:num_phases], ",")
                    println(fid, row_str)
                end
                @printf(fid, "  </Y>\n")
            end
        end

        # --- Footer ---
        println(fid, "</ZY>")
    end
    try
        # Use pretty print option for debugging comparisons if needed
        # open(filename, "w") do io; prettyprint(io, doc); end
        if isfile(file_name)
            @info "XML file saved to: $(file_name)"
        end
        return file_name
    catch e
        @error "Failed to write XML file '$(file_name)': $(e)"
        isa(e, SystemError) && println("SystemError details: ", e.extrainfo)
        return nothing
        rethrow(e) # Rethrow to indicate failure clearly
    end
end