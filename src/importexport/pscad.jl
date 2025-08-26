#=
Generates sequential IDs, used for simulation element identification (e.g., PSCAD).
Starts from 100,000,000 and increments.
=#
let current_id = 100000000
    global _next_id = () -> (id = current_id; current_id += 1; string(id))
end

export_data(format::Symbol, args...; kwargs...) = export_data(Val(format), args...; kwargs...)

"""
$(TYPEDSIGNATURES)

Exports a [`LineCableSystem`](@ref) to a PSCAD-compatible file format.

# Arguments

- `cable_system`: A [`LineCableSystem`](@ref) object representing the cable system to be exported.
- `earth_props`: An [`EarthModel`](@ref) object containing the earth properties.
- `base_freq`: The base frequency \\[Hz\\] used for the PSCAD export.
- `file_name`: The path to the output file (default: "*_export.pscx")

# Returns

- The absolute path of the saved file, or `nothing` on failure.

# Examples

```julia
cable_system = LineCableSystem(...)
earth_model = EarthModel(...)
$(FUNCTIONNAME)(cable_system, earth_model, base_freq=50)
```

# See also

- [`LineCableSystem`](@ref)
"""
function export_data(::Val{:pscad},
    cable_system::LineCableSystem,
    earth_props::EarthModel;
    base_freq=f₀,
    file_name::String="$(cable_system.system_id)_export.pscx",
)::Union{String,Nothing}

    file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)

    # Sets attributes on an existing EzXML.Node from a dictionary.
    function _set_attributes!(element::EzXML.Node, attrs::Dict{String,String})
        # Loop through the dictionary and set each attribute on the element
        for (k, v) in attrs
            element[k] = v
        end
    end

    # Adds <param> child elements to an existing <paramlist> EzXML.Node
    # from a vector of ("name", "value") tuples.
    function _add_params_to_list!(
        list_element::EzXML.Node,
        params::Vector{Tuple{String,String}},
    )
        # Ensure the target element is actually a paramlist for clarity, though not strictly necessary for EzXML
        # if nodename(list_element) != "paramlist"
        #     @warn "Attempting to add params to a non-paramlist node: $(nodename(list_element))"
        # end
        # Loop through the vector and add each parameter as a child element
        for (name, value) in params
            param = addelement!(list_element, "param")
            param["name"] = name
            param["value"] = value
        end
    end

    # --- Initial Setup (Identical to original) ---
    # Local Ref for ID generation ensures it's unique to this function call if nested
    current_id = Ref(100000000)
    _next_id() = string(current_id[] += 1)

    # Formatting function (ensure to_nominal is defined or handle types appropriately)
    format_nominal =
        (X; sigdigits=4, minval=-1e30, maxval=1e30) -> begin

            local_value = round(to_nominal(X), sigdigits=sigdigits)

            local_value = max(min(local_value, maxval), minval)
            if abs(local_value) < eps(Float64)
                local_value = 0.0
            end
            return string(local_value)
        end

    id_map = Dict{String,String}() # Stores IDs needed for linking (Instance IDs in this case)
    doc = XMLDocument()
    project = ElementNode("project")
    setroot!(doc, project)
    project_id = cable_system.system_id

    # --- Project Attributes (Identical) ---
    project["name"] = project_id
    project["version"] = "5.0.2"
    project["schema"] = ""
    project["Target"] = "EMTDC"

    # --- Settings (Use Helper for Params) ---
    settings = addelement!(project, "paramlist")
    settings["name"] = "Settings" # Set name attribute directly as in original
    timestamp = string(round(Int, datetime2unix(now())))
    settings_params = [
        ("creator", "LineCableModels.jl,$timestamp"), ("time_duration", "0.5"),
        ("time_step", "5"), ("sample_step", "250"), ("chatter_threshold", ".001"),
        ("branch_threshold", ".0005"), ("StartType", "0"),
        ("startup_filename", "\$(Namespace).snp"), ("PlotType", "0"),
        ("output_filename", "\$(Namespace).out"), ("SnapType", "0"),
        ("SnapTime", "0.3"), ("snapshot_filename", "\$(Namespace).snp"),
        ("MrunType", "0"), ("Mruns", "1"), ("Scenario", ""), ("Advanced", "14335"),
        ("sparsity_threshold", "200"), ("Options", "16"), ("Build", "18"),
        ("Warn", "0"), ("Check", "0"),
        (
            "description",
            "Created with LineCableModels.jl (https://github.com/Electa-Git/LineCableModels.jl)",
        ),
        ("Debug", "0"),
    ]
    _add_params_to_list!(settings, settings_params) # Use helper to add <param> children

    # --- Empty Elements (Identical) ---
    addelement!(project, "Layers")
    addelement!(project, "List")["classid"] = "Settings"
    addelement!(project, "bookmarks")

    # --- GlobalSubstitutions (Identical Structure) ---
    global_subs = addelement!(project, "GlobalSubstitutions")
    global_subs["name"] = "Default"
    addelement!(global_subs, "List")["classid"] = "Sub"
    addelement!(global_subs, "List")["classid"] = "ValueSet"
    global_pl = addelement!(global_subs, "paramlist") # No name attribute
    # Add the single parameter directly as in original
    global_param = addelement!(global_pl, "param")
    global_param["name"] = "Current"
    global_param["value"] = ""

    # --- Definitions Section (Identical Start) ---
    definitions = addelement!(project, "definitions")

    # --- StationDefn (Use Helpers for Attrs/Params) ---
    station = addelement!(definitions, "Definition")
    station_id = _next_id()
    id_map["DS_Defn"] = station_id # Map Definition ID
    station_attrs = Dict(
        "classid" => "StationDefn", "name" => "DS", "id" => station_id,
        "group" => "", "url" => "", "version" => "", "build" => "",
        "crc" => "-1", "view" => "false",
    )
    _set_attributes!(station, station_attrs) # Use helper

    station_pl = addelement!(station, "paramlist")
    station_pl["name"] = "" # Keep empty name attribute exactly as original
    # Add Description param directly as original
    desc_param_st = addelement!(station_pl, "param")
    desc_param_st["name"] = "Description"
    desc_param_st["value"] = ""

    schematic = addelement!(station, "schematic")
    schematic["classid"] = "StationCanvas"
    schematic_pl = addelement!(schematic, "paramlist") # No name attribute
    schematic_params = [
        ("show_grid", "0"), ("size", "0"), ("orient", "1"), ("show_border", "0"),
        ("monitor_bus_voltage", "0"), ("show_signal", "0"), ("show_virtual", "0"),
        ("show_sequence", "0"), ("auto_sequence", "1"), ("bus_expand_x", "8"),
        ("bus_expand_y", "8"), ("bus_length", "4"),
    ]
    _add_params_to_list!(schematic_pl, schematic_params) # Use helper

    addelement!(schematic, "grouping") # Identical

    # --- Station Schematic: Wire/User Instance for "Main" (Use Helpers) ---
    wire = addelement!(schematic, "Wire")
    wire_id = _next_id()
    wire_attrs = Dict(
        "classid" => "Branch", "id" => wire_id, "name" => "Main", "x" => "180",
        "y" => "180",
        "w" => "66", "h" => "82", "orient" => "0", "disable" => "false",
        "defn" => "Main",
        "recv" => "-1", "send" => "-1", "back" => "-1",
    )
    _set_attributes!(wire, wire_attrs) # Use helper

    # Keep vertex loop identical
    for (x, y) in [(0, 0), (0, 18), (54, 54), (54, 72)]
        vertex = addelement!(wire, "vertex")
        vertex["x"] = string(x)
        vertex["y"] = string(y)
    end

    user = addelement!(wire, "User") # User instance nested in Wire
    user_id = _next_id()
    id_map["Main"] = user_id # Original maps the *instance* ID here for hierarchy link
    user_attrs = Dict(
        "classid" => "UserCmp", "id" => user_id, "name" => "$project_id:Main",
        "x" => "0", "y" => "0", "w" => "0", "h" => "0", "z" => "-1", "orient" => "0",
        "defn" => "$project_id:Main", # Links to definition named "Main" (implicitly in same project)
        "link" => "-1", "q" => "4", "disable" => "false",
    )
    _set_attributes!(user, user_attrs) # Use helper

    user_pl = addelement!(user, "paramlist")
    # Original sets attributes directly on paramlist and adds no <param> children - replicate exactly:
    user_pl["name"] = ""
    user_pl["link"] = "-1"
    user_pl["crc"] = "-1"

    # --- UserCmpDefn "Main" (Use Helpers) ---
    user_cmp = addelement!(definitions, "Definition")
    user_cmp_id = _next_id() # This is the definition ID
    id_map["Main_Defn"] = user_cmp_id # Map Definition ID separately
    user_cmp_attrs = Dict(
        "classid" => "UserCmpDefn", "name" => "Main", "id" => user_cmp_id,
        "group" => "",
        "url" => "", "version" => "", "build" => "", "crc" => "-1", "view" => "false",
        "date" => timestamp,
    )
    _set_attributes!(user_cmp, user_cmp_attrs) # Use helper

    user_cmp_pl = addelement!(user_cmp, "paramlist")
    user_cmp_pl["name"] = "" # Empty name attribute
    # Add Description param directly
    desc_param_ucmp = addelement!(user_cmp_pl, "param")
    desc_param_ucmp["name"] = "Description"
    desc_param_ucmp["value"] = ""

    # Form (Identical)
    form = addelement!(user_cmp, "form")
    form["name"] = ""
    form["w"] = "320"
    form["h"] = "400"
    form["splitter"] = "60"

    # Graphics (Identical Structure)
    graphics = addelement!(user_cmp, "graphics")
    graphics["viewBox"] = "-200 -200 200 200"
    graphics["size"] = "2"

    # Graphics Rectangle (Use Helpers)
    rect = addelement!(graphics, "Gfx")
    rect_id = _next_id()
    rect_attrs = Dict(
        "classid" => "Graphics.Rectangle", "id" => rect_id, "x" => "-36", "y" => "-36",
        "w" => "72", "h" => "72",
    )
    _set_attributes!(rect, rect_attrs) # Use helper
    rect_pl = addelement!(rect, "paramlist") # No name attribute
    rect_params = [
        ("color", "Black"), ("dasharray", "0"), ("thickness", "0"), ("port", ""),
        ("fill_style", "0"), ("fill_fg", "Black"), ("fill_bg", "Black"),
        ("cond", "true"),
    ]
    _add_params_to_list!(rect_pl, rect_params) # Use helper

    # Graphics Text (Use Helpers)
    text = addelement!(graphics, "Gfx")
    text_id = _next_id()
    text_attrs = Dict("classid" => "Graphics.Text", "id" => text_id, "x" => "0", "y" => "0")
    _set_attributes!(text, text_attrs) # Use helper
    text_pl = addelement!(text, "paramlist") # No name attribute
    text_params = [
        ("text", "%:Name"), ("anchor", "0"), ("full_font", "Tahoma, 13world"),
        ("angle", "0"), ("color", "Black"), ("cond", "true"),
    ]
    _add_params_to_list!(text_pl, text_params) # Use helper

    # --- UserCmpDefn "Main" Schematic (Use Helpers) ---
    user_schematic = addelement!(user_cmp, "schematic")
    user_schematic["classid"] = "UserCanvas"
    user_sch_pl = addelement!(user_schematic, "paramlist") # No name attribute
    user_sch_params = [
        ("show_grid", "0"), ("size", "0"), ("orient", "1"), ("show_border", "0"),
        ("monitor_bus_voltage", "0"), ("show_signal", "0"), ("show_virtual", "0"),
        ("show_sequence", "0"), ("auto_sequence", "1"), ("bus_expand_x", "8"),
        ("bus_expand_y", "8"), ("bus_length", "4"), ("show_terminals", "0"),
        ("virtual_filter", ""), ("animation_freq", "500"),
    ]
    _add_params_to_list!(user_sch_pl, user_sch_params) # Use helper

    addelement!(user_schematic, "grouping") # Identical

    # --- UserCmpDefn "Main" Schematic: CableSystem Instance (Use Helpers) ---
    cable = addelement!(user_schematic, "Wire") # Wire instance
    cable_id = _next_id()
    cable_attrs = Dict(
        "classid" => "Cable", "id" => cable_id, "name" => "$project_id:CableSystem",
        "x" => "72", "y" => "36", "w" => "107", "h" => "128", "orient" => "0",
        "disable" => "false", "defn" => "$project_id:CableSystem", # Links to definition named "CableSystem"
        "recv" => "-1", "send" => "-1", "back" => "-1", "crc" => "-1",
    )
    _set_attributes!(cable, cable_attrs) # Use helper

    # Keep vertex loop identical
    for (x, y) in [(0, 0), (0, 18), (54, 54), (54, 72)]
        vertex = addelement!(cable, "vertex")
        vertex["x"] = string(x)
        vertex["y"] = string(y)
    end

    cable_user = addelement!(cable, "User") # User instance nested in Wire
    cable_user_id = _next_id()
    id_map["CableSystem"] = cable_user_id # Original maps this *instance* ID for hierarchy link
    cable_user_attrs = Dict(
        "classid" => "UserCmp", "id" => cable_user_id,
        "name" => "$project_id:CableSystem",
        "x" => "0", "y" => "0", "w" => "0", "h" => "0", "z" => "-1", "orient" => "0",
        "defn" => "$project_id:CableSystem", # Links to definition named "CableSystem"
        "link" => "-1", "q" => "4", "disable" => "false",
    )
    _set_attributes!(cable_user, cable_user_attrs) # Use helper

    cable_pl = addelement!(cable_user, "paramlist")
    # Original sets attributes on paramlist AND adds params - replicate exactly
    cable_pl["name"] = ""
    cable_pl["link"] = "-1"
    cable_pl["crc"] = "-1"
    cable_params = [ # Instance parameters
        ("Name", "LineCableSystem_1"), ("R", "#NaN"), ("X", "#NaN"), ("B", "#NaN"),
        ("Freq", format_nominal(base_freq)),
        ("Length", format_nominal(cable_system.line_length / 1000)), # Assumes field exists
        ("Dim", "0"), ("Mode", "0"), ("CoupleEnab", "0"), ("CoupleName", "row"),
        ("CoupleOffset", "0.0 [m]"), ("CoupleRef", "0"), ("tname", "tandem_segment"),
        ("sfault", "0"), ("linc", "10.0 [km]"), ("steps", "3"), ("gen_cnst", "1"),
        ("const_path", "%TEMP%\\my_constants_file.tlo"), ("Date", timestamp),
    ]
    _add_params_to_list!(cable_pl, cable_params) # Use helper

    # --- RowDefn "CableSystem" (Use Helpers) ---
    row = addelement!(definitions, "Definition")
    row_id = _next_id()
    id_map["CableSystem_Defn"] = row_id # Map definition ID separately
    row_attrs = Dict(
        "id" => row_id, "classid" => "RowDefn", "name" => "CableSystem", "group" => "",
        "url" => "", "version" => "RowDefn", "build" => "RowDefn", "crc" => "-1",
        "key" => "", "view" => "false", "date" => timestamp,
    )
    _set_attributes!(row, row_attrs) # Use helper

    row_pl = addelement!(row, "paramlist") # No name attribute
    row_params = [("Description", ""), ("type", "Cable")]
    _add_params_to_list!(row_pl, row_params) # Use helper

    row_schematic = addelement!(row, "schematic")
    row_schematic["classid"] = "RowCanvas"
    row_sch_pl = addelement!(row_schematic, "paramlist") # No name attribute
    row_sch_params =
        [("show_grid", "0"), ("size", "0"), ("orient", "1"), ("show_border", "0")]
    _add_params_to_list!(row_sch_pl, row_sch_params) # Use helper

    # --- Components in RowDefn "CableSystem" Schematic ---

    # FrePhase Component (Use Helpers)
    fre_phase = addelement!(row_schematic, "User")
    fre_phase_id = _next_id()
    fre_phase_attrs = Dict(
        "id" => fre_phase_id, "name" => "master:Line_FrePhase_Options",
        "classid" => "UserCmp",
        "x" => "576", "y" => "180", "w" => "460", "h" => "236", "z" => "-1",
        "orient" => "0",
        "defn" => "master:Line_FrePhase_Options", "link" => "-1", "q" => "4",
        "disable" => "false",
    )
    _set_attributes!(fre_phase, fre_phase_attrs) # Use helper

    fre_pl = addelement!(fre_phase, "paramlist")
    # Original sets crc attribute only on paramlist, replicate exactly
    fre_pl["crc"] = "-1"
    fre_params = [ # Actual params
        ("Interp1", "1"), ("Output", "0"), ("Inflen", "0"), ("FS", "0.5"),
        ("FE", "1.0E6"),
        ("Numf", "100"), ("YMaxP", "20"), ("YMaxE", "0.2"), ("AMaxP", "20"),
        ("AMaxE", "0.2"),
        ("MaxRPtol", "2.0e6"), ("W1", "1.0"), ("W2", "1000.0"), ("W3", "1.0"),
        ("CPASS", "0"),
        ("NFP", "1000"), ("FSP", "0.001"), ("FEP", "1000.0"), ("DCenab", "0"),
        ("DCCOR", "1"),
        ("ECLS", "1"), ("shntcab", "1.0E-9"), ("ET_PE", "1E-10"), ("MER_PE", "2"),
        ("MIT_PE", "5"), ("FDIS", "3"), ("enablf", "1"),
    ]
    _add_params_to_list!(fre_pl, fre_params) # Use helper

    addelement!(row_schematic, "grouping") # Identical

    # --- Coaxial Cables Loop (Use Helpers, keep logic identical) ---
    num_cables = cable_system.num_cables
    dx = 400
    for i in 1:num_cables
        cable_position = cable_system.cables[i]
        coax1 = addelement!(row_schematic, "User")
        coax1_id = _next_id()
        coax1_attrs = Dict(
            "classid" => "UserCmp", "name" => "master:Cable_Coax", "id" => coax1_id,
            "x" => "$(234+(i-1)*dx)", "y" => "612", "w" => "311", "h" => "493",
            "z" => "-1",
            "orient" => "0", "defn" => "master:Cable_Coax", "link" => "-1", "q" => "4",
            "disable" => "false",
        )
        _set_attributes!(coax1, coax1_attrs) # Use helper

        coax1_pl = addelement!(coax1, "paramlist")
        # Original sets attributes on paramlist AND adds params - replicate exactly
        coax1_pl["link"] = "-1"
        coax1_pl["name"] = ""
        coax1_pl["crc"] = "-1"

        # --- Parameter Calculation (Identical Logic from Original) ---
        component_ids = collect(keys(cable_position.design_data.components))
        num_cable_parts = length(component_ids)
        if num_cable_parts > 4
            error(
                "Cable $(cable_position.design_data.cable_id) has $num_cable_parts parts, exceeding the limit of 4 (core/sheath/armor/outer).",
            )
        end
        conn = cable_position.conn
        elim1 = length(conn) >= 2 && conn[2] == 0 ? "1" : "0"
        elim2 = length(conn) >= 3 && conn[3] == 0 ? "1" : "0"
        elim3 = length(conn) >= 4 && conn[4] == 0 ? "1" : "0"
        cable_x = cable_position.horz
        cable_y = cable_position.vert

        # Build the parameter list exactly as in the original's logic
        coax1_params_vector = Vector{Tuple{String,String}}() # Renamed variable
        # Base parameters
        push!(coax1_params_vector, ("CABNUM", "$i"))
        push!(coax1_params_vector, ("Name", "$(cable_position.design_data.cable_id)"))
        push!(coax1_params_vector, ("X", format_nominal(cable_x)))
        push!(coax1_params_vector, ("OHC", "$(cable_y < 0 ? 0 : 1)"))
        push!(
            coax1_params_vector,
            ("Y", (cable_y < 0 ? format_nominal(abs(cable_y)) : "0.0")),
        )
        push!(coax1_params_vector, ("Y2", (cable_y > 0 ? format_nominal(cable_y) : "0.0")))
        push!(coax1_params_vector, ("ShuntA", "1.0e-11 [mho/m]"))
        push!(coax1_params_vector, ("FLT", format_nominal(base_freq)))
        push!(coax1_params_vector, ("RorT", "0"))
        push!(coax1_params_vector, ("LL", "$(2*num_cable_parts-1)"))
        push!(coax1_params_vector, ("CROSSBOND", "0"))
        push!(coax1_params_vector, ("GROUPNO", "1"))
        push!(
            coax1_params_vector,
            ("CBC1", "1"),
            ("CBC2", "0"),
            ("CBC3", "0"),
            ("CBC4", "0"),
        )
        push!(coax1_params_vector, ("SHRad", "1"))
        push!(coax1_params_vector, ("LC", "3"))

        # Component parameters (Keep identical logic)
        ω = 2 * π * base_freq
        for (idx, component_id) in enumerate(component_ids)
            component = cable_position.design_data.components[component_id]
            C_eq = component.insulator_group.shunt_capacitance
            G_eq = component.insulator_group.shunt_conductance
            loss_factor = C_eq > 1e-18 ? G_eq / (ω * C_eq) : 0.0 # Avoid NaN/Inf

            sig_digits_props = 6
            max_loss_tangent = 10.0

            if idx == 1 # Core
                push!(coax1_params_vector, ("CONNAM1", uppercasefirst(component.id)))
                push!(
                    coax1_params_vector,
                    ("R1", format_nominal(component.conductor_group.radius_in)),
                )
                push!(
                    coax1_params_vector,
                    ("R2", format_nominal(component.conductor_group.radius_ext)),
                )
                push!(
                    coax1_params_vector,
                    (
                        "RHOC",
                        format_nominal(
                            component.conductor_props.rho,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERMC",
                        format_nominal(
                            component.conductor_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    ("R3", format_nominal(component.insulator_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T3", "0.0000"))
                push!(coax1_params_vector, ("SemiCL", "0"))
                push!(coax1_params_vector, ("SL2", "0.0000"))
                push!(coax1_params_vector, ("SL1", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "EPS1",
                        format_nominal(
                            component.insulator_props.eps_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERM1",
                        format_nominal(
                            component.insulator_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "LT1",
                        format_nominal(
                            loss_factor,
                            sigdigits=sig_digits_props,
                            maxval=max_loss_tangent,
                        ),
                    ),
                )
            elseif idx == 2 # Sheath
                push!(coax1_params_vector, ("CONNAM2", uppercasefirst(component.id)))
                push!(
                    coax1_params_vector,
                    ("R4", format_nominal(component.conductor_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T4", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "RHOS",
                        format_nominal(
                            component.conductor_props.rho,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERMS",
                        format_nominal(
                            component.conductor_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(coax1_params_vector, ("elim1", elim1))
                push!(
                    coax1_params_vector,
                    ("R5", format_nominal(component.insulator_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T5", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "EPS2",
                        format_nominal(
                            component.insulator_props.eps_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERM2",
                        format_nominal(
                            component.insulator_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "LT2",
                        format_nominal(
                            loss_factor,
                            sigdigits=sig_digits_props,
                            maxval=max_loss_tangent,
                        ),
                    ),
                )
            elseif idx == 3 # Armor
                push!(coax1_params_vector, ("CONNAM3", uppercasefirst(component.id)))
                push!(
                    coax1_params_vector,
                    ("R6", format_nominal(component.conductor_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T6", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "RHOA",
                        format_nominal(
                            component.conductor_props.rho,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERMA",
                        format_nominal(
                            component.conductor_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(coax1_params_vector, ("elim2", elim2))
                push!(
                    coax1_params_vector,
                    ("R7", format_nominal(component.insulator_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T7", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "EPS3",
                        format_nominal(
                            component.insulator_props.eps_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERM3",
                        format_nominal(
                            component.insulator_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "LT3",
                        format_nominal(
                            loss_factor,
                            sigdigits=sig_digits_props,
                            maxval=max_loss_tangent,
                        ),
                    ),
                )
            elseif idx == 4 # Outer
                push!(coax1_params_vector, ("CONNAM4", uppercasefirst(component.id)))
                push!(
                    coax1_params_vector,
                    ("R8", format_nominal(component.conductor_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T8", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "RHOO",
                        format_nominal(
                            component.conductor_props.rho,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERMO",
                        format_nominal(
                            component.conductor_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(coax1_params_vector, ("elim3", elim3))
                push!(
                    coax1_params_vector,
                    ("R9", format_nominal(component.insulator_group.radius_ext)),
                )
                push!(coax1_params_vector, ("T9", "0.0000"))
                push!(
                    coax1_params_vector,
                    (
                        "EPS4",
                        format_nominal(
                            component.insulator_props.eps_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "PERM4",
                        format_nominal(
                            component.insulator_props.mu_r,
                            sigdigits=sig_digits_props,
                        ),
                    ),
                )
                push!(
                    coax1_params_vector,
                    (
                        "LT4",
                        format_nominal(
                            loss_factor,
                            sigdigits=sig_digits_props,
                            maxval=max_loss_tangent,
                        ),
                    ),
                )
            end
        end

        # Default empty values (Keep identical logic)
        if num_cable_parts < 2
            append!(
                coax1_params_vector,
                [
                    ("CONNAM2", "none"),
                    ("R4", "0.0"),
                    ("T4", "0.0000"),
                    ("RHOS", "0.0"),
                    ("PERMS", "0.0"),
                    ("elim1", "0"),
                    ("R5", "0.0"),
                    ("T5", "0.0000"),
                    ("EPS2", "0.0"),
                    ("PERM2", "0.0"),
                    ("LT2", "0.0000"),
                ],
            )
        end
        if num_cable_parts < 3
            append!(
                coax1_params_vector,
                [
                    ("CONNAM3", "none"),
                    ("R6", "0.0"),
                    ("T6", "0.0000"),
                    ("RHOA", "0.0"),
                    ("PERMA", "0.0"),
                    ("elim2", "0"),
                    ("R7", "0.0"),
                    ("T7", "0.0000"),
                    ("EPS3", "0.0"),
                    ("PERM3", "0.0"),
                    ("LT3", "0.0000"),
                ],
            )
        end
        if num_cable_parts < 4
            append!(
                coax1_params_vector,
                [
                    ("CONNAM4", "none"),
                    ("R8", "0.0"),
                    ("T8", "0.0000"),
                    ("RHOO", "0.0"),
                    ("PERMO", "0.0"),
                    ("elim3", "0"),
                    ("R9", "0.0"),
                    ("T9", "0.0000"),
                    ("EPS4", "0.0"),
                    ("PERM4", "0.0"),
                    ("LT4", "0.0000"),
                ],
            )
        end

        # Add all collected parameters to the paramlist created earlier
        _add_params_to_list!(coax1_pl, coax1_params_vector) # Use helper

    end # End Coax cable loop

    # --- Line_Ground Component (Use Helpers) ---
    ground = addelement!(row_schematic, "User")
    ground_id = _next_id()
    ground_attrs = Dict(
        "classid" => "UserCmp", "name" => "master:Line_Ground", "id" => ground_id,
        "x" => "504", "y" => "288", "w" => "793", "h" => "88", "z" => "-1",
        "orient" => "0",
        "defn" => "master:Line_Ground", "link" => "-1", "q" => "4", "disable" => "false",
    )
    _set_attributes!(ground, ground_attrs) # Use helper

    ground_pl = addelement!(ground, "paramlist")
    # Original sets attributes on paramlist AND adds params - replicate exactly
    ground_pl["link"] = "-1"
    ground_pl["name"] = ""
    ground_pl["crc"] = "-1"

    earth_layer = earth_props.layers[end]
    ground_params_vector = [ # Renamed variable
        ("EarthForm2", "0"), ("EarthForm", "3"), ("EarthForm3", "2"), ("GrRho", "0"),
        ("GRRES", format_nominal(earth_layer.base_rho_g)),
        ("GPERM", format_nominal(earth_layer.base_mur_g)),
        ("K0", "0.001"), ("K1", "0.01"), ("alpha", "0.7"),
        ("GRP", format_nominal(earth_layer.base_epsr_g)),
    ]
    _add_params_to_list!(ground_pl, ground_params_vector) # Use helper

    # --- Resource List and Hierarchy (Identical Nested Structure from Original) ---
    addelement!(project, "List")["classid"] = "Resource"

    hierarchy = addelement!(project, "hierarchy")
    # Nested calls exactly as in the original, linking to INSTANCE IDs from id_map
    call1 = addelement!(hierarchy, "call")
    # The link should be to the Station Definition ID, not an instance
    call1["link"] = id_map["DS_Defn"] # Corrected link
    call1["name"] = "$project_id:DS"
    call1["z"] = "-1"
    call1["view"] = "false"
    call1["instance"] = "0"

    call2 = addelement!(call1, "call")
    call2["link"] = id_map["Main"] # Links to Main User INSTANCE ID (as per original id_map usage)
    call2["name"] = "$project_id:Main"
    call2["z"] = "-1"
    call2["view"] = "false"
    call2["instance"] = "0"

    call3 = addelement!(call2, "call")
    call3["link"] = id_map["CableSystem"] # Links to CableSystem User INSTANCE ID (as per original id_map usage)
    call3["name"] = "$project_id:CableSystem"
    call3["z"] = "-1"
    call3["view"] = "true"
    call3["instance"] = "0"

    try
        # Use pretty print option for debugging comparisons if needed
        # open(filename, "w") do io; prettyprint(io, doc); end
        write(file_name, doc) # Standard write
        if isfile(file_name)
            @info "PSCAD file saved to: $(display_path(file_name))"
        end
        return file_name
    catch e
        @error "Failed to write PSCAD file '$(display_path(file_name))': $(e)"
        isa(e, SystemError) && println("SystemError details: ", e.extrainfo)
        return nothing
        rethrow(e) # Rethrow to indicate failure clearly
    end
end