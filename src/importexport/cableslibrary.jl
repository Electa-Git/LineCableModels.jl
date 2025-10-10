"""
$(TYPEDSIGNATURES)

Saves a [`CablesLibrary`](@ref) to a file.
The format is determined by the file extension:
- `.json`: Saves using the custom JSON serialization.
- `.jls`: Saves using Julia native binary serialization.

# Arguments
- `library`: The [`CablesLibrary`](@ref) instance to save.
- `file_name`: The path to the output file (default: "cables_library.json").

# Returns
- The absolute path of the saved file, or `nothing` on failure.
"""
function save(
    library::CablesLibrary;
    file_name::String="cables_library.json",
)::Union{String,Nothing}

    file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)

    _, ext = splitext(file_name)
    ext = lowercase(ext)

    try
        if ext == ".jls"
            return _save_cableslibrary_jls(library, file_name)
        elseif ext == ".json"
            return _save_cableslibrary_json(library, file_name)
        else
            @warn "Unrecognized file extension '$ext' for CablesLibrary. Defaulting to .json format."
            # Ensure filename has .json extension if defaulting
            if ext != ".json"
                file_name = file_name * ".json"
            end
            return _save_cableslibrary_json(library, file_name)
        end
    catch e
        @error "Error saving CablesLibrary to '$(display_path(file_name))': $e"
        showerror(stderr, e, catch_backtrace())
        println(stderr)
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

Saves the [`CablesLibrary`](@ref) using Julia native binary serialization.
This format is generally not portable across Julia versions or machine architectures
but can be faster and preserves exact types.

# Arguments
- `library`: The [`CablesLibrary`](@ref) instance.
- `file_name`: The output file path (should end in `.jls`).

# Returns
- The absolute path of the saved file.
"""
function _save_cableslibrary_jls(library::CablesLibrary, file_name::String)::String
    # Note: Serializing the whole library object directly might be problematic
    # if the library struct itself changes. Serializing the core data (designs) is safer.
    serialize(file_name, library.data)
    @info "Cables library saved using Julia serialization to: $(display_path(file_name))"
    return abspath(file_name)
end

"""
$(TYPEDSIGNATURES)

Saves the [`CablesLibrary`](@ref) to a JSON file using the custom serialization logic.

# Arguments
- `library`: The [`CablesLibrary`](@ref) instance.
- `file_name`: The output file path (should end in `.json`).

# Returns
- The absolute path of the saved file.
"""
function _save_cableslibrary_json(library::CablesLibrary, file_name::String)::String
    # Use the generic _serialize_value, which will delegate to _serialize_obj
    # for the library object, which in turn uses _serializable_fields(::CablesLibrary)
    serialized_library = _serialize_value(library)

    open(file_name, "w") do io
        # Use JSON3.pretty for human-readable output
        # allow_inf=true is needed if Measurements or other fields might contain Inf
        JSON3.pretty(io, serialized_library, allow_inf=true)
    end
    if isfile(file_name)
        @info "Cables library saved to: $(display_path(file_name))"
    end
    return abspath(file_name)
end

"""
$(TYPEDSIGNATURES)

Loads cable designs from a file into an existing [`CablesLibrary`](@ref) object.
Modifies the library in-place.
The format is determined by the file extension:
- `.json`: Loads using the custom JSON deserialization and reconstruction.
- `.jls`: Loads using Julia's native binary deserialization.

# Arguments
- `library`: The [`CablesLibrary`](@ref) instance to populate (modified in-place).
- `file_name`: Path to the file to load (default: "cables_library.json").

# Returns
- The modified [`CablesLibrary`](@ref) instance.
"""
function load!(
    library::CablesLibrary; # Type annotation ensures it's the correct object
    file_name::String="cables_library.json",
)::CablesLibrary # Return the modified library
    if !isfile(file_name)
        throw(ErrorException("Cables library file not found: '$(display_path(file_name))'")) # make caller receive an Exception
    end

    _, ext = splitext(file_name)
    ext = lowercase(ext)

    try
        if ext == ".jls"
            _load_cableslibrary_jls!(library, file_name)
        elseif ext == ".json"
            _load_cableslibrary_json!(library, file_name)
        else
            @warn "Unrecognized file extension '$ext' for CablesLibrary. Attempting to load as .json."
            _load_cableslibrary_json!(library, file_name)
        end
    catch e
        @error "Error loading CablesLibrary from '$(display_path(file_name))': $e"
        showerror(stderr, e, catch_backtrace())
        println(stderr)
        # Optionally clear the library or leave it partially loaded depending on desired robustness
        # empty!(library.data)
    end
    return library # Return the modified library
end

"""
$(TYPEDSIGNATURES)

Loads cable designs from a Julia binary serialization file (`.jls`)
into the provided library object.

# Arguments
- `library`: The [`CablesLibrary`](@ref) instance to modify.
- `file_name`: The path to the `.jls` file.

# Returns
- Nothing. Modifies `library` in-place.
"""
function _load_cableslibrary_jls!(library::CablesLibrary, file_name::String)
    loaded_data = deserialize(file_name)

    if isa(loaded_data, Dict{String,CableDesign})
        # Replace the existing designs
        library.data = loaded_data
        println(
            "Cables library successfully loaded via Julia deserialization from: ",
            display_path(file_name),
        )
    else
        # This indicates the .jls file did not contain the expected dictionary structure
        @error "Invalid data format in '$(display_path(file_name))'. Expected Dict{String, CableDesign}, got $(typeof(loaded_data)). Library not loaded."
        # Ensure library.data exists if it was potentially wiped before load attempt
        if !isdefined(library, :data) || !(library.data isa AbstractDict)
            library.data = Dict{String,CableDesign}()
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Loads cable designs from a JSON file into the provided library object
using the detailed, sequential reconstruction logic.

# Arguments
- `library`: The [`CablesLibrary`](@ref) instance to modify.
- `file_name`: The path to the `.json` file.

# Returns
- Nothing. Modifies `library` in-place.
"""
function _load_cableslibrary_json!(library::CablesLibrary, file_name::String)
    # Ensure library structure is initialized
    if !isdefined(library, :data) || !(library.data isa AbstractDict)
        @warn "Library 'data' field was not initialized or not a Dict. Initializing."
        library.data = Dict{String,CableDesign}()
    else
        # Clear existing designs before loading (common behavior)
        empty!(library.data)
    end

    # Load the entire JSON structure
    json_data = open(file_name, "r") do io
        JSON3.read(io, Dict{String,Any}) # Read the top level as a Dict
    end

    # The JSON might store designs directly under "data" key,
    # or the top level might be the dictionary of designs itself.
    local designs_to_process::Dict
    if haskey(json_data, "data") && json_data["data"] isa AbstractDict
        # Standard case: designs are under the "data" key
        designs_to_process = json_data["data"]
    elseif haskey(json_data, "__julia_type__") &&
           occursin("CablesLibrary", json_data["__julia_type__"]) &&
           haskey(json_data, "data")
        # Case where the entire library object was serialized
        designs_to_process = json_data["data"]
    elseif all(
        v ->
            v isa AbstractDict && haskey(v, "__julia_type__") &&
                occursin("CableDesign", v["__julia_type__"]),
        values(json_data),
    )
        # Fallback: Assume the top-level dict *is* the designs dict
        @info "Assuming top-level JSON object in '$(display_path(file_name))' is the dictionary of cable designs."
        designs_to_process = json_data
    else
        @error "JSON file '$(display_path(file_name))' does not contain a recognizable 'data' dictionary or structure."
        return nothing # Exit loading process
    end

    @info "Loading cable designs from JSON: '$(display_path(file_name))'..."
    num_loaded = 0
    num_failed = 0

    # Process each cable design entry using manual reconstruction
    for (cable_id, design_data) in designs_to_process
        if !(design_data isa AbstractDict)
            @warn "Skipping entry '$cable_id': Invalid data format (expected Dictionary, got $(typeof(design_data)))."
            num_failed += 1
            continue
        end
        try
            # Reconstruct the design using the dedicated function
            reconstructed_design =
                _reconstruct_cabledesign(string(cable_id), design_data)
            # Store the fully reconstructed design in the library
            library.data[string(cable_id)] = reconstructed_design
            num_loaded += 1
        catch e
            num_failed += 1
            @error "Failed to reconstruct cable design '$cable_id': $e"
            # Show stacktrace for detailed debugging, especially for MethodErrors during construction
            showerror(stderr, e, catch_backtrace())
            println(stderr) # Add newline for clarity
        end
    end

    @info "Finished loading from '$(display_path(file_name))'. Successfully loaded $num_loaded cable designs, failed to load $num_failed."
    return nothing
end

"""
$(TYPEDSIGNATURES)

Helper function to reconstruct a [`ConductorGroup`](@ref) or [`InsulatorGroup`](@ref) object with the first layer of the respective [`AbstractCablePart`](@ref). Subsequent layers are added using `add!` methods.

# Arguments
- `layer_data`: Dictionary containing the data for the first layer, parsed from JSON.

# Returns
- A reconstructed [`ConductorGroup`](@ref) object with the initial [`AbstractCablePart`](@ref).

# Throws
- Error if essential data is missing or the layer type is unsupported.
"""
function _reconstruct_partsgroup(layer_data::Dict)
    if !haskey(layer_data, "__julia_type__")
        Base.error("Layer data missing '__julia_type__' key: $layer_data")
    end
    type_str = layer_data["__julia_type__"]
    LayerType = _resolve_type(type_str)

    # Use generic deserialization for the whole layer data first.
    # _deserialize_value now returns Dict{Symbol, Any} for plain dicts
    local deserialized_layer_dict::Dict{Symbol,Any}
    try
        # Temporarily remove type key to avoid recursive loop in _deserialize_value -> _deserialize_obj
        temp_data = filter(p -> p.first != "__julia_type__", layer_data)
        deserialized_layer_dict = _deserialize_value(temp_data) # Should return Dict{Symbol, Any}
    catch e
        # This fallback might not be strictly needed anymore if _deserialize_value is robust,
        # but kept for safety. It also needs to produce Dict{Symbol, Any}.
        @error "Initial deserialization failed for first layer data ($type_str): $e. Trying manual field extraction."
        deserialized_layer_dict = Dict{Symbol,Any}()
        for (k_str, v) in layer_data # k_str is String from JSON parsing
            if k_str != "__julia_type__"
                deserialized_layer_dict[Symbol(k_str)] = _deserialize_value(v) # Deserialize value, use Symbol key
            end
        end
    end

    # Ensure the result is Dict{Symbol, Any}
    if !(deserialized_layer_dict isa Dict{Symbol,Any})
        error(
            "Internal error: deserialized_layer_dict is not Dict{Symbol, Any}, but $(typeof(deserialized_layer_dict))",
        )
    end

    # Extract necessary fields using get with Symbol keys
    radius_in = get_as(deserialized_layer_dict, :radius_in, missing, BASE_FLOAT)
    material_props = get_as(deserialized_layer_dict, :material_props, missing, BASE_FLOAT)
    temperature = get_as(deserialized_layer_dict, :temperature, T₀, BASE_FLOAT)

    # Check for essential properties common to most first layers
    ismissing(radius_in) &&
        Base.error("Missing 'radius_in' for first layer type $LayerType in data: $layer_data")
    ismissing(material_props) && error(
        "Missing 'material_props' for first layer type $LayerType in data: $layer_data",
    )
    !(material_props isa Material) && error(
        "'material_props' did not deserialize to a Material object for first layer type $LayerType. Got: $(typeof(material_props))",
    )


    # Type-specific reconstruction using POSITIONAL constructors + Keywords
    # This requires knowing the exact constructor signatures.
    try

        if LayerType == WireArray
            radius_wire = get_as(deserialized_layer_dict, :radius_wire, missing, BASE_FLOAT)
            num_wires = get_as(deserialized_layer_dict, :num_wires, missing, Int)
            lay_ratio = get_as(deserialized_layer_dict, :lay_ratio, missing, BASE_FLOAT)
            lay_direction = get_as(deserialized_layer_dict, :lay_direction, 1, Int) # Default lay_direction
            # Validate required fields
            any(ismissing, (radius_wire, num_wires, lay_ratio)) && error(
                "Missing required field(s) (radius_wire, num_wires, lay_ratio) for WireArray first layer.",
            )
            # Ensure num_wires is Int
            num_wires_int = isa(num_wires, Int) ? num_wires : Int(num_wires)
            lay_direction_int = isa(lay_direction, Int) ? lay_direction : Int(lay_direction)
            return WireArray(
                radius_in,
                radius_wire,
                num_wires_int,
                lay_ratio,
                material_props;
                temperature=temperature,
                lay_direction=lay_direction_int,
            )
        elseif LayerType == Tubular
            radius_ext = get_as(deserialized_layer_dict, :radius_ext, missing, BASE_FLOAT)
            ismissing(radius_ext) && Base.error("Missing 'radius_ext' for Tubular first layer.")
            return Tubular(
                radius_in, radius_ext, material_props; temperature=temperature)
        elseif LayerType == Strip
            radius_ext = get_as(deserialized_layer_dict, :radius_ext, missing, BASE_FLOAT)
            width = get_as(deserialized_layer_dict, :width, missing, BASE_FLOAT)
            lay_ratio = get_as(deserialized_layer_dict, :lay_ratio, missing, BASE_FLOAT)
            lay_direction = get(deserialized_layer_dict, :lay_direction, 1)
            any(ismissing, (radius_ext, width, lay_ratio)) && error(
                "Missing required field(s) (radius_ext, width, lay_ratio) for Strip first layer.",
            )
            lay_direction_int = isa(lay_direction, Int) ? lay_direction : Int(lay_direction)

            return Strip(
                radius_in,
                radius_ext,
                width,
                lay_ratio,
                material_props;
                temperature=temperature,
                lay_direction=lay_direction_int,
            )
        elseif LayerType == Insulator
            radius_ext = get_as(deserialized_layer_dict, :radius_ext, missing, BASE_FLOAT)
            ismissing(radius_ext) &&
                Base.error("Missing 'radius_ext' for Insulator first layer.")
            return Insulator(
                radius_in,
                radius_ext,
                material_props;
                temperature=temperature,
            )
        elseif LayerType == Semicon
            radius_ext = get_as(deserialized_layer_dict, :radius_ext, missing, BASE_FLOAT)
            ismissing(radius_ext) && Base.error(
                "Missing 'radius_ext' for first layer type $LayerType in data: $layer_data",
            )
            return Semicon(radius_in, radius_ext, material_props; temperature=temperature)
        elseif LayerType == Sector
            params = get_as(deserialized_layer_dict, :params, missing, BASE_FLOAT)
            rotation_angle_deg = get_as(deserialized_layer_dict, :rotation_angle_deg, missing, BASE_FLOAT)

            ismissing(params) && Base.error("Missing 'params' for Sector in data: $layer_data")
            !(params isa SectorParams) && error("'params' did not deserialize to a SectorParams object. Got: $(typeof(params))")
            ismissing(rotation_angle_deg) && Base.error("Missing 'rotation_angle_deg' for Sector in data: $layer_data")

            return Sector(params, rotation_angle_deg, material_props; temperature=temperature)
        elseif LayerType == SectorInsulator
            inner_sector = get_as(deserialized_layer_dict, :inner_sector, missing, BASE_FLOAT)
            thickness = get_as(deserialized_layer_dict, :thickness, missing, BASE_FLOAT)

            ismissing(inner_sector) && Base.error("Missing 'inner_sector' for SectorInsulator in data: $layer_data")
            !(inner_sector isa Sector) && error("'inner_sector' did not deserialize to a Sector object. Got: $(typeof(inner_sector))")
            ismissing(thickness) && Base.error("Missing 'thickness' for SectorInsulator in data: $layer_data")

            return SectorInsulator(inner_sector, thickness, material_props; temperature=temperature)
        else
            Base.error("Unsupported layer type for first layer reconstruction: $LayerType")
        end
    catch e
        @error "Construction failed for first layer of type $LayerType with data: $deserialized_layer_dict. Error: $e"
        rethrow(e)
    end
end

"""
$(TYPEDSIGNATURES)

Reconstructs a complete [`CableDesign`](@ref) object from its dictionary representation (parsed from JSON).
This function handles the sequential process of building cable designs:
 1. Deserialize [`NominalData`](@ref).
 2. Iterate through components.
 3. For each component:
    a. Reconstruct the first layer of the conductor group.
    b. Create the [`ConductorGroup`](@ref) with the first layer.
    c. Add subsequent conductor layers using [`add!`](@ref).
    d. Repeat a-c for the [`InsulatorGroup`](@ref).
    e. Create the [`CableComponent`](@ref).
 4. Create the [`CableDesign`](@ref) with the first component.
 5. Add subsequent components using [`add!`](@ref).

# Arguments
- `cable_id`: The identifier string for the cable design.
- `design_data`: Dictionary containing the data for the cable design.

# Returns
- A fully reconstructed [`CableDesign`](@ref) object.

# Throws
- Error if reconstruction fails at any step.
"""
function _reconstruct_cabledesign(
    cable_id::String,
    design_data::Dict,
)::CableDesign
    @info "Reconstructing CableDesign: $cable_id"

    # 1. Reconstruct NominalData using generic deserialization
    local nominal_data::NominalData
    if haskey(design_data, "nominal_data")
        # Ensure the input to _deserialize_value is the Dict for NominalData
        nominal_data_dict = design_data["nominal_data"]
        if !(nominal_data_dict isa AbstractDict)
            error(
                "Invalid format for 'nominal_data' in $cable_id: Expected Dictionary, got $(typeof(nominal_data_dict))",
            )
        end
        nominal_data_val = _deserialize_value(nominal_data_dict)
        if !(nominal_data_val isa NominalData)
            # This error check relies on _deserialize_value returning the original dict on failure
            error(
                "Field 'nominal_data' did not deserialize to a NominalData object for $cable_id. Got: $(typeof(nominal_data_val))",
            )
        end
        nominal_data = nominal_data_val
        @info "  Reconstructed NominalData"
    else
        @warn "Missing 'nominal_data' for $cable_id. Using default NominalData()."
        nominal_data = NominalData() # Use default if missing
    end

    # 2. Process Components Sequentially
    components_data = get(design_data, "components", [])
    if isempty(components_data) || !(components_data isa AbstractVector)
        Base.error("Missing or invalid 'components' array in design data for $cable_id")
    end

    reconstructed_components = CableComponent[] # Store fully built components


    for (idx, comp_data) in enumerate(components_data)
        if !(comp_data isa AbstractDict)
            @warn "Component data at index $idx for $cable_id is not a dictionary. Skipping."
            continue
        end
        comp_id = get(comp_data, "id", "UNKNOWN_COMPONENT_ID_$idx")
        @info "  Processing Component $idx: $comp_id"

        # --- 2.1 Build Conductor Group ---
        local conductor_group::ConductorGroup
        conductor_group_data = get(comp_data, "conductor_group", Dict())
        cond_layers_data = get(conductor_group_data, "layers", [])

        if isempty(cond_layers_data) || !(cond_layers_data isa AbstractVector)
            Base.error("Component '$comp_id' has missing or invalid conductor group layers.")
        end

        # - Create the FIRST layer object
        # Ensure the input to _reconstruct_partsgroup is the Dict for the layer
        first_layer_dict = cond_layers_data[1]
        if !(first_layer_dict isa AbstractDict)
            error(
                "Invalid format for first conductor layer in component '$comp_id': Expected Dictionary, got $(typeof(first_layer_dict))",
            )
        end
        first_cond_layer = _reconstruct_partsgroup(first_layer_dict)

        # - Initialize ConductorGroup using its constructor with the first layer
        conductor_group = ConductorGroup(first_cond_layer)
        @info "    Created ConductorGroup with first layer: $(typeof(first_cond_layer))"

        # - Add remaining layers using add!
        for i in 2:lastindex(cond_layers_data)
            layer_data = cond_layers_data[i]
            if !(layer_data isa AbstractDict)
                @warn "Conductor layer data at index $i for component $comp_id is not a dictionary. Skipping."
                continue
            end

            # Extract Type and necessary arguments for add!
            LayerType = _resolve_type(layer_data["__julia_type__"])
            material_props = get_as(layer_data, "material_props", missing, BASE_FLOAT)
            material_props isa Material || Base.error("'material_props' must deserialize to Material, got $(typeof(material_props))")

            # Prepare args and kwargs based on LayerType for add!
            args = []
            kwargs = Dict{Symbol,Any}()
            kwargs[:temperature] = get_as(layer_data, "temperature", T₀, BASE_FLOAT)
            if haskey(layer_data, "lay_direction") # Only add if present
                kwargs[:lay_direction] = get_as(layer_data, "lay_direction", 1, Int)
            end

            # Extract type-specific arguments needed by add!
            try
                if LayerType == WireArray
                    radius_wire = get_as(layer_data, "radius_wire", missing, BASE_FLOAT)
                    num_wires = get_as(layer_data, "num_wires", missing, Int)
                    lay_ratio = get_as(layer_data, "lay_ratio", missing, BASE_FLOAT)
                    any(ismissing, (radius_wire, num_wires, lay_ratio)) && error(
                        "Missing required field(s) for WireArray layer $i in $comp_id",
                    )
                    args = [radius_wire, num_wires, lay_ratio, material_props]
                elseif LayerType == Tubular
                    radius_ext = get_as(layer_data, "radius_ext", missing, BASE_FLOAT)
                    ismissing(radius_ext) &&
                        Base.error("Missing 'radius_ext' for Tubular layer $i in $comp_id")
                    args = [radius_ext, material_props]
                elseif LayerType == Strip
                    radius_ext = get_as(layer_data, "radius_ext", missing, BASE_FLOAT)
                    width = get_as(layer_data, "width", missing, BASE_FLOAT)
                    lay_ratio = get_as(layer_data, "lay_ratio", missing, BASE_FLOAT)
                    any(ismissing, (radius_ext, width, lay_ratio)) &&
                        Base.error("Missing required field(s) for Strip layer $i in $comp_id")
                    args = [radius_ext, width, lay_ratio, material_props]
                else
                    Base.error("Unsupported layer type '$LayerType' for add!")
                end

                # Call add! with Type, args..., and kwargs...
                add!(conductor_group, LayerType, args...; kwargs...)
                @info "      Added conductor layer $i: $LayerType"
            catch e
                @error "Failed to add conductor layer $i ($LayerType) to component $comp_id: $e"
                println(stderr, "      Layer Data: $layer_data")
                println(stderr, "      Args: $args")
                println(stderr, "      Kwargs: $kwargs")
                rethrow(e)
            end
        end # End loop for conductor layers

        # --- 2.2 Build Insulator Group (Analogous logic) ---
        local insulator_group::InsulatorGroup
        insulator_group_data = get(comp_data, "insulator_group", Dict())
        insu_layers_data = get(insulator_group_data, "layers", [])

        if isempty(insu_layers_data) || !(insu_layers_data isa AbstractVector)
            Base.error("Component '$comp_id' has missing or invalid insulator group layers.")
        end

        # - Create the FIRST layer object
        first_layer_dict_insu = insu_layers_data[1]
        if !(first_layer_dict_insu isa AbstractDict)
            error(
                "Invalid format for first insulator layer in component '$comp_id': Expected Dictionary, got $(typeof(first_layer_dict_insu))",
            )
        end
        first_insu_layer = _reconstruct_partsgroup(first_layer_dict_insu)

        # - Initialize InsulatorGroup
        insulator_group = InsulatorGroup(first_insu_layer)
        @info "    Created InsulatorGroup with first layer: $(typeof(first_insu_layer))"

        # - Add remaining layers using add!
        for i in 2:lastindex(insu_layers_data)
            layer_data = insu_layers_data[i]
            if !(layer_data isa AbstractDict)
                @warn "Insulator layer data at index $i for component $comp_id is not a dictionary. Skipping."
                continue
            end

            LayerType = _resolve_type(layer_data["__julia_type__"])
            material_props = get_as(layer_data, "material_props", missing, BASE_FLOAT)
            material_props isa Material || Base.error("'material_props' must deserialize to Material, got $(typeof(material_props))")


            args = []
            kwargs = Dict{Symbol,Any}()
            kwargs[:temperature] = get_as(layer_data, "temperature", T₀, BASE_FLOAT)

            try
                # All insulator types (Semicon, Insulator) take radius_ext, material_props
                # for the add! method.
                if LayerType in [Semicon, Insulator]
                    radius_ext = get_as(layer_data, "radius_ext", missing, BASE_FLOAT)
                    ismissing(radius_ext) &&
                        Base.error("Missing 'radius_ext' for $LayerType layer $i in $comp_id")
                    args = [radius_ext, material_props]
                else
                    Base.error("Unsupported layer type '$LayerType' for add!")
                end

                # Call add! with Type, args..., and kwargs...
                add!(insulator_group, LayerType, args...; kwargs...)
                @info "      Added insulator layer $i: $LayerType"
            catch e
                @error "Failed to add insulator layer $i ($LayerType) to component $comp_id: $e"
                println(stderr, "      Layer Data: $layer_data")
                println(stderr, "      Args: $args")
                println(stderr, "      Kwargs: $kwargs")
                rethrow(e)
            end
        end # End loop for insulator layers

        # --- 2.3 Create the CableComponent object ---
        component = CableComponent(comp_id, conductor_group, insulator_group)
        push!(reconstructed_components, component)
        @info "    Created CableComponent: $comp_id"

    end # End loop through components_data

    # 3. Create the final CableDesign object using the first component
    if isempty(reconstructed_components)
        Base.error("Failed to reconstruct any valid components for cable design '$cable_id'")
    end
    # Use the CableDesign constructor which takes the first component
    cable_design =
        CableDesign(cable_id, reconstructed_components[1]; nominal_data=nominal_data)
    @info "  Created initial CableDesign with component: $(reconstructed_components[1].id)"

    # 4. Add remaining components to the design sequentially using add!
    for i in 2:lastindex(reconstructed_components)
        try
            add!(cable_design, reconstructed_components[i])
            @info "  Added component $(reconstructed_components[i].id) to CableDesign '$cable_id'"
        catch e
            @error "Failed to add component '$(reconstructed_components[i].id)' to CableDesign '$cable_id': $e"
            rethrow(e)
        end
    end

    @info "Finished Reconstructing CableDesign: $cable_id"
    return cable_design
end