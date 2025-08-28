"""
$(TYPEDSIGNATURES)

Saves a [`MaterialsLibrary`](@ref) to a JSON file.

# Arguments
- `library`: The [`MaterialsLibrary`](@ref) instance to save.
- `file_name`: The path to the output JSON file (default: "materials_library.json").

# Returns
- The absolute path of the saved file, or `nothing` on failure.
"""
function save(
    library::MaterialsLibrary;
    file_name::String="materials_library.json",
)::Union{String,Nothing}
    # TODO: Add jls serialization to materials library.
    # Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/3
    file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)


    _, ext = splitext(file_name)
    ext = lowercase(ext)
    if ext != ".json"
        @warn "MaterialsLibrary only supports .json saving. Forcing extension for file '$file_name'."
        file_name = first(splitext(file_name)) * ".json"
    end

    try

        return _save_materialslibrary_json(library, file_name)

    catch e
        @error "Error saving MaterialsLibrary to '$(display_path(file_name))': $e"
        showerror(stderr, e, catch_backtrace())
        println(stderr)
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

Internal function to save the [`MaterialsLibrary`](@ref) to JSON.

# Arguments
- `library`: The [`MaterialsLibrary`](@ref) instance.
- `file_name`: The output file path.

# Returns
- The absolute path of the saved file.
"""
function _save_materialslibrary_json(library::MaterialsLibrary, file_name::String)::String
    # Check if the library has the data field initialized correctly
    if !isdefined(library, :data) || !(library.data isa AbstractDict)
        error("MaterialsLibrary does not have a valid 'data' dictionary. Cannot save.")
    end

    # Use the generic _serialize_value, which handles the dictionary and its Material contents
    serialized_library_data = _serialize_value(library) # Serialize the dict directly

    open(file_name, "w") do io
        JSON3.pretty(io, serialized_library_data, allow_inf=true)
    end
    if isfile(file_name)
        @info "Materials library saved to: $(display_path(file_name))"
    end

    return abspath(file_name)
end

"""
$(TYPEDSIGNATURES)

Loads materials from a JSON file into an existing [`MaterialsLibrary`](@ref) object.
Modifies the library in-place.

# Arguments
- `library`: The [`MaterialsLibrary`](@ref) instance to populate (modified in-place).
- `file_name`: Path to the JSON file to load (default: \"materials_library.json\").

# Returns
- The modified [`MaterialsLibrary`](@ref) instance.

# See also
- [`MaterialsLibrary`](@ref)
"""
function load!(
    library::MaterialsLibrary;
    file_name::String="materials_library.json",
)::MaterialsLibrary

    if !isfile(file_name)
        throw(ErrorException("Materials library file not found: '$(display_path(file_name))'")) # make caller receive an Exception 

    end

    # Only JSON format is supported now
    _, ext = splitext(file_name)
    ext = lowercase(ext)
    if ext != ".json"
        @error "MaterialsLibrary loading only supports .json files. Cannot load '$(display_path(file_name))'."
        return library
    end

    try
        _load_materialslibrary_json!(library, file_name)
    catch e
        @error "Error loading MaterialsLibrary from '$(display_path(file_name))': $e"
        showerror(stderr, e, catch_backtrace())
        println(stderr)
        # Optionally clear or leave partially loaded
        # empty!(library)
    end
    return library
end

"""
$(TYPEDSIGNATURES)

Internal function to load materials from JSON into the library.

# Arguments
- `library`: The [`MaterialsLibrary`](@ref) instance to modify.
- `file_name`: The path to the JSON file.

# Returns
- Nothing. Modifies `library` in-place.

# See also
- [`MaterialsLibrary`](@ref)
- [`Material`](@ref)
- [`add!`](@ref)
- [`_deserialize_value`](@ref)
"""
function _load_materialslibrary_json!(library::MaterialsLibrary, file_name::String)
    # Ensure library structure is initialized
    if !isdefined(library, :data) || !(library.data isa AbstractDict)
        @warn "Library 'data' field was not initialized or not a Dict. Initializing."
        library.data = Dict{String,Material}()
    else
        # Clear existing materials before loading
        empty!(library.data)
    end

    # Load and parse the JSON data (expecting a Dict of material_name => material_data)
    json_data = open(file_name, "r") do io
        JSON3.read(io, Dict{String,Any})
    end


    @info "Loading materials from JSON: '$(display_path(file_name))'..."
    num_loaded = 0
    num_failed = 0

    # Process each material entry
    for (name::String, material_data::Any) in json_data
        if !(material_data isa AbstractDict)
            @warn "Skipping material '$name': Invalid data format (expected Dictionary, got $(typeof(material_data)))."
            num_failed += 1
            continue
        end
        try
            # Use the generic _deserialize_value function.
            # It will detect __julia_type__ and call _deserialize_obj for Material.
            deserialized_material = _deserialize_value(material_data)

            # **Crucial Check:** Verify the deserialized object is actually a Material
            if deserialized_material isa Material
                add!(library, name, deserialized_material) # Assumes this function exists
                num_loaded += 1
            else
                # This path is taken if _deserialize_obj failed and returned the original Dict
                @warn "Skipping material '$name': Failed to deserialize into Material object. Data received: $material_data"
                # The error from _deserialize_obj inside _deserialize_value would have already been logged.
                num_failed += 1
            end
        catch e
            # Catch errors that might occur outside _deserialize_value (e.g., in add!)
            num_failed += 1
            @error "Error processing material entry '$name': $e"
            showerror(stderr, e, catch_backtrace())
            println(stderr)
        end
    end

    @info "Finished loading materials from '$(display_path(file_name))'. Successfully loaded $num_loaded materials, failed to load $num_failed."
    return nothing
end