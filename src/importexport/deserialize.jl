"""
$(TYPEDSIGNATURES)

Resolves a fully qualified type name string (e.g., \"Module.Type\")
into a Julia `Type` object. Assumes the type is loaded in the current environment.

# Arguments
- `type_str`: The string representation of the type.

# Returns
- The corresponding Julia `Type` object.

# Throws
- `Error` if the type cannot be resolved.
"""
function _resolve_type(type_str::String)
    try
        return Core.eval(@__MODULE__, Meta.parse(type_str))
    catch e
        @error "Could not resolve type '$type_str'. Ensure module structure is correct and type is loaded in Main."
        rethrow(e)
    end
end

"""
$(TYPEDSIGNATURES)

Deserializes a value from its JSON representation back into a Julia value.
Handles special type markers for `Measurements`, `Inf`/`NaN`, and custom structs
identified by `__julia_type__`. Ensures plain dictionaries use Symbol keys.

# Arguments
- `value`: The JSON-parsed value (Dict, Vector, Number, String, Bool, Nothing).

# Returns
- The deserialized Julia value.
"""
function _deserialize_value(value)
    if value isa Dict
        # Check for special type markers first
        if haskey(value, "__type__")
            type_marker = value["__type__"]
            if type_marker == "Measurement"
                # Reconstruct Measurement
                uncval = get_as(value, "uncertainty", nothing, Measurement)
                if isa(uncval, Measurement)
                    return uncval
                else
                    @warn "Could not reconstruct Measurement from input: value=$(typeof(get_as(value, "value", nothing, BASE_FLOAT))), uncertainty=$(typeof(get_as(value, "uncertainty", nothing, BASE_FLOAT))). Returning original Dict."
                    return value # Return original dict if parts are invalid
                end

            elseif type_marker == "SpecialFloat"
                # Reconstruct Inf/NaN
                val_str = get(value, "value", "")
                if val_str == "Inf"
                    return Inf
                end
                if val_str == "-Inf"
                    return -Inf
                end
                if val_str == "NaN"
                    return NaN
                end
                @warn "Unknown SpecialFloat value: '$val_str'. Returning original Dict."
                return value

            elseif type_marker == "Float"
                return get_as(value, "value", nothing, BASE_FLOAT)

            elseif type_marker == "Int"
                return get_as(value, "value", nothing, Int)

            elseif type_marker == "Complex"
                return get_as(value, "value", nothing, Complex)

            else
                @warn "Unknown __type__ marker: '$type_marker'. Processing as regular dictionary."
                # Fall through to regular dictionary processing
            end
        end

        # Check for Julia object marker
        if haskey(value, "__julia_type__")
            type_str = value["__julia_type__"]
            try
                T = _resolve_type(type_str)
                # Delegate object construction to _deserialize_obj
                return _deserialize_obj(value, T)
            catch e
                # Catch errors specifically from _deserialize_obj or _resolve_type
                @error "Failed to resolve or deserialize type '$type_str': $e. Returning original Dict."
                showerror(stderr, e, catch_backtrace())
                println(stderr)
                return value # Return original dict on error
            end
        end
        return Dict(Symbol(k) => _deserialize_value(v) for (k, v) in value)

    elseif value isa Vector
        # Recursively deserialize array elements
        return [_deserialize_value(v) for v in value]

    else
        # Basic JSON types (Number, String, Bool, Nothing) pass through
        return value

    end
end

"""
$(TYPEDSIGNATURES)

Retrieves a value from a dictionary by key, deserializes it, and coerces it to the specified type `T`. Returns a default value if the key is missing or the value is `missing`.

# Arguments

- `d`: The dictionary to query.
- `key`: The key to look up (symbol or string).
- `default`: The value to return if the key is not present.
- `T`: The target type for coercion.

# Returns

- The value associated with `key` in `d`, deserialized and coerced to type `T`, or `default` if the key is missing or the value is `missing`.

# Examples

```julia
result = $(FUNCTIONNAME)(Dict(:a => 1), :a, 0, Int)  # Returns 1
result = $(FUNCTIONNAME)(Dict(), :b, 42, Int)        # Returns 42
```
"""
get_as(d::AbstractDict, key::Union{Symbol,AbstractString}, default, ::Type{T}) where {T} = begin
    v = get(d, key, default)
    v === missing ? missing : coerce_to_T(_deserialize_value(v), T)
end

"""
$(TYPEDSIGNATURES)

Deserializes a dictionary (parsed from JSON) into a Julia object of type `T`.
Attempts keyword constructor first, then falls back to positional constructor
if the keyword attempt fails with a specific `MethodError`.

# Arguments
- `dict`: Dictionary containing the serialized object data. Keys should match field names.
- `T`: The target Julia `Type` to instantiate.

# Returns
- An instance of type `T`.

# Throws
- `Error` if construction fails by both methods.
"""
function _deserialize_obj(dict::Dict, ::Type{T}) where {T}
    # Prepare a dictionary mapping field symbols to deserialized values
    deserialized_fields = Dict{Symbol,Any}()
    for (key_str, val) in dict
        # Skip metadata keys
        if key_str == "__julia_type__" || key_str == "__type__"
            continue
        end
        key_sym = Symbol(key_str)
        # Ensure value is deserialized before storing
        deserialized_fields[key_sym] = _deserialize_value(val)

    end

    # --- Attempt 1: Keyword Constructor ---
    try
        # Convert Dict{Symbol, Any} to pairs for keyword constructor T(; pairs...)
        # Ensure kwargs only contain keys that are valid fieldnames for T
        # This prevents errors if extra keys were present in JSON
        valid_keys = fieldnames(T)
        kwargs = pairs(filter(p -> p.first in valid_keys, deserialized_fields))

        # @info "Attempting keyword construction for $T with kwargs: $(collect(kwargs))" # Debug logging
        if !isempty(kwargs) || hasmethod(T, Tuple{}, Symbol[]) # Check if kw constructor exists or if kwargs are empty
            return T(; kwargs...)
        else
            # If no kwargs and no zero-arg kw constructor, trigger fallback
            error(
                "No keyword arguments provided and no zero-argument keyword constructor found for $T.",
            )
        end
    catch e
        # Check if the error is specifically a MethodError for the keyword call
        is_kw_meth_error =
            e isa MethodError && (e.f === Core.kwcall || (e.f === T && isempty(e.args))) # Check for kwcall or zero-arg method error

        if is_kw_meth_error
            # @info "Keyword construction failed for $T (as expected for types without kw constructor). Trying positional." # Debug logging
            # Fall through to positional attempt
        else
            # Different error during keyword construction (e.g., type mismatch inside constructor)
            @error "Keyword construction failed for type $T with unexpected error: $e"
            println(stderr, "Input dictionary: $dict")
            println(stderr, "Deserialized fields (kwargs used): $(deserialized_fields)")
            rethrow(e) # Rethrow unexpected errors
        end
    end

    # --- Attempt 2: Positional Constructor (Fallback) ---
    # @info "Attempting positional construction for $T" # Debug logging
    fields_in_order = fieldnames(T)
    positional_args = []

    try
        # Check if the number of deserialized fields matches the number of struct fields
        # This is a basic check for suitability of positional constructor
        # It might be too strict if optional fields were omitted in JSON for keyword constructor types
        # but for true positional types, all fields should generally be present.
        # if length(deserialized_fields) != length(fields_in_order)
        #     error("Number of fields in JSON ($(length(deserialized_fields))) does not match number of fields in struct $T ($(length(fields_in_order))). Cannot use positional constructor.")
        # end

        for field_sym in fields_in_order
            if haskey(deserialized_fields, field_sym)
                push!(positional_args, deserialized_fields[field_sym])
            else
                # If a field is missing, positional construction will fail.
                error(
                    "Cannot attempt positional construction for $T: Missing required field '$field_sym' in input data.",
                )
            end
        end

        # @info "Positional args for $T: $positional_args" # Debug logging
        return T(positional_args...)
    catch e
        # Catch errors during positional construction (e.g., MethodError, TypeError)
        @error "Positional construction failed for type $T with args: $positional_args. Error: $e"
        println(stderr, "Input dictionary: $dict")
        println(
            stderr,
            "Deserialized fields used for positional args: $(deserialized_fields)",
        )
        # Check argument count mismatch again, although the loop above should ensure it if no error occurred there
        if length(positional_args) != length(fields_in_order)
            println(
                stderr,
                "Mismatch between number of args provided ($(length(positional_args))) and fields expected ($(length(fields_in_order))).",
            )
        end
        # Rethrow the error after providing context. This indicates neither method worked.
        rethrow(e)
    end

    # This line should ideally not be reached
    error(
        "Failed to construct object of type $T using both keyword and positional methods.",
    )
end