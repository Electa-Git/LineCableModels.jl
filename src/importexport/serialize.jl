"""
$(TYPEDSIGNATURES)

Defines which fields of a given object should be serialized to JSON.
This function acts as a trait. Specific types should overload this method
to customize which fields are needed for reconstruction.

# Arguments

- `obj`: The object whose serializable fields are to be determined.

# Returns

- A tuple of symbols representing the fields of `obj` that should be serialized.

# Methods

$(_CLEANMETHODLIST)
"""
function _serializable_fields end

# Default fallback: Serialize all fields. This might include computed fields
# that are not needed for reconstruction. Overload for specific types.
_serializable_fields(obj::T) where {T} = fieldnames(T)

# Define exactly which fields are needed to reconstruct each object.
# These typically match the constructor arguments or the minimal set
# required by the reconstruction logic (e.g., for groups).

# Core Data Types
_serializable_fields(::Material) = (:rho, :eps_r, :mu_r, :T0, :alpha)
_serializable_fields(::NominalData) = (
    :designation_code,
    :U0,
    :U,
    :conductor_cross_section,
    :screen_cross_section,
    :armor_cross_section,
    :resistance,
    :capacitance,
    :inductance,
)

# Layer Types (Conductor Parts)
_serializable_fields(::WireArray) = (
    :radius_in, # Needed for first layer reconstruction
    :radius_wire,
    :num_wires,
    :lay_ratio,
    :material_props,
    :temperature,
    :lay_direction,
)
_serializable_fields(::Tubular) = (
    :radius_in, # Needed for first layer reconstruction
    :radius_ext,
    :material_props,
    :temperature,
)
_serializable_fields(::Strip) = (
    :radius_in, # Needed for first layer reconstruction
    :radius_ext,
    :width,
    :lay_ratio,
    :material_props,
    :temperature,
    :lay_direction,
)

# Layer Types (Insulator Parts)
_serializable_fields(::Insulator) = (
    :radius_in, # Needed for first layer reconstruction
    :radius_ext,
    :material_props,
    :temperature,
)
_serializable_fields(::Semicon) = (
    :radius_in, # Needed for first layer reconstruction
    :radius_ext,
    :material_props,
    :temperature,
)

# Group Types - Only serialize the layers needed for reconstruction.
_serializable_fields(::ConductorGroup) = (:layers,)
_serializable_fields(::InsulatorGroup) = (:layers,)

# Component & Design Types
_serializable_fields(::CableComponent) = (:id, :conductor_group, :insulator_group)
# For CableDesign, we need components and nominal data. ID is handled as the key.
_serializable_fields(::CableDesign) = (:cable_id, :nominal_data, :components)

# Library Types
_serializable_fields(::CablesLibrary) = (:data,)
_serializable_fields(::MaterialsLibrary) = (:data,)


#=
Serializes a Julia value into a JSON-compatible representation.
Handles special types like Measurements, Inf/NaN, Symbols, and custom structs
using the `_serializable_fields` trait.

# Arguments
- `value`: The Julia value to serialize.

# Returns
- A JSON-compatible representation (Dict, Vector, Number, String, Bool, Nothing).
=#
# Helper: only used in serialization, never leaks to core math.
function _serialize_value(value)

    if isnothing(value)
        return nothing

    elseif value isa Measurements.Measurement
        v = Measurements.value(value)
        u = Measurements.uncertainty(value)
        return Dict(
            "__type__" => "Measurement",
            "value" => _serialize_value(v),
            "uncertainty" => _serialize_value(u),
        )

    elseif value isa Number && !isfinite(value)
        # Inf / -Inf / NaN stay tagged
        local val_str
        if isinf(value)
            val_str = value > 0 ? "Inf" : "-Inf"
        else
            # NaN
            val_str = "NaN"
        end
        return Dict("__type__" => "SpecialFloat", "value" => val_str)

    elseif value isa AbstractFloat
        return Dict("__type__" => "Float", "value" => value)

    elseif value isa Integer
        return Dict("__type__" => "Int", "value" => value)

    elseif value isa Complex
        return Dict("__type__" => "Complex",
            "re" => _serialize_value(real(value)),
            "im" => _serialize_value(imag(value)),
        )

    elseif value isa Number || value isa String || value isa Bool
        return value

    elseif value isa Symbol
        return string(value)

    elseif value isa AbstractDict
        return Dict(string(k) => _serialize_value(v) for (k, v) in value)

    elseif value isa Union{AbstractVector,Tuple}
        return [_serialize_value(v) for v in value]
    else
        !isprimitivetype(typeof(value)) && fieldcount(typeof(value)) > 0
        # Custom structs
        return _serialize_obj(value)
    end
end


"""
$(TYPEDSIGNATURES)

Serializes a Julia value into a JSON-compatible representation.
Handles special types like Measurements, Inf/NaN, Symbols, and custom structs
using the [`_serializable_fields`](@ref) trait.

# Arguments
- `value`: The Julia value to serialize.

# Returns
- A JSON-compatible representation (Dict, Vector, Number, String, Bool, Nothing).
"""
function _serialize_obj(obj)
    T = typeof(obj)
    # Get fully qualified type name (e.g., Main.MyModule.MyType)
    try
        mod = parentmodule(T)
        typeName = nameof(T)
        type_str = string(mod, ".", typeName)

        result = Dict{String,Any}()
        result["__julia_type__"] = type_str

        # Get the fields to serialize using the trait function
        fields_to_include = _serializable_fields(obj)

        # Iterate only through the fields specified by the trait
        for field in fields_to_include
            if hasproperty(obj, field)
                value = getfield(obj, field)
                result[string(field)] = _serialize_value(value) # Recursively serialize
            else
                # This indicates an issue with the _serializable_fields definition for T
                @warn "Field :$field specified by _serializable_fields(::$T) not found in object. Skipping."
            end
        end
        return result
    catch e
        error("Error determining module or type name for object of type $T: $e. Cannot serialize.")
        # Return a representation indicating the error
        return Dict(
            "__error__" => "Serialization failed for type $T",
            "__details__" => string(e),
        )
    end
end