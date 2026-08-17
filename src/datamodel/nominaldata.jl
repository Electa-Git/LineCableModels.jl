"""
$(TYPEDEF)

Store optional cable-datasheet values.

$(TYPEDFIELDS)
"""
struct NominalData{T <: Real}
    designation_code::Union{Nothing, String}
    U0::Union{Nothing, T}
    U::Union{Nothing, T}
    conductor_cross_section::Union{Nothing, T}
    screen_cross_section::Union{Nothing, T}
    armor_cross_section::Union{Nothing, T}
    resistance::Union{Nothing, T}
    capacitance::Union{Nothing, T}
    inductance::Union{Nothing, T}
end

Base.eltype(::NominalData{T}) where {T} = T
Base.eltype(::Type{NominalData{T}}) where {T} = T

function NominalData(;
        designation_code::Union{Nothing, AbstractString} = nothing,
        U0::Union{Nothing, Real} = nothing,
        U::Union{Nothing, Real} = nothing,
        conductor_cross_section::Union{Nothing, Real} = nothing,
        screen_cross_section::Union{Nothing, Real} = nothing,
        armor_cross_section::Union{Nothing, Real} = nothing,
        resistance::Union{Nothing, Real} = nothing,
        capacitance::Union{Nothing, Real} = nothing,
        inductance::Union{Nothing, Real} = nothing
)
    source = (
        U0, U, conductor_cross_section, screen_cross_section,
        armor_cross_section, resistance, capacitance, inductance
    )
    present = filter(!isnothing, source)
    T = isempty(present) ? Float64 : promote_type(typeof.(float.(present))...)
    values = map(value -> value === nothing ? nothing : convert(T, float(value)), source)
    return NominalData{T}(
        designation_code === nothing ? nothing : String(designation_code),
        values...
    )
end

function Base.convert(::Type{NominalData{T}}, data::NominalData) where {T <: Real}
    values = map(
        value -> value === nothing ? nothing : convert(T, value),
        (
            data.U0, data.U, data.conductor_cross_section,
            data.screen_cross_section, data.armor_cross_section,
            data.resistance, data.capacitance, data.inductance
        )
    )
    return NominalData{T}(data.designation_code, values...)
end

Base.convert(::Type{NominalData{T}}, data::NominalData{T}) where {T <: Real} = data
