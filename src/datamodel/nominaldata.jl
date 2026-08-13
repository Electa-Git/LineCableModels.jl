
"""
$(TYPEDEF)

Stores nominal electrical and geometric parameters for a cable design.

$(TYPEDFIELDS)
"""
struct NominalData{T<:REALSCALAR}
    "Cable designation as per DIN VDE 0271/0276."
    designation_code::Union{Nothing,String}
    "Rated phase-to-earth voltage \\[kV\\]."
    U0::Union{Nothing,T}
    "Rated phase-to-phase voltage \\[kV\\]."
    U::Union{Nothing,T}
    "Cross-sectional area of the conductor \\[mm²\\]."
    conductor_cross_section::Union{Nothing,T}
    "Cross-sectional area of the screen \\[mm²\\]."
    screen_cross_section::Union{Nothing,T}
    "Cross-sectional area of the armor \\[mm²\\]."
    armor_cross_section::Union{Nothing,T}
    "Base (DC) resistance of the cable core \\[Ω/km\\]."
    resistance::Union{Nothing,T}
    "Capacitance of the main insulation \\[μF/km\\]."
    capacitance::Union{Nothing,T}
    "Inductance of the cable (trifoil formation) \\[mH/km\\]."
    inductance::Union{Nothing,T}

    # --- Tight / typed kernel: assumes values already coerced to T (or nothing)
    @inline function NominalData{T}(;
        designation_code::Union{Nothing,String}=nothing,
        U0::Union{Nothing,T}=nothing,
        U::Union{Nothing,T}=nothing,
        conductor_cross_section::Union{Nothing,T}=nothing,
        screen_cross_section::Union{Nothing,T}=nothing,
        armor_cross_section::Union{Nothing,T}=nothing,
        resistance::Union{Nothing,T}=nothing,
        capacitance::Union{Nothing,T}=nothing,
        inductance::Union{Nothing,T}=nothing,
    ) where {T<:REALSCALAR}
        new{T}(
            designation_code,
            U0,
            U,
            conductor_cross_section,
            screen_cross_section,
            armor_cross_section,
            resistance,
            capacitance,
            inductance,
        )
    end
end

"""
$(TYPEDSIGNATURES)

Weakly-typed constructor that infers the target scalar type `T` from the **provided numeric kwargs** (ignoring `nothing` and the string designation), coerces numerics to `T`, and calls the strict kernel.

If no numeric kwargs are provided, it defaults to `Float64`.
"""
@inline function NominalData(;
    designation_code::Union{Nothing,String}=nothing,
    U0::Union{Nothing,Number}=nothing,
    U::Union{Nothing,Number}=nothing,
    conductor_cross_section::Union{Nothing,Number}=nothing,
    screen_cross_section::Union{Nothing,Number}=nothing,
    armor_cross_section::Union{Nothing,Number}=nothing,
    resistance::Union{Nothing,Number}=nothing,
    capacitance::Union{Nothing,Number}=nothing,
    inductance::Union{Nothing,Number}=nothing,
)
    # collect provided numerics (skip `nothing`)
    nums = Tuple(x for x in
                 (U0, U, conductor_cross_section, screen_cross_section, armor_cross_section,
        resistance, capacitance, inductance) if x !== nothing)

    # infer T from numerics, fallback to Float64 if none
    T = isempty(nums) ? Float64 : resolve_T(nums...)

    return NominalData{T}(;
        designation_code=designation_code,
        U0=(U0 === nothing ? nothing : coerce_to_T(U0, T)),
        U=(U === nothing ? nothing : coerce_to_T(U, T)),
        conductor_cross_section=(conductor_cross_section === nothing ? nothing : coerce_to_T(conductor_cross_section, T)),
        screen_cross_section=(screen_cross_section === nothing ? nothing : coerce_to_T(screen_cross_section, T)),
        armor_cross_section=(armor_cross_section === nothing ? nothing : coerce_to_T(armor_cross_section, T)),
        resistance=(resistance === nothing ? nothing : coerce_to_T(resistance, T)),
        capacitance=(capacitance === nothing ? nothing : coerce_to_T(capacitance, T)),
        inductance=(inductance === nothing ? nothing : coerce_to_T(inductance, T)),
    )
end

include("nominaldata/base.jl")


