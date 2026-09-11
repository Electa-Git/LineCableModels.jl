"""
$(TYPEDEF)

Frequency-independent description of dielectric materials combined radially
in series. Constituent kinds retain independent insulation and semicon
formulation selection. No constitutive law or reference frequency is selected.

For logarithmic radial weights ``q_i = \\log(r_{i+1}/r_i)``, the effective
admittivity is ``\\kappa_{eq} = \\sum_i q_i / \\sum_i(q_i/\\kappa_i)``.
The engine evaluates each ``\\kappa_i`` before combining it.

The displayed `rho` is the DC radial resistivity and `eps_r` is the lossless
radial permittivity, not a lossy fit. There is no frequency-independent
equivalent loss tangent in general.

$(TYPEDFIELDS)
"""
struct RadialDielectric{T <: Real} <: AbstractMaterial
    "Original insulation and semicon materials."
    materials::Vector{Material{T}}
    "Positive logarithmic radial weights, dimensionless."
    weights::Vector{T}
    "Passive region classification; constituent kinds remain authoritative."
    kind::Symbol
    "Series DC resistivity, in Ω·m."
    rho::T
    "Lossless relative permittivity, dimensionless."
    eps_r::T
    "Equivalent relative permeability, dimensionless."
    mu_r::T
    "Shared constituent reference temperature, in °C."
    T0::T
    function RadialDielectric(materials::Vector{Material{T}}, weights::Vector{T},
            mu_r::T) where {T <: Real}
        isempty(materials) && throw(ArgumentError("RadialDielectric requires at least one constituent"))
        length(materials) == length(weights) || throw(DimensionMismatch(
            "RadialDielectric requires one radial weight per constituent"))
        total = sum(weights)
        rho = sum(w * m.rho for (w, m) in zip(weights, materials)) / total
        eps_r = total / sum(w / m.eps_r for (w, m) in zip(weights, materials))
        return validate(new{T}(copy(materials), copy(weights), :insulator,
            rho, eps_r, mu_r, first(materials).T0))
    end
end

"""
$(TYPEDSIGNATURES)

Describe a radial series composition without selecting its dielectric law.

# Arguments

- `materials`: Physical [`Material`](@ref) constituents, each classified as
  `:insulator` or `:semicon`.
- `weights`: Positive logarithmic radius ratios, dimensionless.

# Keywords

- `mu_r`: Equivalent relative permeability. Defaults to the radial weighted
  mean; homogenization may supply its helical-solenoid correction.

# Returns

- A [`RadialDielectric`](@ref) with promoted concrete scalar storage.
"""
function RadialDielectric(materials, weights;
        mu_r = sum(w * m.mu_r for (w, m) in zip(weights, materials)) / sum(weights))
    T = promote_type(typeof(float(mu_r)), map(eltype, materials)...,
        map(x -> typeof(float(x)), weights)...)
    return RadialDielectric(Material{T}[convert(Material{T}, m) for m in materials],
        T[weights...], convert(T, mu_r))
end

function validate(material::RadialDielectric)
    all(w -> isfinite(w) && w > zero(w), material.weights) || throw(ArgumentError(
        "RadialDielectric.weights must be finite and positive logarithmic radius ratios"))
    for (index, constituent) in enumerate(material.materials)
        validate(constituent)
        constituent.kind in (:insulator, :semicon) || throw(ArgumentError(
            "RadialDielectric.materials[$index] must be insulation or semicon, not :$(constituent.kind)"))
        isapprox(constituent.T0, material.T0) || throw(ArgumentError(
            "RadialDielectric constituents must share one reference temperature"))
    end
    isfinite(material.mu_r) && material.mu_r > zero(material.mu_r) ||
        throw(ArgumentError("RadialDielectric.mu_r must be finite and positive"))
    return material
end

Base.eltype(::RadialDielectric{T}) where {T} = T
Base.eltype(::Type{RadialDielectric{T}}) where {T} = T
function Base.convert(::Type{RadialDielectric{T}}, material::RadialDielectric) where {T <: Real}
    return RadialDielectric(Material{T}[convert(Material{T}, m) for m in material.materials],
        T[material.weights...], convert(T, material.mu_r))
end
Base.convert(::Type{RadialDielectric{T}}, material::RadialDielectric{T}) where {T <: Real} = material
Base.:(==)(a::RadialDielectric, b::RadialDielectric) =
    a.materials == b.materials && a.weights == b.weights && a.mu_r == b.mu_r
Base.isequal(a::RadialDielectric, b::RadialDielectric) =
    isequal(a.materials, b.materials) && isequal(a.weights, b.weights) && isequal(a.mu_r, b.mu_r)
Base.hash(material::RadialDielectric, h::UInt) =
    hash(material.mu_r, hash(material.weights, hash(material.materials, hash(:RadialDielectric, h))))
TextDisplay.name(::Type{<:RadialDielectric}) = "RadialDielectric"
Base.summary(io::IO, material::RadialDielectric) =
    print(io, "RadialDielectric · ", length(material.materials), " constituents")
Base.show(io::IO, material::RadialDielectric) = summary(io, material)
function Base.show(io::IO, ::MIME"text/plain", material::RadialDielectric)
    get(io, :compact, false) && return show(io, material)
    return TextDisplay.fields(io, "RadialDielectric", (
        constituents = length(material.materials),
        ρ_DC = TextDisplay.engineering(material.rho, :ohm_meter),
        εᵣ_lossless = TextDisplay.value(material.eps_r),
        μᵣ = TextDisplay.value(material.mu_r)); multiline = true)
end
