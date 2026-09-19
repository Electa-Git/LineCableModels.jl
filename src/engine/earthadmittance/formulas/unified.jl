function assumptions(::Val{:unified})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Circumferentially averaged, two-half-space earth return
with the complete enclosed-current constraint. Rows are receivers and columns
are sources. Air, earth and both ordered mixed directions are supported, with
independent permeabilities and caller-prescribed Γ (zero by default).

**Expression.** Assemble the source kernels and solve the complete system:

```math
P_e L=H,\\qquad Z_e L=K+\\Gamma^2 H/s,\\qquad Y_e H=sL,
\\quad L=A_r^{-1}-F_rK,\\quad K=\\mathcal Z-\\Gamma^2\\mathcal P_\\phi/s.
```

Here s=jω, Ze has units \\[Ω/m\\], Pe has units \\[m/F\\], and Ye has units
\\[S/m\\]. The shared implementation forms physical matrices before selecting
entries. Air receivers use the interface voltage; earth receivers use deep-earth
voltage. These references are fixed by receiving layer, not selectable parameters.

**Reference.** User-supplied manuscript, *Unified circumferentially averaged
framework for overhead, buried, and mixed conductor systems*. Current closure, reference conventions and the
equal-medium limit are exercised in `test/unit/engine/unified_earth_return.jl`.
These scoped controls do not establish acceptance of arbitrary complete matrices.
"""
function description(::Type{<:Formula{:unified}}; compact::Bool = false)
    compact ? "Unified" :
    "Unified circumferential earth potential with full current closure"
end

function earth_bindings(selected::Formula{:unified},
        physical::AbstractVector{<:EarthPair}, homogeneous, indices)
    binding = invoke(earth_bindings,
        Tuple{EarthAdmittanceFormulation, AbstractVector{<:EarthPair}, Any, Any},
        selected, physical, homogeneous, collect(eachindex(physical)))
    EarthImpedance._unified_geometry(physical)
    options = first(binding.equations).declaration.options
    all(g -> isequal(g.declaration.options, options), binding.equations) ||
        throw(ArgumentError("the unified current closure requires common integration controls"))
    equations = [(declaration = group.declaration,
                     indices = intersect(group.indices, indices))
                 for group in binding.equations if !isdisjoint(group.indices, indices)]
    return merge(binding, (; equations))
end

function computation_type(::Type{T},
        selected::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        frequencies) where {T <: Real}
    argument = selected.parameters.Γ
    argument isa AbstractVector && length(argument) != length(frequencies) &&
        throw(DimensionMismatch("unified Γ must contain one value per frequency sample"))
    argument isa Number && return promote_type(T, typeof(real(argument)), typeof(imag(argument)))
    return foldl(argument; init = T) do scalar, value
        promote_type(scalar, typeof(real(value)), typeof(imag(value)))
    end
end

function initialize_buffers(
        selected::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        ::Type{T}, input, invariants, buffers) where {T}
    buffers = initialize_buffers(selected.equivalent_earth, T, input, invariants, buffers)
    buffers = initialize_buffers(Val(:quad), T, input, invariants, buffers)
    R=typeof(float(nominal(one(T))))
    haskey(buffers, :unified) && return buffers
    n = length(invariants.geometry.radius)
    arrays = ntuple(_ -> Matrix{Complex{T}}(undef, n, n), 7)
    current = (x = Vector{Complex{T}}(undef, n), scaling = Vector{T}(undef, n),
        A = Vector{Complex{T}}(undef, n), F = Vector{Complex{T}}(undef, n))
    return merge(buffers,
        (unified = (
            K = arrays[1], H = arrays[2], L = arrays[3], Ze = arrays[4], Pe = arrays[5],
            factor = arrays[6], rhs = arrays[7], current = current,
            points = sizehint!(R[], 128), seeds = sizehint!(R[], 128),
            scales = sizehint!(R[], 16)),))
end

function earth!(Z, P, ::Nothing, selected::Formula{:unified},
        ::Nothing, binding, materials, workspace, frequency)
    EarthImpedance._unified_current!(workspace, materials, binding, frequency)
    earth!(P, selected, binding, materials, workspace.input.jω[frequency],
        workspace, materials.thickness)
    return workspace
end

function earth!(destination, ::Formula{:unified},
        group::NamedTuple{(:declaration, :indices)}, interactions,
        materials, jω, workspace, thickness)
    for index in group.indices
        pair=interactions[index].pair
        value=group.declaration.equation(nothing, pair, workspace)
        isfinite(value) ||
            throw(DomainError(value, "earth potential coefficient must be finite"))
        destination[pair.row, pair.column]=oftype(jω, value)
    end
    return destination
end

function earth_potential_coefficient(
        ::Formula{:unified}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Pe[pair.row, pair.column]
end

function earth_potential_coefficient(
        ::Formula{:unified}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Pe[pair.row, pair.column]
end

function earth_potential_coefficient(::Formula{:unified}, ::Val{:mutual},
        ::Val{1}, ::Val{2}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Pe[pair.row, pair.column]
end

function earth_potential_coefficient(::Formula{:unified}, ::Val{:mutual},
        ::Val{2}, ::Val{1}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Pe[pair.row, pair.column]
end

function formulation_options(::FormulaMethod{<:Formula{:unified},
        typeof(earth_potential_coefficient),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return FormulationOptions((integration = (method = :quad, options = (;)),))
end

function validate(
        binding::FormulaMethod{<:Formula{:unified}, typeof(earth_potential_coefficient)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    return binding
end

function earth_bindings(z::EarthImpedance.Formula{:unified}, p::Formula{:unified}, impedance, admittance)
    same_physical_state(z.parameters, p.parameters) &&
    same_physical_state(z.equivalent_earth, p.equivalent_earth) &&
    isequal(first(impedance.equations).declaration.options,
        first(admittance.equations).declaration.options) ||
        return nothing
    return (; impedance, admittance)
end

function earth!(Z, P, z::EarthImpedance.Formula{:unified}, p::Formula{:unified},
        impedance, admittance, materials, workspace, frequency)
    EarthImpedance._unified_current!(workspace, materials, impedance, frequency)
    earth!(Z, z, impedance, materials, workspace.input.jω[frequency],
        workspace, materials.thickness)
    earth!(P, p, admittance, materials, workspace.input.jω[frequency],
        workspace, materials.thickness)
    return workspace
end

:unified
