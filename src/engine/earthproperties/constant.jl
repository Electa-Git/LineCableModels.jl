defaults(::Type{Formula{:constant}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the properties of the earth `layer` at the frequency `f` by the formula
`Formula{:constant}`, whose equation is in the section "Formulas" of
[`EarthProperties`](@ref).

# Arguments

- `Formula{:constant}`: frequency-independent properties law.
- `layer`: layer of the earth or of the air.
- `f`: frequency \\[Hz\\].

# Returns

- `rho`: resistivity \\[Ω·m\\].
- `eps_r`: relative permittivity \\[dimensionless\\].
- `mu_r`: relative permeability \\[dimensionless\\].
"""
function properties(::Formula{:constant}, layer::EarthLayer{T}, f::T) where {T <: Real}
    return (rho = layer.rho, eps_r = layer.eps_r, mu_r = layer.mu_r)
end
