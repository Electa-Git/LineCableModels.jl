function routes(identifier::Val{:Deri1981})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

assumptions(::Val{:Deri1981}) = (
    air = _lossless, earth = _conductive, permeability = vacuum_permeability
)
propagation(::Val{:Deri1981}) = Val(:zero)
media(::Formula{:Deri1981}) = Val(:stratified)

"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Ideal infinitely long parallel thin wires; conductor internal impedance and insulation are outside the formulation. The plane-field construction uses current per unit width to derive ``p``. |
| Calculated quantities | Effective complex depth for horizontally layered earth, inserted into the source's overhead-wire self and mutual impedance expressions |
| Earth structure | ``n`` finite horizontal layers of thickness ``\\Delta_K`` over an ``(n+1)``th layer extending to infinite depth. |
| Model and approximation | The one-dimensional layer fields and their flux are used to construct an effective homogeneous return plane. Retaining that plane as an image surface for thin conductors is explicitly heuristic; the paper says equal complex depth is expected to produce approximately equal impedances. This is not an exact multilayer Sommerfeld solution. |
| Main source | A. Déri, G. Tevan, A. Semlyen, and A. Castanheira (1981) |
| Citation key(s) | `:Deri1981` |
| Evidence status | PDF page images checked; a printed coefficient identity is internally inconsistent and is retained verbatim below rather than repaired |

**Expression.** Integrating the source field equation
``dE/dx=-j\\omega\\mu_0H`` from the surface to infinity gives
``p=E_0/(j\\omega\\mu_0H_0)``. Starting in the bottom half-space,
evaluate the field ratio upward through each finite layer:

```math
\\begin{aligned}
p_{\\mathrm{eff},N}&=p_N, \\\\
p_{\\mathrm{eff},k}&=
p_k\\frac{p_{\\mathrm{eff},k+1}+p_k\\tanh(\\Delta_k/p_k)}
{p_k+p_{\\mathrm{eff},k+1}\\tanh(\\Delta_k/p_k)}, \\\\
p_k&=(j\\omega\\mu_0\\sigma_k)^{-1/2}.
\\end{aligned}
```

This follows equations (9), (15), and (22)–(24), with the decaying bottom
field. It avoids products of growing hyperbolic functions. The resulting
depth enters the self and mutual image-distance logarithms of Dubanton1969;
it does not replace an exact layered-earth spectral integral.

**Reference.** [Deri1981](@cite), equations (9), (15), (22)–(24), and the
equivalent surface field ratio in equation (51).
"""
description(::Formula{:Deri1981}) =
    "Deri et al. multilayer complex ground-return plane (1981)"

function propagation_constant(::Val{:Deri1981}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Deri1981})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    all(value -> isfinite(value) && value > zero(value), rho[2:end]) ||
        throw(DomainError(rho, ":Deri1981 requires positive finite earth resistivities"))
    iszero(jω) && throw(DomainError(jω, ":Deri1981 requires nonzero frequency"))
    functor = _stratified_functor(
        Val(:Deri1981), formula, rho, epsilon, mu, jω, Γ, segments, thickness
    )
    all(value -> isfinite(value) && value >= zero(value), thickness[2:(end-1)]) ||
        throw(DomainError(thickness, "finite earth-layer thicknesses must be nonnegative"))
    p = inv(functor.state.gamma[end])
    for layer in (length(rho)-1):-1:2
        local_depth = inv(functor.state.gamma[layer])
        t = tanh(thickness[layer] / local_depth)
        p = local_depth * (p + local_depth * t) / (local_depth + p * t)
    end
    state = merge(functor.state, (return_plane_depth = p,))
    return Functor{:Deri1981, typeof(formula.routes), typeof(state)}(formula.routes, state)
end

function (formula::Formula{:Deri1981})(rho, epsilon, mu, jω, Γ, segments = nothing)
    length(rho) == 2 || throw(DimensionMismatch(
        ":Deri1981 requires thicknesses for more than one earth layer"
    ))
    return formula(rho, epsilon, mu, jω, Γ, segments, fill(oftype(first(rho), Inf), 2))
end

function earth_impedance(::Val{:Deri1981}, ::Val{:self}, functor, pair)
    return _image_plane_impedance(functor, pair, functor.state.return_plane_depth, Val(:self))
end

function earth_impedance(::Val{:Deri1981}, ::Val{:mutual}, functor, pair)
    return _image_plane_impedance(functor, pair, functor.state.return_plane_depth, Val(:mutual))
end

:Deri1981
