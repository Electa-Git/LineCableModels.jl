function description(::Type{<:Formula{:pollaczek1926}}; compact::Bool = false)
    compact ? "Pollaczek" : "Pollaczek homogeneous-earth underground impedance (1926)"
end

"""
$(TYPEDSIGNATURES)

Evaluate the classical Pollaczek underground earth-return impedance \\[Ω/m\\]
for parallel filaments in homogeneous conductive earth below non-conducting air:

```math
Z_{e,ij}=\\frac{s\\mu_0}{2\\pi}\\left[
K_0(md)-K_0(mD)+2\\int_0^\\infty
\\frac{e^{-\\ell\\sqrt{\\lambda^2+m^2}}}
{\\lambda+\\sqrt{\\lambda^2+m^2}}\\cos(x\\lambda)\\,d\\lambda\\right],
\\qquad m=\\sqrt{\\frac{s\\mu_0\\mu_{r,g}}{\\rho_g}}.
```

Here ``s=\\sigma+j\\omega`` \\[1/s\\] has nonnegative real part and is nonzero.
The principal square roots give the analytic continuation to complex frequency.
Earth resistivity ``\\rho_g`` \\[Ω·m\\] and absolute permeability
``\\mu_0\\mu_{r,g}`` \\[H/m\\] enter ``m`` \\[1/m\\]. Earth permittivity does
not enter because Pollaczek neglects displacement current in earth. The prefactor
uses vacuum permeability ``\\mu_0`` \\[H/m\\]. The expression sets the longitudinal propagation constant to zero.

Positive burial depths ``h_i,h_j`` \\[m\\] give ``\\ell=h_i+h_j``,
``d=\\sqrt{x^2+(h_i-h_j)^2}``, and ``D=\\sqrt{x^2+\\ell^2}`` \\[m\\].
For a mutual term, ``x`` is horizontal separation \\[m\\]. For a self term,
``x=d=r_i`` is the exterior cable radius \\[m\\] and ``h_i=h_j``.
``K_0`` is the modified Bessel function of the second kind, order zero, and
``\\lambda`` is the Fourier integration variable \\[1/m\\].

The integral uses the shared adaptive Gauss–Kronrod method on `[0, Inf)`.
The formula's `integration=(method=:quad, options=(rtol=1e-8, atol=0.0,
maxevals=10^7))` controls its estimated error. Subdivisions resolve both
``|m|`` and ``1/\\ell`` without a finite spectral cutoff.

# Arguments

- `formula`: Pollaczek formulation and its integration controls.
- `kind`: self or mutual interaction.
- `source`, `target`: buried source and target layer tags, both `Val(2)`.
- `functor`: evaluated material quantities, complex frequency and pair geometry.
- `workspace`: computation buffers, or `nothing` for a standalone evaluation.

# Returns

- Earth-return impedance \\[Ω/m\\]. Separate formulas compute internal cable contributions.

# Reference

F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926.

The underground parent expressions appear in Saad et al., Eqs. (1)-(4):

O. Saad, G. Gaba, and M. Giroux, “A Closed-Form Approximation for Ground Return
Impedance of Underground Cables,” 1996.
*IEEE Transactions on Power Delivery*, 11(3), 1536–1545.

Their full-line ``J_m`` and ``J_s``
integrands have an even kernel in ``\\gamma`` through ``|\\gamma|`` and
``\\gamma^2``. The odd sine part vanishes.
The result is twice the half-line cosine integral evaluated here.
"""
function earth_impedance(
        formula::Formula{:pollaczek1926}, kind::Union{Val{:self}, Val{:mutual}},
        source::Val{2}, target::Val{2}, functor, workspace)
    pair = functor.input.pair
    s = functor.input.jω
    μ0 = vacuum_permeability(typeof(real(s)))
    m2 = s * functor.input.mu[2] / functor.input.rho[2]
    m = sqrt(m2)
    hi, hj = abs.(pair.heights)
    x = something(pair.radius, pair.separation)
    ell = hi + hj
    d, D = hypot(x, hi - hj), hypot(x, ell)
    integral = SpectralIntegral(λ -> begin
        a = sqrt(λ^2 + m2)
        exp(-ell * a) * cos(x * λ) / (λ + a)
    end)
    integration = functor.input.options.data.integration
    buffers = workspace === nothing ? nothing : workspace.buffers
    R = typeof(float(nominal(real(s))))
    context = (formula = :pollaczek1926, frequency = imag(s) / (2π),
        receiver = pair.row, source = pair.column)
    value, _ = integrate(integral, integration.method, integration.options, buffers;
        points = (R(abs(nominal(m))), inv(R(nominal(ell)))), coordinate_type = R,
        context, observations = buffers === nothing ? nothing : buffers.observations)
    return s * μ0 / (2 * (one(real(s)) * π)) *
           (special_besselk(0, m * d) - special_besselk(0, m * D) + 2value)
end

function formulation_options(::Expression{<:Formula{:pollaczek1926}, typeof(earth_impedance)})
    return FormulationOptions((integration = (method = :quad, options = (;)),))
end

function initialize_buffers(selected::Formula{:pollaczek1926}, ::Type{T}, input, plan,
        buffers) where {T}
    buffers = initialize_buffers(selected.equivalent_earth, T, input, plan, buffers)
    return initialize_buffers(SpectralIntegral, Val(:quad), T, input, plan, buffers)
end

:pollaczek1926
