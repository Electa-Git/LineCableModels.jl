function routes(identifier::Val{:Wise1948})
    return (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Wise1948})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Wise1948}) = Val(:zero)
function Formula(::Val{:Wise1948};approximation::Symbol=:integral,kwargs...)
    approximation in (:integral,:rational,:coarse,:small_g) ||
        throw(ArgumentError(":Wise1948 approximation must be :integral, :rational, :coarse, or :small_g"))
    identifier=Val(:Wise1948)
    defaults=routes(identifier)
    if approximation!==:integral
        defaults=merge(defaults,(
            self=FormulaMethod(identifier,earth_potential_coefficient,Val(:approximation)),
            mutual=FormulaMethod(identifier,earth_potential_coefficient,Val(:approximation))))
    end
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) ||
        throw(ArgumentError("unknown routes for earth-admittance formula :Wise1948"))
    selected=merge(defaults,overrides)
    values=merge(assumptions(identifier),(;approximation))
    return Formula{:Wise1948,typeof(selected),typeof(values)}(selected,values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinitely long, straight, thin parallel wire(s); self example uses wire radius ``a`` and height ``h``; mutual geometry uses direct distance ``\\rho'`` and image distance ``\\rho''``. No insulation layer is included. |
| Calculated quantities | Self and mutual scalar potential coefficients for one or two parallel overhead wires; source-provided capacitance/admittance assembly statements |
| Earth structure | Homogeneous flat half-space. |
| Model and approximation | Equations (7)–(8) are the source's integral representation after assigning ``\\gamma=jk`` as a first approximation to longitudinal propagation. They are not full-wave expressions for an independently solved ``\\gamma``. Wise's separate analytical approximation to this integral is recorded in [the companion approximation record]. |
| Main source | W. Howard Wise (1948) |
| Citation key(s) | `:Wise1948` |
| Evidence status | Original PDF page images checked |

**Expression.**

```math
P_{e,ij}=\\frac{P_{0,ij}+C_{ij}}{2\\pi\\varepsilon_0},
```

```math
C_{ij}=2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{(\\gamma_1^2/\\gamma_0^2)\\lambda+
\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}}d\\lambda,
\\quad P_{0,ij}=\\ln(D_{ij}/d_{ij}).
```

**Reference.** W. H. Wise, “Potential Coefficients for Ground Return
Circuits,” *Bell System Technical Journal*, 27, 365–371, 1948.

## Analytical approximations

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite straight thin parallel wire(s); the approximation variables use ``h+z`` and ``y``. |
| Calculated quantities | Analytical approximations to the complex earth correction ``M+jN`` used in overhead self and mutual potential coefficients |
| Earth structure | Homogeneous flat half-space. |
| Model and approximation | The parent integral is [Wise's equation (8)]. The source expands ``\\sqrt{t^2+u^2}=t+u-tu/(t+u)+\\cdots``, retains the displayed truncated expression, performs partial fractions with roots ``r_1,r_2``, and obtains (9). Equation (10) imposes ``y=0``. The final unnumbered form further assumes ``g`` is very small. The source gives no formal remainder bound; its reported sub-one-percent discrepancy applies only to checked points. |
| Main source | W. Howard Wise (1948) |
| Citation key(s) | Primary: `:Wise1948`; later MKS applications: `:Nakagawa1981a`, `:Nakagawa1981b` |
| Evidence status | Original PDF page images checked; later applications retain the same correction |

**Expression.** The default `approximation=:integral` retains (8).
Select `approximation=:rational` for (9), including its `y=0`
specialization (10). In SI variables, set
``a=\\gamma_g^2/\\gamma_a^2``,
``\\beta=\\sqrt{\\gamma_g^2-\\gamma_a^2}``,
``q_{1,2}=\\beta[1\\mp\\sqrt{1-4/(a+1)}]/2``, and
``F(q,z)=e^{qz}E_1(qz)``. The correction added to the ideal-ground
logarithm is

```math
C=\\sum_{z\\in\\{H-jx,H+jx\\}}
\\frac{q_2F(q_1,z)-q_1F(q_2,z)}
{(a+1)(q_2-q_1)}.
```

Thus ``P=[\\ln(D/d)+C]/(2\\pi\\varepsilon_a)``, with
``C=2(M+jN)`` in the source's notation. The scaled exponential
integral is continued from its Laplace integral as the separation angle
changes; principal values of the product ``qz`` alone would introduce
a spurious branch jump for widely separated wires. The smaller root is
obtained from the root product to avoid cancellation.

For zero horizontal separation only, `approximation=:coarse` selects
the first unnumbered approximation on printed p. 371,
``C=2F(\\beta/a,H)/a``. The still coarser
`approximation=:small_g` selects
``C=2\\ln(1+a/u)/a+2E_1(g)/(1+a)``, where
``u=\\beta/|\\beta|`` and ``g=H|\\beta|\\ll1``.
The latter rejects ``g\\ge1``; this necessary restriction does not
establish an accuracy bound within the accepted interval.

**Geometry and limits.** All three analytical selections require
nonmagnetic media, lossless air, and positive earth conductivity.
For a self coefficient ``x=0``, ``D=2h``, and ``d=r``.
The conductor radius does not become a horizontal source displacement.
The two coarse selections also apply to vertically aligned mutual pairs;
they reject nonzero horizontal separation. The phase-domain potential
matrix is assembled before inversion. The source's reported accuracy
is not a global bound.

**Reference.** [Wise1948](@cite), equations (7)–(10) and the two
unnumbered zero-separation approximations on printed p. 371.
"""
function description(::Formula{:Wise1948})
    "Wise wideband homogeneous-earth overhead potential coefficient (1948)"
end

function propagation_constant(
        ::Val{:Wise1948}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Wise1948})(rho, epsilon, mu, jω, Γ, segments = nothing)
    if hasproperty(formula.assumptions,:approximation) && formula.assumptions.approximation!==:integral
        length(rho)==length(epsilon)==length(mu)==2 ||
            throw(ArgumentError("Wise1948 analytical expressions require two homogeneous half-spaces"))
        all(x->isapprox(x,vacuum_permeability(x)),mu) ||
            throw(DomainError(mu,"Wise1948 assumes nonmagnetic media"))
        all(x->isfinite(x)&&x>0,epsilon) ||
            throw(DomainError(epsilon,"Wise1948 requires positive finite permittivity"))
        isfinite(jω) && !iszero(jω) && iszero(real(jω)) ||
            throw(DomainError(jω,"Wise1948 requires nonzero real frequency"))
    end
    return _homogeneous_functor(Val(:Wise1948), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Wise's wideband overhead earth potential coefficient:

```math
P_{e,ij}=\frac{P_{0,ij}+C_{ij}}{2\pi\varepsilon_0},
```

```math
C_{ij}=2\int_0^\infty
\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{(\gamma_1^2/\gamma_0^2)\lambda+
\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}}d\lambda,
\qquad P_{0,ij}=\ln(D_{ij}/d_{ij}).
```
"""
function earth_potential_coefficient(
        ::Val{:Wise1948}, ::Val{:mutual}, functor, pair
)
    state = functor.state
    geometry = earth_potential_coefficient(Val(:Wise1948),Val(:geometry),pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_1_squared / gamma_0_squared
    integral = _quadrature(state) do t
        lambda=t/geometry.H
        radial = sqrt(lambda^2 + gamma_1_squared - gamma_0_squared)
        exp(-t) * cos(geometry.x * lambda) /
        ((ratio * lambda + radial)*geometry.H)
    end
    return _complex_result(state.jω,(log(geometry.D / geometry.d) + 2 * integral) /
           (2*(one(geometry.H)*π) * state.epsilon[1]))
end

function earth_potential_coefficient(::Val{:Wise1948},::Val{:geometry},pair)
    _require(pair,Val(:overhead))
    H=sum(pair.heights)
    self=pair.row==pair.column
    x=self ? zero(pair.separation) : abs(pair.separation)
    d=self ? pair.separation : hypot(x,pair.heights[1]-pair.heights[2])
    all(isfinite,(H,x,d)) && minimum(pair.heights)>0 && d>0 ||
        throw(DomainError((pair.heights,pair.separation),":Wise1948 requires positive overhead geometry"))
    return (;H,x,d,D=hypot(H,x))
end

function earth_potential_coefficient(
        ::Val{:Wise1948},::Val{:laplace},q::Complex{T},z::Complex{T}
) where {T <: AbstractFloat}
    product=q*z
    if real(product)<0 && abs(imag(product))<=8eps(T)*abs(product)
        H=real(z)
        return quadgk(t->exp(-z*t/H)/(t/H+q)/H,zero(T),T(Inf);
            rtol=max(sqrt(eps(T)),T(1e-30)))[1]
    end
    value=scaled_expint_negative(-product)
    angle=atan(imag(q),real(q))+atan(imag(z),real(z))
    if angle>T(π)
        value-=complex(zero(T),2T(π))*exp(product)
    elseif angle < -T(π)
        value+=complex(zero(T),2T(π))*exp(product)
    end
    return value
end

function earth_potential_coefficient(
        ::Val{:Wise1948},::Val{:analytical},a::Complex{T},beta::Complex{T},
        H::T,x::T,method::Val
) where {T <: AbstractFloat}
    if T <: Union{Float32,Float64} && precision(BigFloat)>=128
        return Complex{T}(earth_potential_coefficient(Val(:Wise1948),Val(:analytical),
            Complex{BigFloat}(a),Complex{BigFloat}(beta),BigFloat(H),BigFloat(x),method))
    end
    if method===Val(:rational)
        q2=beta*(1+sqrt(1-4/(a+1)))/2
        q1=beta^2/((a+1)*q2)
        if abs(q2-q1)<=sqrt(eps(T))*max(abs(q1),abs(q2))
            # Evaluate the same rational function at coalescing roots.
            integrand=t->begin
                lambda=t/H
                2exp(-t)*cos(x*lambda)*(lambda+beta)/
                    ((a+1)*lambda^2+(a+1)*beta*lambda+beta^2)/H
            end
            return quadgk(integrand,zero(T),T(Inf);rtol=max(sqrt(eps(T)),T(1e-30)))[1]
        end
        total=zero(a)
        for z in (complex(H,-x),complex(H,x))
            first=earth_potential_coefficient(Val(:Wise1948),Val(:laplace),q1,z)
            second=earth_potential_coefficient(Val(:Wise1948),Val(:laplace),q2,z)
            total+=(q2*first-q1*second)/((a+1)*(q2-q1))
        end
        return total
    end
    iszero(x) || throw(ArgumentError("Wise1948 coarse approximations require zero horizontal separation"))
    if method===Val(:coarse)
        return 2earth_potential_coefficient(Val(:Wise1948),Val(:laplace),
            beta/a,complex(H,zero(T)))/a
    end
    method===Val(:small_g) || throw(ArgumentError("unknown Wise1948 analytical selection"))
    scale=abs(beta); g=H*scale
    g<1 || throw(DomainError(g,"Wise1948 small-g expression requires H*abs(beta) much smaller than one"))
    u=beta/scale
    E1=exp(-g)*scaled_expint_negative(complex(-g,zero(T)))
    return 2log1p(a/u)/a+2*E1/(1+a)
end

function earth_potential_coefficient(::Val{:Wise1948},::Val{:approximation},functor,pair)
    state=functor.state
    isinf(state.rho[1]) && isfinite(state.rho[2]) && state.rho[2]>0 ||
        throw(DomainError(state.rho,"Wise1948 analytical expressions require lossless air and conducting earth"))
    geometry=earth_potential_coefficient(Val(:Wise1948),Val(:geometry),pair)
    g0,g1=state.gamma_medium_squared
    beta=sqrt(g1-g0)
    correction=earth_potential_coefficient(Val(:Wise1948),Val(:analytical),
        g1/g0,beta,geometry.H,geometry.x,Val(state.formula.assumptions.approximation))
    return _complex_result(state.jω,(log(geometry.D/geometry.d)+correction)/
        (2*(one(geometry.H)*π)*state.epsilon[1]))
end

:Wise1948
