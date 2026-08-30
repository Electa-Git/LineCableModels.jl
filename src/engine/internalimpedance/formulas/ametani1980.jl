function routes(::Val{:Ametani1980})
    (
        inner = schelkunoff_inner,
        outer = schelkunoff_outer,
        mutual = schelkunoff_mutual
    )
end

assumptions(::Val{:Ametani1980}) = (;)

function description(::Formula{:Ametani1980})
    "Ametani single-core coaxial cable impedance assembly (1980)"
end

raw"""
Construct the cylindrical surface terms used by Ametani's single-core coaxial
cable impedance assembly:

```math
\mathbf Z_{ij}^{(c,s,a)}=
\begin{bmatrix}Z_{ccj}&Z_{csj}&Z_{caj}\\
Z_{csj}&Z_{ssj}&Z_{saj}\\Z_{caj}&Z_{saj}&Z_{aaj}\end{bmatrix},
```

```math
Z_{cc}=z_{cs}+z_{sa}+z_{a4}-2z_{2m}-2z_{3m},\quad
Z_{ss}=z_{sa}+z_{a4}-2z_{3m},\quad Z_{aa}=z_{a4},
```

```math
Z_{cs}=z_{sa}+z_{a4}-z_{2m}-2z_{3m},\qquad
Z_{ca}=Z_{sa}=z_{a4}-z_{3m}.
```

The engine's outward conductor loop performs these additions and mutual-term
subtractions for any number of concentric conductors. Each ``z_{ki}``,
``z_{ko}``, and ``z_{km}`` is evaluated with the exact Schelkunoff cylindrical
surface expressions; intervening ``z_{12}``, ``z_{23}``, and ``z_{34}`` are
provided by the selected insulation-impedance route.
"""
function ametani1980(
        formula::Formula{:Ametani1980},
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    mu_c = vacuum_permeability(r_in) * mur_c
    sigma_c = conductivity(rho_c)
    m = sqrt(jω * mu_c * sigma_c)
    w_ex = m * r_ex
    state = (; r_in, r_ex, rho_c, mur_c, jω, mu_c, sigma_c, m, w_ex)
    return Functor{:Ametani1980, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

@inline function (formula::Formula{:Ametani1980})(r_in, r_ex, rho_c, mur_c, jω)
    return ametani1980(formula, r_in, r_ex, rho_c, mur_c, jω)
end

@inline function (functor::Functor{:Ametani1980})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:Ametani1980})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:Ametani1980})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

:Ametani1980
