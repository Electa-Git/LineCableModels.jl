function routes(::Val{:Schelkunoff})
    (
        inner = schelkunoff_inner,
        outer = schelkunoff_outer,
        mutual = schelkunoff_mutual
    )
end

assumptions(::Val{:Schelkunoff}) = (;)

description(::Formula{:Schelkunoff}) = "Schelkunoff"

function (formula::Formula{:Schelkunoff})(
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
    state = (;
        r_in,
        r_ex,
        rho_c,
        mur_c,
        jω,
        mu_c,
        sigma_c,
        m,
        w_ex
    )
    return Functor{:Schelkunoff, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

@inline function (functor::Functor{:Schelkunoff})(::Val{:inner})
    return functor.routes.inner(functor.state)
end

@inline function (functor::Functor{:Schelkunoff})(::Val{:outer})
    return functor.routes.outer(functor.state)
end

@inline function (functor::Functor{:Schelkunoff})(::Val{:mutual})
    return functor.routes.mutual(functor.state)
end

@inline function (functor::Functor{:Schelkunoff})(form::Symbol)
    Base.@nospecialize form
    return form === :inner ? functor(Val(:inner)) :
           form === :outer ? functor(Val(:outer)) :
           form === :mutual ? functor(Val(:mutual)) :
           throw(ArgumentError("unknown Schelkunoff interaction: $form"))
end

@inline function schelkunoff_inner(state)
    T = typeof(state.r_in)
    if isapprox(state.r_in, zero(T); atol = eps(T))
        return zero(Complex{T})
    end

    w_in = state.m * state.r_in
    sc_in = exp(abs(real(w_in)) - state.w_ex)
    sc_ex = exp(abs(real(state.w_ex)) - w_in)
    sc = sc_in / sc_ex
    numerator = special_besselkx(0, w_in) * special_besselix(1, state.w_ex) +
                sc * special_besselix(0, w_in) * special_besselkx(1, state.w_ex)
    denominator = special_besselkx(1, w_in) * special_besselix(1, state.w_ex) -
                  sc * special_besselix(1, w_in) * special_besselkx(1, state.w_ex)
    return Complex{T}(
        (state.jω * state.mu_c / 2π) * (1 / w_in) * (numerator / denominator)
    )
end

@inline function schelkunoff_outer(state)
    T = typeof(state.r_in)
    if isapprox(state.r_in, zero(T); atol = eps(T))
        numerator = special_besselix(0, state.w_ex)
        denominator = special_besselix(1, state.w_ex)
    else
        w_in = state.m * state.r_in
        sc_in = exp(abs(real(w_in)) - state.w_ex)
        sc_ex = exp(abs(real(state.w_ex)) - w_in)
        sc = sc_in / sc_ex
        numerator = special_besselix(0, state.w_ex) * special_besselkx(1, w_in) +
                    sc * special_besselkx(0, state.w_ex) * special_besselix(1, w_in)
        denominator = special_besselix(1, state.w_ex) * special_besselkx(1, w_in) -
                      sc * special_besselkx(1, state.w_ex) * special_besselix(1, w_in)
    end
    return Complex{T}(
        (state.jω * state.mu_c / 2π) * (1 / state.w_ex) *
        (numerator / denominator)
    )
end

@inline function schelkunoff_mutual(state)
    T = typeof(state.r_in)
    if isapprox(state.r_in, zero(T); atol = eps(T))
        return zero(Complex{T})
    end

    w_in = state.m * state.r_in
    sc_in = exp(abs(real(w_in)) - state.w_ex)
    sc_ex = exp(abs(real(state.w_ex)) - w_in)
    sc = sc_in / sc_ex
    numerator = one(sc_ex) / sc_ex
    denominator = special_besselix(1, state.w_ex) * special_besselkx(1, w_in) -
                  sc * special_besselix(1, w_in) * special_besselkx(1, state.w_ex)
    return Complex{T}(
        (1 / (2π * state.r_in * state.r_ex * state.sigma_c)) *
        (numerator / denominator)
    )
end

:Schelkunoff
