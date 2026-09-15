"""
$(TYPEDSIGNATURES)

**Identification.** Exact cylindrical surface impedances for a solid or hollow
round conductor, after S. A. Schelkunoff (1934).

**Reference.** S. A. Schelkunoff, “The Electromagnetic Theory of Coaxial
Transmission Lines and Cylindrical Shields,” *Bell System Technical Journal*,
13, 532–579, 1934.
"""
function description(::Type{<:Formula{:schelkunoff1934}}; compact::Bool=false)
    compact ? "Schelkunoff" : "Schelkunoff exact round-conductor surface impedances (1934)"
end

"Delegate the normalized Schelkunoff route to the package default implementation."
surface_impedance_state(
        ::Val{:schelkunoff1934}, args...
) = surface_impedance_state(Val(:default), args...)

@inline function (formula::Formula{:schelkunoff1934})(
        r_in::T,
        r_ex::T,
        rho_c::T,
        mur_c::T,
        jω::Complex{T}
) where {T <: Real}
    state = surface_impedance_state(
        Val(:schelkunoff1934), r_in, r_ex, rho_c, mur_c, jω
    )
    return Functor{:schelkunoff1934, typeof(formula.binding), typeof(formula.hooks),
        typeof(state), typeof(formula.options)}(
        formula.binding,
        formula.hooks,
        state,
        formula.options
    )
end

@inline function internal_impedance(
        ::Val{:schelkunoff1934}, kind::Val, functor, workspace
)
    internal_impedance(Val(:default), kind, functor, workspace)
end

computation_options(::FormulaMethod{:schelkunoff1934, typeof(internal_impedance)}) = (;)

:schelkunoff1934
