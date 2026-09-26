"""
$(TYPEDEF)

Represent the integral of one complete scalar callable on `[0, Inf)`.
The formula supplies the entire integrand, including any weights and Jacobians.
Construction does not evaluate the callable.

$(TYPEDFIELDS)
"""
struct SpectralIntegral{F}
    "Complete integrand, evaluated at nonnegative real coordinates."
    f::F
end

"""
$(TYPEDSIGNATURES)

Allocate reusable QuadGK segments for real coordinate type `R` and scalar
integrand value type `V`. No physical model or sampled values are retained.
`size` is the initial segment capacity; zero leaves all numerical arrays empty.
"""
function integration_workspace(::Type{R}, ::Type{V} = Complex{R};
        size::Integer = 128) where {
        R <: AbstractFloat, V <: Number}
    return (
        segments = alloc_segbuf(R, V, R; size),
        seeds = alloc_segbuf(R, V, R; size),
        seed = alloc_segbuf(R, V, R; size = min(size, 1)),
        prototype = zero(V),
        points = sizehint!(R[], size),
        mapped = sizehint!(R[], size))
end

function initialize_buffers(::Val{:quad}, ::Type{T}, input, invariants, buffers) where {T}
    isempty(buffers.quadrature.segments) || return buffers
    R = typeof(float(nominal(one(T))))
    return merge(buffers, (quadrature = integration_workspace(R, Complex{T}),))
end

"""
$(TYPEDSIGNATURES)

Evaluate the complete callable over `[0, Inf)` with QuadGK. Return
`(value, estimated_error)`, where the error is the estimated absolute
quadrature error in the same units as the integral.

# Keywords

- `points`: Finite nonnegative subdivisions in the callable's coordinate;
  the supplied collection is not mutated.
- `coordinate_type`: Real coordinate type when no workspace is supplied;
  defaults to `Float64`. A workspace supplies its own coordinate type.
- `context`: Optional caller-provided description included in numerical warnings.
- `observations`: Optional trace vector receiving the native result and context.

An unmet requested target produces a warning and returns the finite result.
There is no outer retry or error-budget controller. Nonfinite integral values
and actual quadrature failures remain errors.
"""
@inline function integrate(
        ::Val{:quad}, integral::SpectralIntegral, controls, workspace = nothing;
        points = (), coordinate_type::Type{C} = Float64, context = nothing,
        observations = nothing) where {C <: AbstractFloat}
    R=workspace===nothing ? coordinate_type : eltype(workspace.points)
    subdivisions=workspace===nothing ? R[] : empty!(workspace.points)
    push!(subdivisions, zero(R))
    for point in points
        point isa Real && isfinite(point) && point>=0 ||
            throw(DomainError(point, "integration points must be finite nonnegative real coordinates"))
        push!(subdivisions, R(point))
    end
    all(isfinite, subdivisions) ||
        throw(DomainError(subdivisions, "integration points must be finite in the coordinate type"))
    unique!(sort!(subdivisions; alg = Base.Sort.QuickSort))
    mapped=workspace===nothing ? R[] : empty!(workspace.mapped)
    for point in subdivisions
        push!(mapped, point/(one(R)+point))
    end
    push!(mapped, one(R))
    unique!(sort!(mapped; alg = Base.Sort.QuickSort))
    f=t->begin
        den=inv(one(R)-t)
        integral.f(t*den)*den^2
    end
    value,
    error=if workspace===nothing
        quadgk(f, mapped; rtol = controls.rtol, atol = controls.atol,
            maxevals = controls.maxevals, norm = numerical_magnitude)
    else
        # Public QuadGK evaluation creates reusable seeds in the same mapped
        # coordinate as f. No dependency-private segment constructors are used.
        seeds=empty!(workspace.seeds)
        sample=workspace.prototype
        for i in 1:(length(mapped) - 1)
            quadgk(_->sample, mapped[i], mapped[i + 1];
                segbuf = workspace.seed, norm = numerical_magnitude)
            append!(seeds, workspace.seed)
        end
        quadgk(f, zero(R), one(R); rtol = controls.rtol, atol = controls.atol,
            maxevals = controls.maxevals, segbuf = workspace.segments,
            eval_segbuf = seeds, norm = numerical_magnitude)
    end
    isfinite(value) || throw(DomainError(value, "QuadGK returned a nonfinite integral"))
    observations===nothing ||
        push!(observations, (; context, value, estimated_error = error))
    target=max(controls.atol, controls.rtol*numerical_magnitude(value))
    if !isfinite(error) || error>target
        @warn "QuadGK returned an estimated error above the requested target" value estimated_error=error target rtol=controls.rtol atol=controls.atol maxevals=controls.maxevals context
    end
    return value, error
end

"""
$(TYPEDSIGNATURES)

Normalize adaptive Gauss–Kronrod integration controls for `method=:quad`.
Accepted controls are `rtol`, `atol`, and `maxevals`.
"""
function formulation_options(::Type{SpectralIntegral}, values::NamedTuple)
    isempty(setdiff(keys(values), (:method, :options))) ||
        throw(ArgumentError("integration accepts only method and options"))
    method = get(values, :method, :quad)
    method === :quad ||
        throw(ArgumentError("integration method must be :quad"))
    supplied = get(values, :options, (;))
    supplied isa NamedTuple ||
        throw(ArgumentError("integration options must be a named tuple"))
    defaults = (rtol = 1e-8, atol = 0.0, maxevals = 10^7)
    unknown_controls = setdiff(keys(supplied), keys(defaults))
    isempty(unknown_controls) ||
        throw(ArgumentError("unknown :$method integration controls: $(collect(unknown_controls))"))
    controls = merge(defaults, supplied)
    for key in (:rtol, :atol)
        value = getproperty(controls, key)
        value isa Real && isfinite(value) && value >= 0 ||
            throw(ArgumentError("$key must be finite and nonnegative"))
    end
    controls.rtol > 0 || controls.atol > 0 ||
        throw(ArgumentError("at least one integration tolerance must be positive"))
    for key in setdiff(keys(controls), (:rtol, :atol))
        value = getproperty(controls, key)
        value isa Integer && !(value isa Bool) && value > 0 ||
            throw(ArgumentError("$key must be a positive integer"))
    end
    return (method = Val(method), options = controls)
end

function formulation_options(::FormulaMethod, ::Val{:integration}, defaults::NamedTuple,
        supplied::NamedTuple)
    isempty(setdiff(keys(supplied), (:method, :options))) ||
        throw(ArgumentError("integration accepts only method and options"))
    get(supplied, :options, (;)) isa NamedTuple ||
        throw(ArgumentError("integration options must be a NamedTuple"))
    method = get(supplied, :method, defaults.method)
    controls = merge(defaults.options, get(supplied, :options, (;)))
    return formulation_options(SpectralIntegral, (; method, options = controls))
end
