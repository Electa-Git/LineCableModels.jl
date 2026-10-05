"""
$(TYPEDEF)

Select a Sommerfeld-type integral over the spatial Fourier variable on `[0, Inf)`.
The formula supplies the complete scalar integrand, including any weights and Jacobians.

$(TYPEDFIELDS)
"""
struct SpectralIntegral{F} <: AbstractFormulation
    "Complete integrand, evaluated at nonnegative real coordinates."
    f::F
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with `quadrature`, the reusable QuadGK segments of the `:quad` solution.
Coordinates have the real type of `T`, and integrand values have the type `Complex{T}`.
`buffers` returns unchanged when it holds `quadrature` already. `integrate` reads the
segments from this record.
"""
function initialize_buffers(::Type{SpectralIntegral}, ::Val{:quad}, ::Type{T}, input, plan,
        buffers) where {T}
    haskey(buffers, :quadrature) && return buffers
    R = typeof(float(nominal(one(T))))
    V = Complex{T}
    size = 128
    return merge(buffers, (quadrature = (
        segments = alloc_segbuf(R, V, R; size),
        seeds = alloc_segbuf(R, V, R; size),
        seed = alloc_segbuf(R, V, R; size = 1),
        prototype = zero(V),
        points = sizehint!(R[], size),
        mapped = sizehint!(R[], size),
        warnings = nothing),))
end

"""
$(TYPEDSIGNATURES)

Evaluate the complete callable over `[0, Inf)` with QuadGK. QuadGK takes its storage from
`buffers.quadrature` when the caller passes `buffers`. Return
`(value, estimated_error)`, where the error is the estimated absolute
quadrature error in the same units as the integral.

# Keywords

- `points`: finite nonnegative subdivisions in the callable's coordinate.
  The supplied collection is not mutated.
- `coordinate_type`: real coordinate type when no buffers are supplied.
  Defaults to `Float64`. Buffers supply their own coordinate type.
- `context`: optional caller-provided description included in numerical warnings.
- `observations`: optional trace vector receiving the native result and context.

An unmet requested target produces a warning and returns the finite result.
There is no outer retry or error-budget controller. Nonfinite integral values
and actual quadrature failures remain errors.
"""
@inline function integrate(
        integral::SpectralIntegral, ::Val{:quad}, controls, buffers = nothing;
        points = (), coordinate_type::Type{C} = Float64, context = nothing,
        observations = nothing) where {C <: AbstractFloat}
    quadrature=buffers===nothing ? nothing : buffers.quadrature
    R=quadrature===nothing ? coordinate_type : eltype(quadrature.points)
    subdivisions=quadrature===nothing ? R[] : empty!(quadrature.points)
    push!(subdivisions, zero(R))
    for point in points
        point isa Real && isfinite(point) && point>=0 ||
            throw(DomainError(point, "integration points must be finite nonnegative real coordinates"))
        push!(subdivisions, R(point))
    end
    all(isfinite, subdivisions) ||
        throw(DomainError(subdivisions, "integration points must be finite in the coordinate type"))
    unique!(sort!(subdivisions; alg = Base.Sort.QuickSort))
    mapped=quadrature===nothing ? R[] : empty!(quadrature.mapped)
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
    error=if quadrature===nothing
        quadgk(f, mapped; rtol = controls.rtol, atol = controls.atol,
            maxevals = controls.maxevals, norm = numerical_magnitude)
    else
        # Public QuadGK evaluation creates reusable seeds in the same mapped
        # coordinate as f. No dependency-private segment constructors are used.
        seeds=empty!(quadrature.seeds)
        sample=quadrature.prototype
        for i in 1:(length(mapped) - 1)
            quadgk(_->sample, mapped[i], mapped[i + 1];
                segbuf = quadrature.seed, norm = numerical_magnitude)
            append!(seeds, quadrature.seed)
        end
        quadgk(f, zero(R), one(R); rtol = controls.rtol, atol = controls.atol,
            maxevals = controls.maxevals, segbuf = quadrature.segments,
            eval_segbuf = seeds, norm = numerical_magnitude)
    end
    isfinite(value) || throw(DomainError(value, "QuadGK returned a nonfinite integral"))
    record_integral!(observations, quadrature === nothing ? nothing : quadrature.warnings,
        value, error, controls, context)
    return value, error
end

# Actual and reused integrals report the same estimates and requested controls.
function record_integral!(observations, warnings, value, error, controls, context)
    observations === nothing ||
        push!(observations, (; context, value, estimated_error = error))
    target=max(controls.atol, controls.rtol*numerical_magnitude(value))
    if !isfinite(error) || error>target
        warnings === nothing ||
            push!(warnings, (; context, value, estimated_error = error, controls))
        @warn "QuadGK returned an estimated error above the requested target" value estimated_error=error target rtol=controls.rtol atol=controls.atol maxevals=controls.maxevals context
    end
    return nothing
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
