"""
$(TYPEDEF)

Supertype of the rules that combine the finite sources of a [`Gridspace`](@ref):
[`Product`](@ref) and [`Zip`](@ref). The user-facing spelling `combine=:product` or
`combine=:zip` selects one through `Composition(name)`.
"""
abstract type Composition end

"""
$(TYPEDEF)

Combine every point of each source with every point of the others. The first source
varies fastest.
"""
struct Product <: Composition end

"""
$(TYPEDEF)

Pair the points of equal-length sources, broadcasting sources with one point. Zip pairs
the direct parameters of one call. It refuses a nested parameter space that varies. Paired
complete objects are a `Grid` of those objects.
"""
struct Zip <: Composition end

"""
$(TYPEDSIGNATURES)

Return the composition named by `combine`: `:product` or `:zip`, or a composition itself.

# Errors

- Throws `ArgumentError` for any other name.
"""
Composition(composition::Composition) = composition
Composition(name::Symbol) = Composition(Val(name))
Composition(::Val{:product}) = Product()
Composition(::Val{:zip}) = Zip()
function Composition(::Val{name}) where {name}
    throw(ArgumentError("combine must be :product or :zip; got :$name"))
end

"""
    Gridpoint{Target}(build, args)

Store one selected but unresolved argument tuple from a finite parameter space.
`Target` preserves the semantic object or problem family that the point will
materialize, allowing computation dispatch to consume a scalar point without
first discarding its target identity.
"""
struct Gridpoint{Target, F, A <: Tuple}
    "Scalar constructor or lowering function selected by the finite space."
    build::F
    "Selected argument tuple, with any nested target-bearing points retained."
    args::A
end

# A type constructor is stored as `Type{B}`, so calls through it specialize.
function Gridpoint{Target}(build::Type{B}, args::A) where {Target, B, A <: Tuple}
    return Gridpoint{Target, Type{B}, A}(build, args)
end
function Gridpoint{Target}(build, args::A) where {Target, A <: Tuple}
    return Gridpoint{Target, typeof(build), A}(build, args)
end

"""
$(TYPEDEF)

Represent a typed finite space assembled from explicit [`Grid`](@ref), nested
`Gridspace` or admitted completed result-space sources. The composition `C`, a
[`Product`](@ref) or a [`Zip`](@ref), is part of the type. `Target` identifies the semantic
result family and is the element type of the space. Each point is built by `build` when
the space is iterated, so its concrete type is the type that `build` returns.

$(TYPEDFIELDS)
"""
struct Gridspace{Target, F, G <: Tuple, C <: Composition}
    "Callable that constructs `Target` from one selected argument tuple."
    build::F

    "Explicit Grid or nested Gridspace sources."
    grids::G

    function Gridspace{Target, F, G, C}(build, grids) where {
            Target, F, G <: Tuple, C <: Composition}
        return validate(new{Target, F, G, C}(build, grids))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a lazy space from explicit finite sources. Product composition uses
Julia's product order, with the first source varying fastest. Zip composition
pairs equal-length sources and broadcasts singleton sources. A completed
result space is admitted only through a result-owner transport method of the
form `Gridspace{Target}(source)`.

# Errors

- Throws `ArgumentError` when `combine` is unsupported, when a source is not a
  `Grid`, a nested `Gridspace` or a completed result space, or when a zip includes
  a nested parameter space.
- Throws `DimensionMismatch` when zip source lengths are incompatible.
"""
function Gridspace{Target}(build, grids::Tuple; combine = :product) where {Target}
    return Gridspace{Target}(build, grids, Composition(combine))
end
function Gridspace{Target}(build::Type{B}, grids::G, ::C) where {
        Target, B, G <: Tuple, C <: Composition}
    return Gridspace{Target, Type{B}, G, C}(build, grids)
end
function Gridspace{Target}(build, grids::G, ::C) where {Target, G <: Tuple, C <: Composition}
    return Gridspace{Target, typeof(build), G, C}(build, grids)
end

function Gridspace{Target}(grids::Tuple; combine = :product) where {Target}
    Gridspace{Target}(Target, grids; combine)
end

"""
$(TYPEDSIGNATURES)

Reject result-to-problem transports whose result owner has not defined how to
construct the requested target-bearing `Gridspace`.
"""
function Gridspace{Target}(source::AbstractResultSpace) where {Target}
    throw(ArgumentError(
        "Gridspace transport from $(nameof(typeof(source))) to problem $(nameof(Target)) " *
        "is not defined. Define Gridspace{$(nameof(Target))}(::$(nameof(typeof(source)))) " *
        "for that result family. See `?Gridspace` and " *
        "https://electa-git.github.io/LineCableModels.jl/dev/gridspace/" *
        "#Transporting-completed-result-spaces",
    ))
end

"""
$(TYPEDSIGNATURES)

Return the finite source that `values` contributes to a space. A `Grid` or a `Gridspace`
is already a source. A tuple, vector, named tuple or pair that holds a source becomes the
space of its reconstructions, combined by `combine` like the call it belongs to. Any
other value is a source with one point: a scalar is a space of cardinality one.
"""
function Gridspace(
        values::Union{AbstractGrid, Gridspace, Tuple, AbstractVector, NamedTuple, Pair};
        combine = :product)
    return Gridspace(values, Composition(combine))
end
Gridspace(source::Union{AbstractGrid, Gridspace}, ::Composition) = source
Gridspace(value, ::Composition) = DeterministicGrid((value,))
function Gridspace(values::Tuple, composition::Composition)
    has_grid(values) || return DeterministicGrid((values,))
    sources = map(Base.Fix2(Gridspace, composition), values)
    return Gridspace{Tuple}(tuple, sources, composition)
end
function Gridspace(values::AbstractVector, composition::Composition)
    has_grid(values) || return DeterministicGrid((values,))
    sources = Tuple(map(Base.Fix2(Gridspace, composition), values))
    return Gridspace{Vector}(Base.vect, sources, composition)
end
# The names of a named tuple are a type parameter of its build, so a point of it is isbits.
function Gridspace(values::NamedTuple{K}, composition::Composition) where {K}
    has_grid(values) || return DeterministicGrid((values,))
    sources = map(Base.Fix2(Gridspace, composition), Tuple(values))
    return Gridspace{NamedTuple{K}}((selected...) -> NamedTuple{K}(selected), sources,
        composition)
end
function Gridspace(pair::Pair, composition::Composition)
    has_grid(pair) || return DeterministicGrid((pair,))
    sources = (Gridspace(first(pair), composition), Gridspace(last(pair), composition))
    return Gridspace{Pair}(Pair, sources, composition)
end

Grid(space::Gridspace) = space

"""
$(TYPEDSIGNATURES)

Return whether `value` holds a `Grid` or a `Gridspace`, directly or at any depth of the
tuples, vectors, named tuples and pairs that contain it. A dictionary is one value.
"""
has_grid(::AbstractGrid) = true
has_grid(::Gridspace) = true
# A tuple is walked head first, so each step has one concrete element type and nested
# named tuples of concrete types fold to a constant.
has_grid(values::Tuple) = has_grid(first(values)) | has_grid(Base.tail(values))
has_grid(::Tuple{}) = false
has_grid(values::NamedTuple) = has_grid(Tuple(values))
has_grid(values::AbstractVector) = any(has_grid, values)
has_grid(pair::Pair) = has_grid(first(pair)) || has_grid(last(pair))
has_grid(::Any) = false

"""
    parameterize(Target, caller, values; combine=:product)

Return a `Gridspace{Target}` that reconstructs the call point by point when `values` contain
a `Grid` or a `Gridspace`, and otherwise invoke `caller` directly. A tuple, vector,
named tuple or pair among `values` that holds a source is varied at its depth, with the
same composition.
"""
Base.@constprop :aggressive function parameterize(
        ::Type{Target}, caller, values::Tuple; combine = :product
) where {Target}
    has_grid(values) || return caller(values...)
    composition = Composition(combine)
    sources = map(Base.Fix2(Gridspace, composition), values)
    return Gridspace{Target}(caller, sources, composition)
end

"""
    parameterize(Target, caller, arguments, keywords; combine=:product)

Lift the call `caller(arguments...; keywords...)` over its arguments and the values of its
keywords, all as direct sources of one space. Each point calls `caller` with its
positional arguments and its keywords again.
"""
Base.@constprop :aggressive function parameterize(::Type{Target}, caller, arguments::Tuple{Vararg{Any, N}},
        keywords::NamedTuple{K}; combine = :product) where {Target, N, K}
    split = (selected...) -> caller(ntuple(index -> selected[index], Val(N))...;
        NamedTuple{K}(ntuple(index -> selected[N + index], Val(length(K))))...)
    return parameterize(Target, split, (arguments..., Tuple(keywords)...); combine)
end

"""
    parameterize(Target, arguments; combine=:product)

Lift a call of the constructor of `Target` whose `arguments` hold a `Grid` or a
`Gridspace`. Each owner of an input family adds one method of its constructors that
calls this, `(::Type{T})(arguments...; combine) where {T <: Family}`. Dispatch selects it
only when no scalar method accepts the arguments.

# Errors

- Throws the `MethodError` of the call when no argument varies.
"""
function parameterize(::Type{Target}, arguments::Tuple; combine = :product) where {Target}
    has_grid(arguments) || throw(MethodError(Target, arguments))
    return parameterize(Target, Target, arguments; combine = Composition(combine))
end

# The constructors of the input families declared in Commons lift through the grammar.
function (::Type{T})(arguments...; combine = :product) where {T <: AbstractProblemDefinition}
    return parameterize(T, arguments; combine)
end
function (::Type{T})(arguments...; combine = :product) where {T <: AbstractFormulation}
    return parameterize(T, arguments; combine)
end
# The record type itself, so a call of the parametric constructor reaches only the record.
function FormulationOptions(arguments...; combine = :product)
    return parameterize(FormulationOptions, arguments; combine)
end

"""
    target(builder, subject, arguments...; keywords...)
    target(builder, ::Type)

Return the type of the object that a call of `builder` builds, the type of the points of the
call's `Gridspace` when it lifts. It is the `Target` of `Builder{Target}` unless the owner of
the builder adds a method for the scalar type of the call's subject, its first argument, as
`target(::AtBuilder, ::Type{<:AbstractCablePart}) = AssemblyMember`. A subject that is a
`Grid` or a `Gridspace` selects by the type of its points. The methods that take a subject
receive the call's keywords, and these defaults ignore them. An owner whose target follows a
keyword adds one method that receives the keyword and passes its value on as the subject.
That method types its own subject, so that the value it passes on does not return to it.
`stranded` does this for its `shape` keyword.
"""
target(::Builder{Target}; keywords...) where {Target} = Target
target(::Builder{Target}, ::Type) where {Target} = Target
target(builder::Builder, subject, arguments...; keywords...) = target(builder, typeof(subject))
function target(builder::Builder, subject::Union{AbstractGrid, Gridspace}, arguments...;
        keywords...)
    return target(builder, eltype(subject))
end
# An empty source has no points to select a target by.
target(::Builder{Target}, ::DeterministicGrid{Tuple{}}, arguments...; keywords...) where {Target} =
    Target

# Every builder lifts here, over all its arguments and keywords, so that `combine` composes
# all of them. Each point builds from scalars.
Base.@constprop :aggressive function (builder::Builder)(arguments...; combine = :product, keywords...)
    caller = Base.Fix1(build, builder)
    return parameterize(target(builder, arguments...; keywords...), caller, arguments,
        NamedTuple(keywords); combine)
end

# Each source is checked for the composition of the space, and a zip checks the source
# cardinalities too.
function validate(space::Gridspace{<:Any, <:Any, <:Tuple, C}) where {C <: Composition}
    foreach(Base.Fix2(validate, C()), space.grids)
    return space
end
function validate(space::Gridspace{<:Any, <:Any, <:Tuple, Zip})
    foreach(Base.Fix2(validate, Zip()), space.grids)
    counts = map(length, space.grids)
    target_count = maximum(counts; init = 0)
    all(count -> count == 1 || count == target_count, counts) || throw(
        DimensionMismatch(
        "zip sources must have equal cardinality or be singletons; got $counts",
    ),
    )
    return space
end

# A source of a space, for its composition. A zip pairs the direct parameters of one call
# and the containers that hold them. A nested space of one point pairs as its only point,
# as a product would. A zip refuses only a nested space that varies.
validate(source::Union{AbstractGrid, Gridspace, AbstractResultSpace}, ::Product) = source
validate(source::Union{AbstractGrid, AbstractResultSpace}, ::Zip) = source
function validate(
        source::Gridspace{<:Union{Tuple, NamedTuple, AbstractVector, Pair}}, ::Zip)
    return source
end
function validate(source::Gridspace{Target}, ::Zip) where {Target}
    isone(length(source)) && return source
    throw(ArgumentError(
        "zip pairs the direct parameters of one call; a nested parameter space of " *
        "$(length(source)) $(nameof(Target)) points varies on its own and multiplies " *
        "with them. To pair complete objects, pass a Grid of them, as in " *
        "Grid((first, second))",
    ))
end
function validate(source, ::Composition)
    throw(ArgumentError(
        "Gridspace sources must be Grid, Gridspace, or completed result-space values; " *
        "got $(nameof(typeof(source)))",
    ))
end

"""
Return the unresolved points of a finite source, in its order. A `Grid` gives its values
and a completed result space its results. A Gridspace gives its points lazily. With an
`index`, return that one point. A source is indexed from 1 to its `length`.
"""
points(grid::AbstractGrid) = grid
points(result::AbstractResultSpace) = result
points(space::Gridspace) = (points(space, index) for index in 1:length(space))
points(source::Union{AbstractGrid, AbstractResultSpace}, index::Integer) = source[index]
# The points of `sources` at `indices`, one per source. The tuple is unrolled, so a nested
# space with more sources than its parent still infers: a recursive tuple `map` at both
# levels would widen the inner one.
@generated function points(sources::Tuple, indices::Tuple)
    return Expr(:tuple, (:(points(sources[$k], indices[$k])) for k in 1:fieldcount(sources))...)
end
# A product point takes from each source the point of its index, the first source fastest.
function points(
        space::Gridspace{Target, <:Any, <:Tuple, Product}, index::Integer) where {Target}
    selected = Tuple(CartesianIndices(map(length, space.grids))[index])
    return Gridpoint{Target}(space.build, points(space.grids, selected))
end
# A zip point takes the same index from every source, and the one point of a singleton.
function points(space::Gridspace{Target, <:Any, <:Tuple, Zip}, index::Integer) where {Target}
    1 <= index <= length(space) || throw(BoundsError(space, index))
    selected = map(source -> isone(length(source)) ? 1 : index, space.grids)
    return Gridpoint{Target}(space.build, points(space.grids, selected))
end

"""
Return a value unchanged during deterministic point materialization.
"""
materialize(value) = value

function materialize(value::UncertainValue)
    throw(ArgumentError(
        "direct materialization of an uncertainty-bearing Gridspace requires " *
        "Measurements.jl; load it with `using Measurements` before " *
        "iteration or MonteCarlo computation",
    ))
end

"""
Recursively materialize a selected Gridspace point.
"""
function materialize(point::Gridpoint)
    point.build(map(materialize, point.args)...)
end

# A problem's constructor validates it, so the point build is the checked problem.
function materialize(
        point::Gridpoint{Target}
) where {Target <: AbstractProblemDefinition}
    return point.build(map(materialize, point.args)...)::Target
end

"""
Return a deterministic value unchanged during stochastic realization.
"""
realize(::Random.AbstractRNG, value, _) = value

function realize(rng::Random.AbstractRNG, value::UncertainValue, distribution)
    rand(rng, value; distribution)
end

"""
Draw the arguments of a selected Gridspace point without invoking its builder.
"""
function realize_arguments(rng::Random.AbstractRNG, point::Gridpoint, distribution)
    return map(value -> realize(rng, value, distribution), point.args)
end

"""
Build a selected Gridspace point from an already realized argument tuple.
"""
realize(point::Gridpoint, arguments::Tuple) = point.build(arguments...)

function realize(
        point::Gridpoint{Target},
        arguments::Tuple
) where {Target <: AbstractProblemDefinition}
    return point.build(arguments...)::Target
end

"""
Recursively realize a selected Gridspace point using the caller's RNG.
"""
function realize(rng::Random.AbstractRNG, point::Gridpoint, distribution)
    return realize(point, realize_arguments(rng, point, distribution))
end

function Base.iterate(space::Gridspace{Target}, state...) where {Target}
    item = iterate(points(space), state...)
    item === nothing && return nothing
    point, next_state = item
    return materialize(point)::Target, next_state
end

#! explicit-imports: off
# Base's iterator trait protocol exposes these values without public bindings.
Base.IteratorSize(::Type{<:Gridspace}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:Gridspace}) = Base.EltypeUnknown()
#! explicit-imports: on
Base.eltype(::Type{<:Gridspace{Target}}) where {Target} = Target
# The source lengths are mapped first, so the length of a nested space infers.
function Base.length(space::Gridspace{<:Any, <:Any, <:Tuple, Product})
    prod(map(length, space.grids); init = 1)
end
Base.length(::Gridspace{<:Any, <:Any, Tuple{}, Zip}) = 1
Base.length(space::Gridspace{<:Any, <:Any, <:Tuple, Zip}) = maximum(map(length, space.grids))
Base.size(space::Gridspace) = (length(space),)

"""
Draw one unresolved point, then realize its arguments through the owned construction path.
"""
function Base.rand(
        rng::Random.AbstractRNG,
        space::Gridspace{Target};
        distribution = :normal
) where {Target}
    iszero(length(space)) && throw(ArgumentError("cannot sample an empty Gridspace"))
    offset = rand(rng, 0:(length(space) - 1))
    return realize(rng, points(space, offset + 1), distribution)::Target
end

function Base.rand(space::Gridspace; distribution = :normal)
    return rand(Random.default_rng(), space; distribution)
end

"""
Return whether a value structurally contains an uncertainty descriptor.
"""
has_uncertainty(::UncertainValue) = true
has_uncertainty(grid::AbstractUncertainGrid) = true
has_uncertainty(grid::DeterministicGrid) = any(has_uncertainty, grid.vals)
has_uncertainty(point::Gridpoint) = any(has_uncertainty, point.args)
has_uncertainty(space::Gridspace) = any(has_uncertainty, space.grids)
has_uncertainty(::Any) = false
