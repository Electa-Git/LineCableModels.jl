function _grid_item(value)
    value isa Real && return TextDisplay.value(value)
    return sprint(show, value; context = :compact => true, sizehint = 64)
end

function _display_grid_values(values)
    count = length(values)
    count == 0 && return "empty"
    visible = min(count, 4)
    items = join((_grid_item(values[index]) for index in 1:visible), ", ")
    return count > visible ? "$items, … ($count values)" : items
end

TextDisplay.name(::Type{<:UncertainValue}) = "UncertainValue"
function Base.summary(io::IO, ::UncertainValue)
    print(io, "Uncertain value")
end
function Base.show(io::IO, value::UncertainValue)
    print(io, "UncertainValue(", _grid_item(value.nominal), " ± ", _grid_item(value.sigma), ")")
end
function Base.show(io::IO, ::MIME"text/plain", value::UncertainValue)
    show(io, value)
end

TextDisplay.name(::Type{<:AbsoluteError}) = "AbsoluteError"
function Base.summary(io::IO, error::AbsoluteError)
    print(io, "Absolute error with $(length(error.vals)) values")
end
function Base.show(io::IO, error::AbsoluteError)
    print(io, "AbsoluteError(", _display_grid_values(error.vals), ")")
end
Base.show(io::IO, ::MIME"text/plain", error::AbsoluteError) = show(io, error)

function _grid_kind(::DeterministicGrid)
    return "Grid"
end
function _grid_kind(::RelativeGrid)
    return "Relative-error Grid"
end
function _grid_kind(::AbsoluteGrid)
    return "Absolute-error Grid"
end

function _grid_nominals(grid::DeterministicGrid)
    return grid.vals
end
function _grid_nominals(grid::Union{RelativeGrid, AbsoluteGrid})
    return grid.vals
end
function _grid_errors(grid::RelativeGrid)
    return "relative error  $(_display_grid_values(grid.rel_err)) %"
end
function _grid_errors(grid::AbsoluteGrid)
    return "absolute error  $(_display_grid_values(grid.abs_err))"
end

TextDisplay.name(::Type{<:DeterministicGrid}) = "Grid"
TextDisplay.name(::Type{<:RelativeGrid}) = "Relative-error Grid"
TextDisplay.name(::Type{<:AbsoluteGrid}) = "Absolute-error Grid"
function Base.summary(io::IO, grid::AbstractGrid)
    print(io, _grid_kind(grid), " with ", length(grid), " points")
end
function Base.show(io::IO, grid::AbstractGrid)
    print(io, _grid_kind(grid), "(", _display_grid_values(_grid_nominals(grid)))
    grid isa DeterministicGrid || print(io, "; ", _grid_errors(grid))
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", grid::AbstractGrid)
    get(io, :compact, false) && return show(io, grid)
    children = Any[(label = "values  $(_display_grid_values(_grid_nominals(grid)))", noun = "values")]
    grid isa DeterministicGrid || push!(children, (label = _grid_errors(grid), noun = "values"))
    return TextDisplay.tree(io, "$(_grid_kind(grid)) · $(length(grid)) points", Tuple(children))
end

function _semantic_target(::Type{Target}) where {Target}
    applicable(TextDisplay.name, Target) && return TextDisplay.name(Target)
    return String(nameof(Target))
end

function _space_axis(source::Gridspace)
    return sprint(show, source; context = :compact => true, sizehint = 96)
end
function _space_axis(source::AbstractGrid)
    return _display_grid_values(_grid_nominals(source))
end

TextDisplay.name(::Type{<:Gridspace{Target}}) where {Target} =
    string(_semantic_target(Target), " parameter space")
function Base.summary(io::IO, space::Gridspace{Target}) where {Target}
    print(io, _semantic_target(Target), " parameter space with ", length(space), " points")
end
function Base.show(io::IO, space::Gridspace{Target}) where {Target}
    print(io, _semantic_target(Target), " parameter space · ", length(space), " points")
end
function Base.show(io::IO, ::MIME"text/plain", space::Gridspace{Target}) where {Target}
    get(io, :compact, false) && return show(io, space)
    varying = Tuple((index, source) for (index, source) in enumerate(space.grids)
        if length(source) > 1)
    fixed = length(space.grids) - length(varying)
    children = Any[
        (label = "axis $index  $(_space_axis(source))", noun = "axes")
        for (index, source) in varying
    ]
    fixed == 0 || push!(children, (
        label = "$fixed fixed $(fixed == 1 ? "input" : "inputs")",
        noun = "axes",
    ))
    target = _semantic_target(Target)
    return TextDisplay.tree(
        io,
        "$target parameter space · $(length(varying)) varying axes · $(length(space)) points",
        Tuple(children);
        noun = "axes"
    )
end

TextDisplay.name(::Type{<:Combinatorial}) = "Combinatorial"
Base.summary(io::IO, ::Combinatorial) = print(io, "Combinatorial formulation")
function Base.show(io::IO, formulation::Combinatorial)
    print(io, "Combinatorial(")
    show(IOContext(io, :compact => true), formulation.inner)
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", formulation::Combinatorial) = show(io, formulation)

TextDisplay.name(::Type{<:ParametricProblem}) = "ParametricProblem"
Base.summary(io::IO, ::ParametricProblem) = print(io, "Parametric problem")
function Base.show(io::IO, problem::ParametricProblem)
    print(io, "ParametricProblem(")
    show(IOContext(io, :compact => true), problem.space)
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", problem::ParametricProblem)
    get(io, :compact, false) && return show(io, problem)
    return TextDisplay.tree(io, "ParametricProblem", (
        (label = "space    $(sprint(show, problem.space; context = :compact => true))", noun = "fields"),
        (label = "options  $(length(problem.options)) entries", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:ParametricResult}) = "ParametricResult"
function Base.summary(io::IO, result::ParametricResult)
    print(io, "Parametric result with $(length(result)) values")
end
function Base.show(io::IO, result::ParametricResult)
    print(io, "ParametricResult($(length(result)) values)")
end
function Base.show(io::IO, ::MIME"text/plain", result::ParametricResult)
    get(io, :compact, false) && return show(io, result)
    value_type = isempty(result.values) ? "none" : _semantic_target(eltype(result.values))
    return TextDisplay.tree(io, "ParametricResult · $(length(result)) values", (
        (label = "result type  $value_type", noun = "fields"),
        (label = "details      $(length(result.details)) entries", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:WireEstimate}) = "WireEstimate"
function Base.summary(io::IO, estimate::WireEstimate)
    print(io, "Wire estimate with $(length(estimate)) candidates")
end
function Base.show(io::IO, estimate::WireEstimate)
    print(io, "WireEstimate(", estimate.status, "; candidates=", length(estimate), ")")
end
function Base.show(io::IO, ::MIME"text/plain", estimate::WireEstimate)
    get(io, :compact, false) && return show(io, estimate)
    candidates = Tuple((
        label = sprint(show, candidate; context = :compact => true, sizehint = 96),
        noun = "candidates",
    ) for candidate in estimate.candidates)
    children = Any[
        (label = "target      $(TextDisplay.value(estimate.target))", noun = "fields"),
        (label = "status      $(estimate.status)", noun = "fields"),
        (label = "candidates", children = candidates, noun = "candidates"),
    ]
    isempty(estimate.reasons) || push!(children, (
        label = "reasons",
        children = Tuple((label = reason, noun = "reasons") for reason in estimate.reasons),
        noun = "reasons",
    ))
    return TextDisplay.tree(io, "Wire estimate", Tuple(children))
end

function Base.summary(io::IO, pattern::WirePatterns.HexaPattern)
    print(io, "Strand pattern with $(pattern.wires) wires")
end
function Base.show(io::IO, pattern::WirePatterns.HexaPattern)
    print(io, "StrandPattern(layers=", pattern.layers, ", wires=", pattern.wires,
        ", d=", TextDisplay.engineering(pattern.wire_diameter_m, :meter), ")")
end
Base.show(io::IO, ::MIME"text/plain", pattern::WirePatterns.HexaPattern) = show(io, pattern)

function Base.summary(io::IO, pattern::WirePatterns.ScreenPattern)
    print(io, "Screen pattern with $(pattern.wires) wires")
end
function Base.show(io::IO, pattern::WirePatterns.ScreenPattern)
    print(io, "ScreenPattern(wires=", pattern.wires,
        ", d=", TextDisplay.engineering(pattern.wire_diameter_m, :meter),
        ", coverage=", TextDisplay.value(pattern.coverage_pct), " %)")
end
Base.show(io::IO, ::MIME"text/plain", pattern::WirePatterns.ScreenPattern) = show(io, pattern)
