TextDisplay.name(::Type{<:LineModelProblem}) = "LineModelProblem"
function Base.summary(io::IO, problem::LineModelProblem)
    print(io, "Line-model fitting problem, ", size(problem.segment.parameters.Z, 1), " terminals")
end
function Base.show(io::IO, problem::LineModelProblem)
    print(io, "LineModelProblem(", size(problem.segment.parameters.Z, 1),
        " terminals, length=", line_length(problem.segment), " m)")
end
function Base.show(io::IO, ::MIME"text/plain", problem::LineModelProblem)
    get(io, :compact, false) && return show(io, problem)
    return TextDisplay.tree(io, "Line-model fitting problem",
        ((label="terminals  $(size(problem.segment.parameters.Z, 1))", noun="fields"),
         (label="length  $(line_length(problem.segment)) m", noun="fields"),
         (label="frequencies  $(length(frequencies(problem.segment))) samples", noun="fields")))
end

TextDisplay.name(::Type{<:LineModelFormulation}) = "LineModelFormulation"
function Base.summary(io::IO, formulation::LineModelFormulation)
    print(io, "Line-model formulation, :", formula_id(formulation.formula))
end
function Base.show(io::IO, formulation::LineModelFormulation)
    print(io, "LineModelFormulation(:", formula_id(formulation.formula), ")")
end
function Base.show(io::IO, ::MIME"text/plain", formulation::LineModelFormulation)
    get(io, :compact, false) && return show(io, formulation)
    return TextDisplay.tree(io, "Line-model formulation",
        ((label="formula  :$(formula_id(formulation.formula))", noun="fields"),
         (label="controls  $(length(formulation_options(formulation.formula).data)) sections", noun="fields")))
end

TextDisplay.name(::Type{<:LineModelWorkspace}) = "LineModelWorkspace"
function Base.summary(io::IO, workspace::LineModelWorkspace)
    print(io, "Line-model workspace, ", length(workspace.input.coordinates), " terminals")
end
function Base.show(io::IO, workspace::LineModelWorkspace)
    print(io, "LineModelWorkspace(", length(workspace.input.coordinates), " terminals)")
end
function Base.show(io::IO, ::MIME"text/plain", workspace::LineModelWorkspace)
    get(io, :compact, false) && return show(io, workspace)
    return TextDisplay.tree(io, "Line-model workspace",
        ((label="terminals  $(length(workspace.input.coordinates))", noun="fields"),
         (label="length  $(workspace.input.length) m", noun="fields")))
end

TextDisplay.name(::Type{<:LineModel}) = "LineModel"
function Base.summary(io::IO, result::LineModel)
    print(io, "Line model, ", length(result.coordinates), " phases, ", length(result.frequencies), " frequencies")
end
function Base.show(io::IO, result::LineModel)
    print(io, "LineModel(", length(result.coordinates), " phases, ", result.line_length, " m)")
end
function Base.show(io::IO, ::MIME"text/plain", result::LineModel)
    get(io, :compact, false) && return show(io, result)
    return TextDisplay.tree(io, "Line model",
        ((label="phases  $(length(result.coordinates))", noun="fields"),
         (label="length  $(result.line_length) m", noun="fields"),
         (label="frequency  $(length(result.frequencies)) samples", noun="fields")))
end
