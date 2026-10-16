TextDisplay.name(::Type{<:TransientAnalysisProblem}) = "TransientAnalysisProblem"
function Base.summary(io::IO, problem::TransientAnalysisProblem)
    print(io, "Transient analysis problem, ", length(problem.line.system.terminal_order), " terminals")
end
function Base.show(io::IO, problem::TransientAnalysisProblem)
    print(io, "TransientAnalysisProblem(", length(problem.line.system.terminal_order),
        " terminals, tmax=", problem.window.tmax, " s)")
end
function Base.show(io::IO, ::MIME"text/plain", problem::TransientAnalysisProblem)
    get(io, :compact, false) && return show(io, problem)
    return TextDisplay.tree(io, "Transient analysis problem",
        ((label="terminals  $(length(problem.line.system.terminal_order))", noun="fields"),
         (label="window  $(problem.window.tmax) s", noun="fields"),
         (label="step  $(problem.window.step) s", noun="fields")))
end

TextDisplay.name(::Type{<:TransientFormulation}) = "TransientFormulation"
function Base.summary(io::IO, formulation::TransientFormulation)
    print(io, "Transient formulation, :", formula_id(formulation.formula))
end
function Base.show(io::IO, formulation::TransientFormulation)
    print(io, "TransientFormulation(:", formula_id(formulation.formula), ")")
end
function Base.show(io::IO, ::MIME"text/plain", formulation::TransientFormulation)
    get(io, :compact, false) && return show(io, formulation)
    return TextDisplay.tree(io, "Transient formulation",
        ((label="formula  :$(formula_id(formulation.formula))", noun="fields"),
         (label="controls  $(length(formulation_options(formulation.formula).data)) sections", noun="fields")))
end

TextDisplay.name(::Type{<:TransientWorkspace}) = "TransientWorkspace"
function Base.summary(io::IO, workspace::TransientWorkspace)
    print(io, "Transient workspace, ", length(workspace.plan.coordinates), " terminals")
end
function Base.show(io::IO, workspace::TransientWorkspace)
    print(io, "TransientWorkspace(", length(workspace.plan.coordinates), " terminals)")
end
function Base.show(io::IO, ::MIME"text/plain", workspace::TransientWorkspace)
    get(io, :compact, false) && return show(io, workspace)
    return TextDisplay.tree(io, "Transient workspace",
        ((label="terminals  $(length(workspace.plan.coordinates))", noun="fields"),
         (label="window  $(workspace.plan.window.tmax) s", noun="fields")))
end

TextDisplay.name(::Type{<:TransientResponse}) = "TransientResponse"
function Base.summary(io::IO, response::TransientResponse)
    print(io, "Transient response, ", length(response.coordinates), " terminals, ", length(response.t), " times")
end
function Base.show(io::IO, response::TransientResponse)
    print(io, "TransientResponse(", length(response.coordinates), " terminals, ", length(response.t), " times)")
end
function Base.show(io::IO, ::MIME"text/plain", response::TransientResponse)
    get(io, :compact, false) && return show(io, response)
    return TextDisplay.tree(io, "Transient response",
        ((label="voltage  $(join(size(response.voltage), '×')) V", noun="fields"),
         (label="current  $(join(size(response.current), '×')) A", noun="fields"),
         (label="time  $(length(response.t)) samples", noun="fields")))
end
