TextDisplay.name(::Type{<:ModalAnalysisProblem}) = "ModalAnalysisProblem"
Base.summary(io::IO, problem::ModalAnalysisProblem) =
    print(io,"Modal analysis problem, ",size(problem.parameters.Z,1)," modes")
Base.show(io::IO, problem::ModalAnalysisProblem) =
    print(io,"ModalAnalysisProblem(",size(problem.parameters.Z,1)," modes, ",
        length(problem.parameters.f)," frequencies)")
function Base.show(io::IO, ::MIME"text/plain", problem::ModalAnalysisProblem)
    get(io,:compact,false) && return show(io,problem)
    return TextDisplay.tree(io,"Modal analysis problem",(
        (label="source  phase-domain LineParameters",noun="fields"),
        (label="modes  $(size(problem.parameters.Z,1))",noun="fields"),
        (label="frequency  $(first(problem.parameters.f))–$(last(problem.parameters.f)) Hz ($(length(problem.parameters.f)) samples)",noun="fields"),
    ))
end

TextDisplay.name(::Type{<:ModalAnalysisFormulation}) = "ModalAnalysisFormulation"
Base.summary(io::IO, formulation::ModalAnalysisFormulation) =
    print(io,"Modal analysis formulation, :",formula_id(formulation.formula))
Base.show(io::IO, formulation::ModalAnalysisFormulation) =
    print(io,"ModalAnalysisFormulation(:",formula_id(formulation.formula),")")
function Base.show(io::IO, ::MIME"text/plain", formulation::ModalAnalysisFormulation)
    get(io,:compact,false) && return show(io,formulation)
    return TextDisplay.tree(io,"Modal analysis formulation",(
        (label="formula  :$(formula_id(formulation.formula))",noun="fields"),
        (label="controls  $(length(formulation_options(formulation.formula).data)) sections",noun="fields"),
    ))
end

TextDisplay.name(::Type{<:ModalAnalysisWorkspace}) = "ModalAnalysisWorkspace"
Base.summary(io::IO, workspace::ModalAnalysisWorkspace) =
    print(io,"Modal analysis workspace, ",size(workspace.Ti,1)," modes")
Base.show(io::IO, workspace::ModalAnalysisWorkspace) =
    print(io,"ModalAnalysisWorkspace(",join(size(workspace.Ti),'×'),")")
function Base.show(io::IO, ::MIME"text/plain", workspace::ModalAnalysisWorkspace)
    get(io,:compact,false) && return show(io,workspace)
    return TextDisplay.tree(io,"Modal analysis workspace",(
        (label="bases  $(join(size(workspace.Ti),'×'))",noun="fields"),
        (label="roots  $(join(size(workspace.roots),'×'))",noun="fields"),
    ))
end

TextDisplay.name(::Type{<:ModalOperators}) = "ModalOperators"
Base.summary(io::IO, maps::ModalOperators) =
    print(io,"Modal operators, ",join(size(maps.Tv),'×'))
Base.show(io::IO, maps::ModalOperators) =
    print(io,"ModalOperators(Tv=",join(size(maps.Tv),'×'),
        ", Ti=",join(size(maps.Ti),'×'),")")
function Base.show(io::IO, ::MIME"text/plain", maps::ModalOperators)
    get(io,:compact,false) && return show(io,maps)
    return TextDisplay.tree(io,"Modal operators",(
        (label="Tv  $(join(size(maps.Tv),'×'))",noun="fields"),
        (label="Ti  $(join(size(maps.Ti),'×'))",noun="fields"),
    ))
end

TextDisplay.name(::Type{<:PropagationParameters}) = "PropagationParameters"
Base.summary(io::IO, source::PropagationParameters) =
    print(io,"Propagation parameters over ",source.line_length," m")
Base.show(io::IO, source::PropagationParameters) =
    print(io,"PropagationParameters(length=",source.line_length,", modes=",
        size(Tv(source),1),")")
function Base.show(io::IO, ::MIME"text/plain", source::PropagationParameters)
    get(io,:compact,false) && return show(io,source)
    return TextDisplay.tree(io,"Propagation parameters · finite segment",(
        (label="length  $(source.line_length) m",noun="fields"),
        (label="modes  $(size(Tv(source),1))",noun="fields"),
        (label="frequency  $(first(frequencies(source)))–$(last(frequencies(source))) Hz ($(length(frequencies(source))) samples)",noun="fields"),
        (label="basis  :pul",noun="fields"),
    ))
end
