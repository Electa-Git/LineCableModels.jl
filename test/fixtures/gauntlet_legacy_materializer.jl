# Checkpoint-era implementation retained solely to exercise native executable
# compatibility. Keep its closure order/captures: this is the old wire format.
function _materialize_case(definition::CaseDefinition, sources::NamedTuple, nominal_problem)
    any(source -> source isa Union{AbstractGrid, Gridspace}, values(sources)) ||
        return Base.invokelatest(definition.build, sources)
    names = keys(sources)
    grids = map(values(sources)) do source
        source isa Union{AbstractGrid, Gridspace} ? source : Grid((source,))
    end
    materializer = function (args...)
        Base.invokelatest(definition.build, NamedTuple{names}(args))
    end
    return Gridspace{Engine.LineParametersProblem}(materializer,grids)
end
