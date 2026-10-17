function fit(selected::Formula{:ulm}, workspace, parameters, options)
    options = merge((yc=(;), h=(;), passivity=(;), dc=(;),
        line=LineParametersFormulation()), options.data)
    isempty(options.dc) || throw(ArgumentError("dc targets are not implemented yet"))
    data = LineFrequencyData(workspace.input.segment)
    model, diagnostics = fit(selected, data, workspace.input.line, options)
    return (WidebandModel(model), measure_fitting_error(model, data), measure_passivity(model),
        (targets=nothing, achieved=nothing), diagnostics)
end

function fit(::Formula{:ulm}, data::LineFrequencyData, ::Nothing, options::NamedTuple)
    model = fit_ulm(data; geom=nothing, yc=options.yc, h=options.h,
        passivity=options.passivity)
    return model, (step_check=false, step_err=NaN, step_err_late=NaN)
end

function fit(::Formula{:ulm}, data::LineFrequencyData, line::LineParametersProblem,
        options::NamedTuple)
    geom = LineCableModelsLine(line, options.line)
    model = fit_ulm(data; geom, yc=options.yc, h=options.h, passivity=options.passivity)
    checked = get(options.h, :step_check, true)
    diagnostics = (step_check=checked,
        step_err=checked ? Float64(model.info[:step_err]) : NaN,
        step_err_late=checked ? Float64(model.info[:step_err_late]) : NaN)
    return model, diagnostics
end

function fit(::Formula{:fdcm}, workspace, parameters, options)
    options = merge((yc=(;), h=(;), passivity=(;), dc=(;)), options.data)
    isempty(options.dc) || throw(ArgumentError("dc targets are not implemented yet"))
    data = LineFrequencyData(workspace.input.segment)
    model = fit_fdcm(data; geom=nothing, yc=options.yc, options.h...)
    if options.passivity !== false
        controls = options.passivity === true ? (;) : options.passivity
        model, _ = enforce_passivity(model; controls...)
    end
    return (WidebandModel(model), measure_fitting_error(model, data), measure_passivity(model),
        (targets=nothing, achieved=nothing),
        (step_check=false, step_err=NaN, step_err_late=NaN))
end
