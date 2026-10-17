function fit(selected::Formula{:ulm}, workspace, parameters, options)
    options = merge((yc=(;), h=(;), passivity=(;), dc=(;),
        line=LineParametersFormulation()), options.data)
    data = LineFrequencyData(workspace.input.segment)
    model, dc, diagnostics = fit(selected, data, workspace.input.line, options)
    return (WidebandModel(model), measure_fitting_error(model, data), measure_passivity(model),
        dc, diagnostics)
end

function fit(::Formula{:ulm}, data::LineFrequencyData, ::Nothing, options::NamedTuple)
    model = fit_ulm(data; geom=nothing, yc=options.yc, h=options.h,
        passivity=options.passivity)
    return model, (targets=nothing, achieved=nothing),
        (step_check=false, step_err=NaN, step_err_late=NaN)
end

function fit(::Formula{:ulm}, data::LineFrequencyData, line::LineParametersProblem,
        options::NamedTuple)
    geom = LineCableModelsLine(line, options.line)
    yc, h = options.yc, options.h
    dc = (targets=nothing, achieved=nothing)
    diagnostics = (;)
    if !isempty(options.dc)
        targets, convergence = dc_targets(geom, data.len)
        yc_dc = merge(options.dc, (value=targets.admittance,))
        h_dc = merge(options.dc, (value=targets.propagation,))
        if options.dc.method == :add_pole && !haskey(options.dc, :k)
            reference = dc_reference(geom, data)
            yc_dc = merge(yc_dc, (reference=reference.admittance,))
            h_dc = merge(h_dc, (reference=reference.propagation,))
        end
        yc = merge(yc, (dc=yc_dc,))
        h = merge(h, (dc=h_dc,))
        dc = (targets=targets, achieved=nothing)
        diagnostics = (dc=(convergence=convergence,),)
    end
    model = fit_ulm(data; geom, yc, h, passivity=options.passivity)
    checked = get(options.h, :step_check, true)
    diagnostics = merge(diagnostics, (step_check=checked,
        step_err=checked ? Float64(model.info[:step_err]) : NaN,
        step_err_late=checked ? Float64(model.info[:step_err_late]) : NaN))
    if !isempty(options.dc)
        dc = (targets=dc.targets,
            achieved=(admittance=eval_yc(model, 0), propagation=eval_h(model, 0)))
        diagnostics = merge(diagnostics, (dc=merge(diagnostics.dc,
            (admittance=model.info[:dc_yc], propagation=model.info[:dc_h])),))
    end
    return model, dc, diagnostics
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
