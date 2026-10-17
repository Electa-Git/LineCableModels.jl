function fit(::Formula{:ulm}, workspace, parameters, options)
    options = merge((yc=(;), h=(;), passivity=(;), dc=(;)), options.data)
    isempty(options.dc) || throw(ArgumentError("dc targets are not implemented yet"))
    data = LineFrequencyData(workspace.input.segment)
    model = fit_ulm(data; geom=nothing, yc=options.yc, h=options.h,
        passivity=options.passivity)
    return (WidebandModel(model), measure_fitting_error(model, data), measure_passivity(model),
        (targets=nothing, achieved=nothing), (;))
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
        (targets=nothing, achieved=nothing), (;))
end
