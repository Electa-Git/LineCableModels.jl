function measure_passivity(model::ULMModel)
    measured = @timed passivity_check(model)
    @info "Line-model passivity check" wall_seconds=measured.time
    checked = measured.value
    return (checked=true, passive=checked.passive,
        bands=[(band.f1, band.f2) for band in checked.bands],
        worst=checked.worst, f_worst=checked.fworst,
        perturbation=get(model.info, :passivity, nothing))
end
