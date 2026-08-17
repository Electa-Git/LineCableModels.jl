"Return a tabular summary of the static air and earth layers."
function DataFrames.DataFrame(model::EarthModel)
    return DataFrames.DataFrame(
        rho = getproperty.(model.layers, :rho),
        eps_r = getproperty.(model.layers, :eps_r),
        mu_r = getproperty.(model.layers, :mu_r),
        thickness = getproperty.(model.layers, :thickness)
    )
end
