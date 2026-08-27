import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert a material library to one row per stored material.

# Arguments

- `library`: Material library.

# Returns

- A `DataFrame` with the stored material class and constitutive properties.

"""
function DataFrame(library::MaterialsLibrary)::DataFrame
    rows = [(
                name = name,
                kind = m.kind,
                rho = m.rho,
                eps_r = m.eps_r,
                mu_r = m.mu_r,
                T0 = m.T0,
                alpha = m.alpha,
                rho_thermal = m.rho_thermal,
                theta_max = m.theta_max,
                tan_delta = m.tan_delta,
                sigma_solar = m.sigma_solar
            )
            for (name, m) in library]
    data = DataFrame(rows)
    return data
end
