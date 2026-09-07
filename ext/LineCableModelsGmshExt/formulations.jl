"Evaluate the shared explicit lossless choice for a passive FEM material region."
function LineCableModels.constitutive(
        ::LineCableModelsFEM, ::Val{:default}, selected,
        material, frequency, temperature
)
    return LineCableModels.constitutive(selected, material, frequency, temperature)
end

"Evaluate the shared Ametani material relation for a passive FEM material region."
function LineCableModels.constitutive(
        ::LineCableModelsFEM, ::Val{:Ametani2004}, selected,
        material, frequency, temperature
)
    return LineCableModels.constitutive(selected, material, frequency, temperature)
end

function LineCableModels.constitutive(
        ::LineCableModelsFEM, ::Val{ID}, selected, material, frequency, temperature
) where {ID}
    throw(ArgumentError(
        "FEM constitutive adaptation for :$ID is not yet implemented; " *
        "select :default for lossless dielectrics or an explicitly supported lossy formula"))
end

function formulation_record(formulation::LineCableModelsFEM)
    return (
        requested = map(formulation.definitions) do value
            value === nothing ? "nothing" :
            value isa Symbol ? string(value) :
            applicable(LineCableModels.formula_id, value) ?
            string(LineCableModels.formula_id(value)) : repr(value)
        end,
        effective = (
            internal_impedance = nothing,
            insulation_impedance = nothing,
            earth_impedance = nothing,
            insulation_admittance = LineCableModels.formula_id(
                formulation.methods.insulation_admittance),
            semicon_admittance = LineCableModels.formula_id(
                formulation.methods.semicon_admittance),
            earth_admittance = nothing,
            earth_properties = nothing,
            equivalent_earth = nothing,
            pipe_impedance = nothing
        ),
        assumptions = (
            impedance = "Fixed quasi-TEM finite-element field equations; analytical impedance kernels are not evaluated",
            earth = "One homogeneous earth half-space; analytical earth kernels and equivalent-earth reductions are not evaluated",
            insulation_admittance = LineCableModels.description(formulation.methods.insulation_admittance),
            semicon_admittance = LineCableModels.description(formulation.methods.semicon_admittance),
            semicon_domain = "Passive material region, without electrical terminal ownership",
            pipe_impedance = "Supported enclosures are included in the FEM domain; no analytical pipe correction is applied"
        )
    )
end
