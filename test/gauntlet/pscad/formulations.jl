"""
    NativeEarthAdmittance

Select PSCAD's native earth-admittance calculation for the exported line-data
model.
"""
struct NativeEarthAdmittance <: EarthAdmittanceFormulation end

"""
    NativeInsulationAdmittance

Select PSCAD's native insulation-admittance calculation for the exported
line-data model.
"""
struct NativeInsulationAdmittance <: InsulationAdmittanceFormulation end

"""
    PSCADFormulation

Store the physical methods selected for a PSCAD line-constants calculation.
"""
struct PSCADFormulation{B, M <: NamedTuple, O <: NamedTuple} <: AbstractFormulation
    backend::B
    methods::M
    options::O
end

function Base.getproperty(formulation::PSCADFormulation, name::Symbol)
    name in fieldnames(typeof(formulation)) && return getfield(formulation, name)
    methods = getfield(formulation, :methods)
    haskey(methods, name) && return getproperty(methods, name)
    return getfield(formulation, name)
end

function _pscad_options(options)
    options isa NamedTuple ||
        throw(ArgumentError("PSCAD formulation options must be a named tuple"))
    unknown = setdiff(Set(keys(options)), Set((:output_stem,)))
    isempty(unknown) || throw(ArgumentError(
        "unknown PSCAD formulation options: $(sort!(collect(unknown)))",
    ))
    stem = String(get(options, :output_stem, "gauntlet"))
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_]{0,19}$", stem) || throw(ArgumentError(
        "PSCAD output_stem must contain 1–20 ASCII letters, digits, or underscores",
    ))
    return (output_stem = stem,)
end

function Formulation(
        ::Val{:pscad};
        earth_impedance::EarthImpedance.ReferenceEarthImpedance =
        EarthImpedance.Wedepohl(),
        earth_admittance::NativeEarthAdmittance = NativeEarthAdmittance(),
        insulation_admittance::NativeInsulationAdmittance =
        NativeInsulationAdmittance(),
        options::NamedTuple = (;)
)
    methods = (; earth_impedance, earth_admittance, insulation_admittance)
    return PSCADFormulation(
        Val(:pscad),
        methods,
        _pscad_options(options)
    )
end

description(::NativeEarthAdmittance) = "PSCAD native earth admittance"
description(::NativeInsulationAdmittance) = "PSCAD native insulation admittance"

pscad_field(::EarthImpedance.Deri) = :EarthForm2
pscad_value(::EarthImpedance.Deri) = 0
pscad_readback(::EarthImpedance.Deri) = "DERISEMLYEN"

pscad_field(::EarthImpedance.DirectNumericalIntegration{:overhead}) = :EarthForm2
pscad_value(::EarthImpedance.DirectNumericalIntegration{:overhead}) = 2
function pscad_readback(::EarthImpedance.DirectNumericalIntegration{:overhead})
    "DIRECT_NUMERICAL_INTEGRATION"
end

pscad_field(::EarthImpedance.Wedepohl) = :EarthForm
pscad_value(::EarthImpedance.Wedepohl) = 0
pscad_readback(::EarthImpedance.Wedepohl) = "WEDEPOHL"

pscad_field(::EarthImpedance.DirectNumericalIntegration{:underground}) = :EarthForm
pscad_value(::EarthImpedance.DirectNumericalIntegration{:underground}) = 2
function pscad_readback(::EarthImpedance.DirectNumericalIntegration{:underground})
    "DIRECT_NUMERICAL_INTEGRATION"
end

pscad_field(::EarthImpedance.Saad) = :EarthForm
pscad_value(::EarthImpedance.Saad) = 3
pscad_readback(::EarthImpedance.Saad) = "SAAD"

pscad_field(::EarthImpedance.Ametani) = :EarthForm3
pscad_value(::EarthImpedance.Ametani) = 0
pscad_readback(::EarthImpedance.Ametani) = "AMETANIL"

pscad_field(::EarthImpedance.Lucca) = :EarthForm3
pscad_value(::EarthImpedance.Lucca) = 2
pscad_readback(::EarthImpedance.Lucca) = "LUCCA"

function formulation_record(formulation::PSCADFormulation)
    earth_impedance = formulation.earth_impedance
    earth_admittance = formulation.earth_admittance
    insulation_admittance = formulation.insulation_admittance
    return (
        type = string(parentmodule(typeof(formulation)), ".", nameof(typeof(formulation))),
        earth_impedance = (
            type = string(
                parentmodule(typeof(earth_impedance)),
                ".",
                nameof(typeof(earth_impedance))
            ),
            description = description(earth_impedance),
            pscad_field = pscad_field(earth_impedance),
            pscad_value = pscad_value(earth_impedance),
            pscad_readback = pscad_readback(earth_impedance)
        ),
        earth_admittance = (
            type = string(
                parentmodule(typeof(earth_admittance)),
                ".",
                nameof(typeof(earth_admittance))
            ),
            description = description(earth_admittance)
        ),
        insulation_admittance = (
            type = string(
                parentmodule(typeof(insulation_admittance)),
                ".",
                nameof(typeof(insulation_admittance))
            ),
            description = description(insulation_admittance)
        ),
        options = (
            output_stem = formulation.options.output_stem,
        )
    )
end
