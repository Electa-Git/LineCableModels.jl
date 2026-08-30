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
    DirectNumericalIntegration(placement)

Select PSCAD's direct numerical integration mode for an `:overhead` or
`:underground` line model.

This selector is local to the PSCAD backend because PSCAD represents the
numerical solver choice as an earth-formulation setting.
"""
struct DirectNumericalIntegration{Placement} <: EarthImpedanceFormulation
    function DirectNumericalIntegration(placement::Symbol)
        placement in (:overhead, :underground) || throw(ArgumentError(
            "PSCAD direct numerical integration placement must be :overhead or :underground",
        ))
        return new{placement}()
    end
end

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

function formulation_options(
        ::Val{PSCADFormulation},
        options::NamedTuple
)::FormulationOptions
    isempty(options) || throw(ArgumentError("PSCAD has no formulation options"))
    return (;)
end

function Formulation(
        ::Val{:pscad};
        earth_impedance = :WedepohlWilcox1973,
        earth_admittance::NativeEarthAdmittance = NativeEarthAdmittance(),
        insulation_admittance::NativeInsulationAdmittance =
        NativeInsulationAdmittance(),
        options::NamedTuple = (;)
)
    selected_earth_impedance = earth_impedance isa Symbol ?
                               EarthImpedance.Formula(earth_impedance) :
                               earth_impedance
    selected_earth_impedance isa EarthImpedanceFormulation || throw(
        ArgumentError(
        "PSCAD earth_impedance must be a literature symbol or " *
        "EarthImpedanceFormulation"
    )
    )
    methods = (
        earth_impedance = selected_earth_impedance,
        earth_admittance,
        insulation_admittance
    )
    return PSCADFormulation(
        Val(:pscad),
        methods,
        formulation_options(Val(PSCADFormulation), options)
    )
end

description(::NativeEarthAdmittance) = "PSCAD native earth admittance"
description(::NativeInsulationAdmittance) = "PSCAD native insulation admittance"
function description(::DirectNumericalIntegration{:overhead})
    "PSCAD direct numerical integration (overhead)"
end
function description(::DirectNumericalIntegration{:underground})
    "PSCAD direct numerical integration (underground)"
end

function pscad_field(formula::EarthImpedance.Formula)
    pscad_field(Val(EarthImpedance.formula_id(formula)))
end
function pscad_value(formula::EarthImpedance.Formula)
    pscad_value(Val(EarthImpedance.formula_id(formula)))
end
function pscad_readback(formula::EarthImpedance.Formula)
    pscad_readback(Val(EarthImpedance.formula_id(formula)))
end

pscad_field(::Val{:DeriSemlyen1981}) = :EarthForm2
pscad_value(::Val{:DeriSemlyen1981}) = 0
pscad_readback(::Val{:DeriSemlyen1981}) = "DERISEMLYEN"

pscad_field(::DirectNumericalIntegration{:overhead}) = :EarthForm2
pscad_value(::DirectNumericalIntegration{:overhead}) = 2
function pscad_readback(::DirectNumericalIntegration{:overhead})
    "DIRECT_NUMERICAL_INTEGRATION"
end

pscad_field(::Val{:WedepohlWilcox1973}) = :EarthForm
pscad_value(::Val{:WedepohlWilcox1973}) = 0
pscad_readback(::Val{:WedepohlWilcox1973}) = "WEDEPOHL"

pscad_field(::DirectNumericalIntegration{:underground}) = :EarthForm
pscad_value(::DirectNumericalIntegration{:underground}) = 2
function pscad_readback(::DirectNumericalIntegration{:underground})
    "DIRECT_NUMERICAL_INTEGRATION"
end

pscad_field(::Val{:Saad1996}) = :EarthForm
pscad_value(::Val{:Saad1996}) = 3
pscad_readback(::Val{:Saad1996}) = "SAAD"

pscad_field(::Val{:Ametani2009}) = :EarthForm3
pscad_value(::Val{:Ametani2009}) = 0
pscad_readback(::Val{:Ametani2009}) = "AMETANIL"

pscad_field(::Val{:Lucca1994}) = :EarthForm3
pscad_value(::Val{:Lucca1994}) = 2
pscad_readback(::Val{:Lucca1994}) = "LUCCA"

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
        options = formulation.options
    )
end
