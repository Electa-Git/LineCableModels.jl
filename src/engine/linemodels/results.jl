"""
$(TYPEDEF)

Store rational admittance in the phase domain and groups for delayed propagation.
Residues have axes phase × phase × pole.

$(TYPEDFIELDS)
"""
struct WidebandModel{T}
    "Poles for admittance \\[1/s\\], residues \\[dimensionless\\] and direct matrix \\[S\\]."
    admittance::NamedTuple{(:poles, :residues, :direct),
        Tuple{Vector{Complex{T}}, Array{Complex{T}, 3}, Matrix{T}}}
    "Propagation groups with delays \\[s\\], poles \\[1/s\\] and residues \\[dimensionless\\]."
    propagation::Vector{NamedTuple{(:delay, :poles, :residues),
        Tuple{T, Vector{Complex{T}}, Array{Complex{T}, 3}}}}
end

"""
$(TYPEDEF)

Retain a real, constant modal-to-phase transformation and real-pole sections.
The transformation maps modal quantities to phase quantities.

$(TYPEDFIELDS)
"""
struct ModalModel{T}
    "Transformation matrix with real entries, phase × mode \\[dimensionless\\]."
    transformation::Matrix{T}
    "Frequency used to select the transformation \\[Hz\\]."
    transformation_frequency::T
    "Sections for modal impedance: constants, real poles, residues and delays \\[s\\]."
    impedance::Vector{NamedTuple{(:constant, :poles, :residues, :delay),
        Tuple{T, Vector{T}, Vector{T}, T}}}
    "Sections for modal propagation: constants, real poles, residues and delays \\[s\\]."
    propagation::Vector{NamedTuple{(:constant, :poles, :residues, :delay),
        Tuple{T, Vector{T}, Vector{T}, T}}}
end

"""
$(TYPEDEF)

Retain one fitted line model and its measured fitting diagnostics.
The errors are relative fractions. The `f_max` field records frequency \\[Hz\\].

$(TYPEDFIELDS)
"""
struct LineModel{M, T, D} <: AbstractCoreResult
    "Fitted wideband or modal payload."
    model::M
    "Phase labels in physical order."
    coordinates::Vector{String}
    "Segment length \\[m\\]."
    line_length::T
    "Frequencies for the fit \\[Hz\\]."
    frequencies::Vector{T}
    "Admittance and propagation errors, each with rms, max and f_max matrices."
    fitting_error::NamedTuple{(:admittance, :propagation), Tuple{
        NamedTuple{(:rms, :max, :f_max), Tuple{Matrix{T}, Matrix{T}, Matrix{T}}},
        NamedTuple{(:rms, :max, :f_max), Tuple{Matrix{T}, Matrix{T}, Matrix{T}}}}}
    "Passivity evidence and perturbation."
    passivity::NamedTuple{(:checked, :passive, :bands, :worst, :f_worst, :perturbation)}
    "Direct-current targets and achieved values, or nothing."
    dc::NamedTuple{(:targets, :achieved)}
    "Supplemental computation record."
    details::D
end

"""Select poles \\[1/s\\] by admittance or impedance and propagation group."""
function poles end
"""Select residues \\[dimensionless\\] by admittance or impedance and propagation group."""
function residues end
"""Select propagation delays \\[s\\]."""
function delays end
"""Select relative errors of the fit and the frequencies of their maxima \\[Hz\\]."""
function fitting_error end
"""Select the retained passivity evidence."""
function passivity end

poles(model::WidebandModel) = [model.admittance.poles,
    [group.poles for group in model.propagation]...]
poles(model::ModalModel) = [[section.poles for section in model.impedance]...,
    [section.poles for section in model.propagation]...]
residues(model::WidebandModel) = [model.admittance.residues,
    [group.residues for group in model.propagation]...]
residues(model::ModalModel) = [[section.residues for section in model.impedance]...,
    [section.residues for section in model.propagation]...]
delays(model::WidebandModel) = [group.delay for group in model.propagation]
delays(model::ModalModel) = [section.delay for section in model.propagation]
poles(result::LineModel) = poles(result.model)
residues(result::LineModel) = residues(result.model)
delays(result::LineModel) = delays(result.model)
fitting_error(result::LineModel) = result.fitting_error
passivity(result::LineModel) = result.passivity
details(result::LineModel) = result.details
