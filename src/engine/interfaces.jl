"Return the series-impedance values of a line-parameter result."
function Z end

"Return the shunt-admittance values of a line-parameter result."
function Y end

function X end
function G end
function B end
function series_impedance end
function shunt_admittance end
function reactance end
function conductance end
function susceptance end
function frequencies end
function nconductors end
function nfrequencies end
"Return absolute numerical errors from an owned comparison result."
function absolute_error end
"Return reference-normalised numerical errors from an owned comparison result."
function relative_error end

"Abstract tag for the physical domain represented by line-parameter matrices."
abstract type LineParamsDomain end

"""
$(TYPEDEF)

Store one earth-return matrix interaction and its resolved physical layers.

$(TYPEDFIELDS)
"""
struct EarthPair{T <: Real}
    "Destination row."
    row::Int
    "Destination column."
    column::Int
    "Conductor heights relative to the air-earth interface \\[m\\]."
    heights::Tuple{T, T}
    "Horizontal conductor separation \\[m\\]."
    separation::T
    "Physical layer indices of the source and target conductors."
    layers::Tuple{Int, Int}
end
"Tag line parameters expressed in the physical phase domain."
struct PhaseDomain <: LineParamsDomain end
"Store the coordinate system of a calculated modal transformation."
struct ModalDomain{O, F} <: LineParamsDomain
    "Frequency-dependent phase-to-modal voltage and current operators."
    operators::O
    "Formula that resolved mode order, scaling, and phase convention."
    formula::F
end

"Return a domain value restricted to selected frequency samples."
selectdomain(domain::LineParamsDomain, _) = domain

@inline domain(::Type{PhaseDomain}) = PhaseDomain
@inline domain(::Type{<:ModalDomain}) = ModalDomain

"Return the domain tag type of a value, or `nothing` when it has no domain."
@inline domain(::Type) = nothing
@inline domain(value) = domain(typeof(value))
