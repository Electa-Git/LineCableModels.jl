module Commons

include("docstringextension.jl")
include("consts.jl")
include("retired.jl")

export get_description, add!, domain, basis,
       Z, Y, R, X, L, G, B, C,
       series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance,
       frequencies, nconductors, nfrequencies,
       LineParamsDomain, PhaseDomain, ModalDomain

function get_description end

function add! end

"""
    basis(value) -> Symbol

Return the physical storage basis of a result container. Supported line-
parameter values return `:per_length` for distributed quantities or `:total`
for quantities integrated over the modeled line length.
"""
function basis end

"""
    Z(parameters[, i, j[, k]])

Return series impedance in the units selected by [`basis`](@ref).
Index selection is defined by the concrete result container.
"""
function Z end

"""
    Y(parameters[, i, j[, k]])

Return shunt admittance in the units selected by [`basis`](@ref).
Index selection is defined by the concrete result container.
"""
function Y end
function R end
function X end
function L end
function G end
function B end
function C end
function series_impedance end
function shunt_admittance end
function resistance end
function reactance end
function inductance end
function conductance end
function susceptance end
function capacitance end
function frequencies end
function nconductors end
function nfrequencies end

abstract type LineParamsDomain end
struct PhaseDomain <: LineParamsDomain end
struct ModalDomain <: LineParamsDomain end

"""
Return the domain tag type for objects that carry one.

Fallback returns `nothing` for domainless objects.
"""
@inline domain(::Type) = nothing
@inline domain(x) = domain(typeof(x))

end
