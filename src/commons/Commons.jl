module Commons

include("docstringextension.jl")
include("consts.jl")


export get_description, add!, domain, LineParamsDomain, PhaseDomain, ModalDomain

function get_description end

function add! end

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