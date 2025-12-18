module Commons

include("docstringextension.jl")
include("consts.jl")


export get_description, add!, domain

function get_description end

function add! end

"""
Return the domain tag type for objects that carry one.

Fallback returns `nothing` for domainless objects.
"""
@inline domain(::Type) = nothing
@inline domain(x) = domain(typeof(x))

end