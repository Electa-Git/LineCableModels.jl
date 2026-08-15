"""
$(TYPEDSIGNATURES)

Report that sector-shaped cable support was removed and identify the final
legacy snapshot.
"""
function SectorParams(args...; kwargs...)
    retired_fem_sector("SectorParams")
end

"""
$(TYPEDSIGNATURES)

Report that sector-shaped cable support was removed and identify the final
legacy snapshot.
"""
function Sector(args...; kwargs...)
    retired_fem_sector("Sector")
end

"""
$(TYPEDSIGNATURES)

Report that sector-shaped cable support was removed and identify the final
legacy snapshot.
"""
function SectorInsulator(args...; kwargs...)
    retired_fem_sector("SectorInsulator")
end
