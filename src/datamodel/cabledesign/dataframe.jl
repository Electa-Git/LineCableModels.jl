import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert an eager cable design to one physical or analytical table.

`format = :regions` returns one row per resolved Region. `:terminals` returns
the current radial analytical equivalence and therefore rejects designs that
the analytical adapter does not support. `:baseparams` returns one row with
the calculated cable constants.
"""
function DataFrame(design::CableDesign, format::Symbol = :terminals)::DataFrame
    if format in (:terminals, :components)
        design.effective === nothing && throw(ArgumentError(
            "terminal equivalence requires the current radial analytical profile"
        ))
        return DataFrame(
            terminal = Symbol[item.name for item in design.effective],
            conductor_r_in = [item.conductor.r_in for item in design.effective],
            conductor_r_ex = [item.conductor.r_ex for item in design.effective],
            conductor_area = [item.conductor.cross_section for item in design.effective],
            conductor_rho = [item.conductor.material.rho for item in design.effective],
            conductor_mu_r = [item.conductor.material.mu_r for item in design.effective],
            conductor_alpha = [item.conductor.alpha for item in design.effective],
            resistance = [item.conductor.resistance for item in design.effective],
            gmr = [item.conductor.gmr for item in design.effective],
            dielectric_r_ex = [item.dielectric.r_ex for item in design.effective],
            dielectric_rho = [item.dielectric.material.rho for item in design.effective],
            dielectric_eps_r = [item.dielectric.material.eps_r for item in design.effective],
            dielectric_mu_r = [item.dielectric.material.mu_r for item in design.effective],
            capacitance = [item.dielectric.shunt_capacitance for item in design.effective],
            conductance = [item.dielectric.shunt_conductance for item in design.effective]
        )
    elseif format in (:regions, :detailed)
        return DataFrame(
            tag = Symbol[source.region.tag for source in design.geometry.regions],
            terminal = Union{Missing, Symbol}[
                source.terminal === nothing ? missing : source.terminal
                for source in design.geometry.regions
            ],
            primitive = Symbol[nameof(typeof(source.region.primitive))
                               for source in design.geometry.regions],
            material_kind = Symbol[source.region.material.kind
                                   for source in design.geometry.regions],
            area = [area(source.shape) for source in design.geometry.regions],
            centroid_x = [first(centroid(source.shape))
                          for source in design.geometry.regions],
            centroid_y = [last(centroid(source.shape))
                          for source in design.geometry.regions],
            overlength = [source.overlength for source in design.geometry.regions]
        )
    elseif format == :baseparams
        constants = base_parameters(design)
        return DataFrame(R = [constants.R], L = [constants.L], C = [constants.C])
    end
    throw(ArgumentError(
        "unsupported cable-design table format :$format; " *
        "use :terminals, :regions, or :baseparams"
    ))
end
