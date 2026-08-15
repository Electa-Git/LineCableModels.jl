function FormulationSet(::Val{:FEM}; kwargs...)
    retired_fem_sector("FormulationSet(:FEM)")
end

module FEM

import ...Commons: retired_fem_sector
import ...Engine: FormulationSet

export Darwin, Electrodynamics, FormulationSet, MeshTransition, calc_domain_size,
       preview_results

Darwin(args...; kwargs...) = retired_fem_sector("Darwin")
Electrodynamics(args...; kwargs...) = retired_fem_sector("Electrodynamics")
MeshTransition(args...; kwargs...) = retired_fem_sector("MeshTransition")
calc_domain_size(args...; kwargs...) = retired_fem_sector("calc_domain_size")
preview_results(args...; kwargs...) = retired_fem_sector("preview_results")

end
