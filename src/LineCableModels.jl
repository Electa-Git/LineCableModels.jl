module LineCableModels

include("commons/Commons.jl")
include("materials/Materials.jl")
include("earth/Earth.jl")
include("datamodel/DataModel.jl")
include("engine/Engine.jl")

using .Commons: build, resolve, validate, compute, defaults, select, AbstractFormulation,
                AbstractBackend, FormulationOptions, ComputationOptions
using .Materials: Material, material
using .Earth: EarthLayer, EarthModel, homogeneous
using .DataModel: Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle,
                  Offset, BentStrip, AbstractLayoutPattern, CircularPattern, Helix,
                  LayRatio, area, centroid, support, distance, boundary, r_in, r_ex,
                  placements, overlength, pitch, DifferenceShape, AbstractCablePart, Region,
                  Stack, Terminal, Repeat, Lay, Enclosure, PlacedRegion, CableDesign,
                  CableSystem, CableBlueprint, BlueprintConductor, BlueprintDielectric,
                  flatten, @cable, @terminal, core, insulation, screen, sheath, bedding,
                  jacket, stranded, armor
using .Engine: vacuum_permittivity, vacuum_permeability, AbstractInternalImpedance,
               AbstractInsulationImpedance, AbstractInsulationAdmittance,
               AbstractSemiconAdmittance, AbstractEarthReturn, AbstractShunt,
               AbstractResistivity, AbstractEarthProperties, admittivity, resistivity,
               properties, impedance, impedance!, potential!, prepare!, initialize_buffers,
               LineParametersProblem, CoaxialFormulation, FEM, LCMBackend,
               GmshBackend, backend, LineParameters, bundle!, eliminate!, transposition!,
               SpectralIntegral, integrate, nominal, ray!, formulas, reductions

export Material, material, EarthLayer, EarthModel, homogeneous
export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Offset
export AbstractLayoutPattern, CircularPattern, Helix, LayRatio
export area, centroid, support, distance, boundary, r_in, r_ex, placements, overlength,
       pitch
export AbstractCablePart, Region, Stack, Terminal, Repeat, Lay, Enclosure, CableDesign,
       CableSystem, build
export @cable, @terminal, core, insulation, screen, sheath, bedding, jacket, stranded, armor
public Commons, Materials, Earth, DataModel, validate, resolve, BentStrip, DifferenceShape
public PlacedRegion, CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten
export compute, LineParametersProblem, CoaxialFormulation, FEM, LineParameters
public Engine, defaults, select, AbstractFormulation, AbstractBackend, FormulationOptions,
       ComputationOptions
public vacuum_permittivity, vacuum_permeability, LCMBackend, GmshBackend, backend
public AbstractInternalImpedance, AbstractInsulationImpedance, AbstractInsulationAdmittance,
       AbstractSemiconAdmittance, AbstractEarthReturn, AbstractShunt, AbstractResistivity,
       AbstractEarthProperties
public admittivity, resistivity, properties, impedance, impedance!, potential!, prepare!,
       initialize_buffers, bundle!, eliminate!, transposition!, SpectralIntegral, integrate,
       nominal, ray!, formulas, reductions

end
