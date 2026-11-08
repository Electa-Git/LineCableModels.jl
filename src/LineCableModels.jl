"""
    LineCableModels

Compute the series impedances and shunt admittances of cables per unit length. The data
model declares the parts of a cable, and a problem places cables above or in the earth. A
formulation selects the formulas, and a backend executes the computation. The Gmsh and
Measurements extensions add the finite-element backend and the propagation of input
uncertainty. The ForwardDiff and IntervalArithmetic extensions add the derivatives with
respect to the inputs and the enclosures of the Bessel functions on an interval argument.
"""
module LineCableModels

include("commons/Commons.jl")
include("materials/Materials.jl")
include("earth/Earth.jl")
include("basegeometry/BaseGeometry.jl")
include("datamodel/DataModel.jl")
include("engine/Engine.jl")

using .Commons: build, resolve, validate, compute, defaults, select, nominal,
                AbstractFormulation,
                AbstractBackend, FormulationOptions, ComputationOptions
using .Materials: Material, material
using .Earth: EarthLayer, EarthModel, homogeneous, @layers
using .BaseGeometry: Placement, area, perimeter, centroid, support, distance
using .DataModel: AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle,
                  Offset, Sector, Ellipse, Polygon, AbstractLayoutPattern, CircularPattern,
                  LinearPattern, Helix,
                  LayRatio, boundary, r_in,
                  r_ex,
                  placements, overlength, pitch, PerforatedShape, ArcPolygon,
                  AbstractCablePart, Region,
                  Stack, Terminal, Repeat, Lay, Enclosure, Assembly, PlacedRegion,
                  CableDesign,
                  CableSystem, CableBlueprint, BlueprintConductor, BlueprintDielectric,
                  flatten, @design, @terminal, core, insulation, screen, sheath, bedding,
                  jacket, stranded, wires, tape, trefoil, flat
using .Engine: ε₀, μ₀, AbstractInternalImpedance,
               AbstractInsulationImpedance, AbstractInsulationAdmittance,
               AbstractSemiconAdmittance, AbstractEarthReturn, AbstractShunt,
               AbstractResistivity, AbstractEarthProperties, admittivity,
               impedance!, potential!, to_buffer!, from_buffer, initialize_buffers,
               LineParametersProblem, CoaxialFormulation, FEM, LCMBackend,
               GmshBackend, backend, LineParameters, bundle!, eliminate!, factorize!,
               transposition!,
               SpectralIntegral, integrate, Λ, reductions,
               @formulation, @air, @earth, @mixed

export Material, material, EarthLayer, EarthModel, homogeneous, @layers
export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Offset
export AbstractLayoutPattern, CircularPattern, LinearPattern, Helix, LayRatio
export area, perimeter, centroid, support, distance, boundary, r_in, r_ex, placements,
       overlength,
       pitch
export AbstractCablePart, Region, Stack, Terminal, Repeat, Lay, Enclosure, Assembly,
       CableDesign,
       CableSystem, build
export @design, @terminal, core, insulation, screen, sheath, bedding, jacket, stranded,
       wires, tape, trefoil, flat
public Commons, Materials, Earth, DataModel, validate, resolve, Sector, Ellipse, Polygon,
       PerforatedShape, ArcPolygon
public PlacedRegion, CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten
export compute, LineParametersProblem, CoaxialFormulation, FEM, LineParameters
export @formulation, @air, @earth, @mixed
public Engine, BaseGeometry, defaults, select, AbstractFormulation, AbstractBackend,
       FormulationOptions,
       ComputationOptions
public ε₀, μ₀, LCMBackend, GmshBackend, backend
public AbstractInternalImpedance, AbstractInsulationImpedance, AbstractInsulationAdmittance,
       AbstractSemiconAdmittance, AbstractEarthReturn, AbstractShunt, AbstractResistivity,
       AbstractEarthProperties
public admittivity, impedance!, potential!, to_buffer!, from_buffer,
       initialize_buffers, bundle!, eliminate!, factorize!, transposition!,
       SpectralIntegral, integrate, nominal, Λ, reductions

end
