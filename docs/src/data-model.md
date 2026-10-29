# Cable data model

```@meta
CurrentModule = LineCableModels
```

The data model describes the cross-section of a cable and the cables of a line. A cable is a
tree of declared parts. [`build`](@ref) resolves the parts into a [`CableDesign`](@ref): its
placed regions and the tables that map them to terminals, repeats, lays and fills. Designs
placed in one cross-section, with a conductor index for each terminal, form a
[`CableSystem`](@ref). The computations read the system through a problem, which
[Line parameters](computing.md) describes.

## Materials

A [`Material`](@ref) stores the properties of one material. Its class is the type parameter
`K`: `:conductor`, `:insulator` or `:semicon`, so that the laws and the formulations dispatch
on it. Its nine properties are, in order:

| Field | Meaning | Unit |
| --- | --- | --- |
| `rho` | electrical resistivity at the reference temperature `T0` | Ω·m |
| `eps_r` | relative permittivity | dimensionless |
| `mu_r` | relative permeability | dimensionless |
| `T0` | reference temperature | °C |
| `alpha` | temperature coefficient of the resistivity | 1/°C |
| `rho_thermal` | thermal resistivity | K·m/W |
| `theta_max` | maximum continuous operating temperature | °C |
| `tan_delta` | polarization loss tangent, apart from the conduction that `rho` gives | dimensionless |
| `sigma_solar` | solar-absorption coefficient | dimensionless |

`Material(kind, rho, eps_r, mu_r, T0, alpha, …)` takes them by position after the class, and
[`material`](@ref) takes the same values as keywords with the same defaults.

## Shapes

A region declares its cross-section with a primitive, a subtype of
[`AbstractPrimitive`](@ref), or with an offset:

- [`Circle`](@ref): a circular cross-section, such as a solid conductor or a round wire.
- [`Annulus`](@ref): an annular cross-section, such as a tubular conductor.
- [`Rectangle`](@ref): a rectangular strand before a circular pattern bends it onto its
  course. Its width lies along the course and its height across it.
- [`Offset`](@ref) declares a layer of thickness `t` around the outline of the part inside it,
  the cross-section of an insulation, a semicon screen, a sheath or a jacket. The offset resolves against that
  outline when the design is built.

The resolution produces cross-sections that no region declares: [`AnnularSector`](@ref), the
bent rectangular strand of the same area, and [`DifferenceShape`](@ref), a cross-section with
holes, such as the fill between the wires of a course. A [`Placement`](@ref) is a translation and a rotation relative to a parent frame.

## Parts

Every part is an [`AbstractCablePart`](@ref):

- [`Region`](@ref): one material in one cross-section, with a tag that names it.
- [`Stack`](@ref): parts layered outward, each resolved against the outline of the parts
  inside it.
- [`Terminal`](@ref): the conductive regions of a part assigned to one terminal, one node of
  the circuit.
- [`Repeat`](@ref): a part repeated at the placements of a pattern. With a
  [`CircularPattern`](@ref), `n` copies lie evenly on a circle about the center of the
  enclosing part, at a declared radius or at an [`Offset`](@ref) from the outline inside them.
  The pattern bends a rectangular strand into an annular sector of the same area.
- [`Lay`](@ref): the copies of a repeat laid along the cable on a [`Helix`](@ref), whose
  [`LayRatio`](@ref) is the ratio of its lay length to the mean diameter of its course.
  [`pitch`](@ref) and [`overlength`](@ref) give the lay length and the length of the helix per
  unit length of cable.
- [`Enclosure`](@ref): a fill of one material around parts, the container without the
  outlines of its outermost contents.

## Builders

The macro [`@cable`](@ref) builds a design from the parts of its block, from the inside out,
and [`@terminal`](@ref) assigns the conductive regions of its block to one terminal. The
builders return the parts of the usual layers:

- [`core`](@ref): a solid conductor of radius `r`.
- [`insulation`](@ref), [`screen`](@ref), [`sheath`](@ref), [`bedding`](@ref) and
  [`jacket`](@ref): a layer of thickness `t` on the outline of the part inside it.
- [`stranded`](@ref): the center wire and the courses of round wires or rectangular strands
  on helices, with the fill of the conductor.
- [`armor`](@ref): one course of round wires on a helix at no offset from the outline inside
  them, with the fill of their layer.

The insulated wire of the [`@cable`](@ref) example:

```@example data-model
using LineCableModels
copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
polyethylene = Material(:insulator, 1.0e14, 2.3)
wire = @cable "insulated_wire" begin
    @terminal :core begin
        core(copper; r = 4.0e-3) # [m]
    end
    insulation(polyethylene; t = 1.0e-3) # [m]
end
r_ex(wire.outer) # [m]
```

## Designs

`build(CableDesign, id, parts...)` resolves the parts into a [`CableDesign`](@ref). The
design stores its declared tree, its placed regions from the inside out, each a
[`PlacedRegion`](@ref) with its declaration and its resolved cross-section, and the tables of
terminals, sub-cables, repeats, lays, compactions and fills. Each table row stores resolved
values in the floating type `T` of the design, to which the numbers of the cross-sections and
the materials promote, and the resolved cross-sections of all regions are in `T`.

## Blueprints

A coaxial formulation reads a cable as concentric tubes. `LineCableModels.flatten(design)`
returns that view, a [`CableBlueprint`](@ref): one [`BlueprintConductor`](@ref) per terminal
and the [`BlueprintDielectric`](@ref) layers of the interval of each conductor. The tube of a
terminal spans its conductive regions and the fills between them, and its material reproduces
the direct-current resistance and the geometric mean radius of those regions. A strand on a
helix is longer than the cable by its overlength, which raises the resistivity of the tube.
The helical path of the current adds the solenoid factor to the relative permeability of the
interval, `mu_r_ins`. The equations are in the docstring of [`flatten`](@ref).

## Systems

`build(CableSystem, designs, placements; connections, id, line_length)` places each design at
its [`Placement`](@ref) in the frame of the system. The connection of a cable maps each
terminal name of its design to a conductor index:

- An index `k ≥ 1` connects the terminal to conductor `k` of the circuit.
- Terminals with one index form one conductor, a bundle.
- The index `0` grounds the terminal, which the Kron reduction then removes.

## The earth

An [`EarthModel`](@ref) is semi-infinite air above the earth, each an [`EarthLayer`](@ref)
with its resistivity, relative permittivity and relative permeability. This release admits
one semi-infinite earth layer. [`homogeneous`](@ref) builds the model of air, of relative
permittivity and permeability 1, above one homogeneous earth.

## The interface between air and earth

The surface of the earth is the line `y = 0` of the system's frame. A cable above it lies in
air and a cable below it in the earth. A cable may touch the interface and never cross it:
[`LineParametersProblem`](@ref) throws an `ArgumentError` for a cable that crosses it. The
unified earth-return formula and the finite-element model refuse a cable that touches the
interface.
