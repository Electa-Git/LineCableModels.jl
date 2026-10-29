# Cable data model

```@meta
CurrentModule = LineCableModels
```

The data model describes the cross-section of each cable and the arrangement of the cables
in a line. A cable is declared as a tree of parts, which [`build`](@ref) resolves into a
[`CableDesign`](@ref) containing the placed regions of the cable and the tables that relate
them to terminals, repeats, lays and fills. The designs placed in one cross-section, with a
conductor index assigned to each terminal, form a [`CableSystem`](@ref). The computations
receive the system through a problem, as described in [Line parameters](computing.md).

## Materials

A [`Material`](@ref) stores the properties of one material. The class of the material,
`:conductor`, `:insulator` or `:semicon`, is the type parameter `K`, so that the material
laws and the formulations dispatch on it. The properties of a material are, in order:

| Field | Meaning | Unit |
| --- | --- | --- |
| `rho` | electrical resistivity at the reference temperature `T0` | Ω·m |
| `eps_r` | relative permittivity | dimensionless |
| `mu_r` | relative permeability | dimensionless |
| `T0` | reference temperature | °C |
| `alpha` | temperature coefficient of the resistivity | 1/°C |
| `rho_thermal` | thermal resistivity | K·m/W |
| `theta_max` | maximum continuous operating temperature | °C |
| `tan_delta` | polarization loss tangent, apart from the conduction given by `rho` | dimensionless |
| `sigma_solar` | solar-absorption coefficient | dimensionless |

The constructor `Material(kind, rho, eps_r, mu_r, T0, alpha, …)` takes the properties by
position after the class, whereas [`material`](@ref) takes the same values as keywords, with
the same defaults.

## Shapes

The cross-section of a region is declared by a primitive, a subtype of
[`AbstractPrimitive`](@ref), or by an offset:

- [`Circle`](@ref): circular cross-section, such as a solid conductor or a round wire.
- [`Annulus`](@ref): annular cross-section, such as a tubular conductor.
- [`Rectangle`](@ref): rectangular strand, before it is bent onto its layer by a circular
  pattern. Its width lies along the layer and its height across it.
- [`Offset`](@ref): layer of thickness `t` around the outline of the part beneath it. An
  offset represents an insulation, a semicon screen, a sheath or a jacket, and is resolved
  against that outline when the design is built.

The resolution also produces cross-sections that are not declared by any region, namely the
[`AnnularSector`](@ref), which is the bent rectangular strand of the same area, and the
[`DifferenceShape`](@ref), which is a cross-section with holes, such as the filler between
the wires of a layer. A [`Placement`](@ref) is a translation and a rotation relative to a
parent frame.

## Cable parts

All parts are subtypes of [`AbstractCablePart`](@ref):

- [`Region`](@ref): one material in one cross-section, identified by a tag.
- [`Stack`](@ref): parts layered outward, each one resolved against the outline of the parts
  beneath it.
- [`Terminal`](@ref): conductive regions of a part assigned to one terminal, which is one
  node of the circuit.
- [`Repeat`](@ref): part repeated at the placements of a pattern. With a
  [`CircularPattern`](@ref), `n` copies are evenly distributed on a circle around the center
  of the enclosing part, either at a declared radius or at an [`Offset`](@ref) from the
  outline beneath them. The pattern bends a rectangular strand into an annular sector of the
  same area.
- [`Lay`](@ref): copies of a repeat laid along the cable on a [`Helix`](@ref), whose
  [`LayRatio`](@ref) is the ratio of the lay length to the mean diameter of the layer. The
  functions [`pitch`](@ref) and [`overlength`](@ref) return the lay length and the length of
  the helix per unit length of cable, respectively.
- [`Enclosure`](@ref): filler of one material around other parts, occupying the container
  without the outlines of its outermost contents.

## Builders

The macro [`@cable`](@ref) builds a design from the parts listed in its block, from the
inside out, and the macro [`@terminal`](@ref) assigns the conductive regions of its block to
one terminal. The builders below return the parts of the usual cable layers:

- [`core`](@ref): solid conductor of radius `r`.
- [`insulation`](@ref), [`screen`](@ref), [`sheath`](@ref), [`bedding`](@ref) and
  [`jacket`](@ref): layer of thickness `t` over the outline of the part beneath it.
- [`stranded`](@ref): central wire and layers of round wires or rectangular strands on
  helices, with the filler of the conductor.
- [`armor`](@ref): one layer of round wires on a helix, placed directly over the outline
  beneath them and completed by the filler of the layer.

Build the insulated wire of the [`@cable`](@ref) example and inspect its outer radius:

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

## Cable designs

The call `build(CableDesign, id, parts...)` resolves the parts into a [`CableDesign`](@ref).
The design stores the declared tree, the placed regions from the inside out and the tables of
terminals, sub-cables, repeats, lays, compactions and fills. Each placed region is a
[`PlacedRegion`](@ref) containing its declaration and its resolved cross-section. The
resolved values are stored in the floating-point type `T` of the design, to which the numbers
of the cross-sections and the materials are promoted.

## Equivalent coaxial representation

The coaxial formulations represent a cable as a set of concentric tubes. The function
`LineCableModels.flatten(design)` returns this representation as a
[`CableBlueprint`](@ref), composed of one [`BlueprintConductor`](@ref) per terminal and of
the [`BlueprintDielectric`](@ref) layers in the insulation interval of each conductor. The
tube of a terminal spans its conductive regions and the fillers between them, and its
material reproduces the DC resistance and the geometric mean radius of those regions. A
strand laid on a helix is longer than the cable by its overlength, which increases the
resistivity of the tube. The helical path of the current also introduces a solenoid
correction in the relative permeability of the interval, `mu_r_ins`. The equations are given
in the docstring of [`flatten`](@ref).

## Cable systems

The call `build(CableSystem, designs, placements; connections, id, line_length)` places each
design at its [`Placement`](@ref) in the frame of the system. The connections of each cable
assign a conductor index to each terminal of its design:

- An index `k ≥ 1` connects the terminal to conductor `k` of the circuit.
- Terminals with the same index are bonded and form one conductor, called a bundle.
- The index `0` grounds the terminal, which is then eliminated by the Kron reduction.

## Earth model

An [`EarthModel`](@ref) consists of a semi-infinite air layer above the earth, each
represented by an [`EarthLayer`](@ref) with its resistivity, relative permittivity and
relative permeability. The earth is represented by one semi-infinite layer. The function [`homogeneous`](@ref) builds the model of an air layer, with relative permittivity
and permeability equal to 1, above a homogeneous earth.

## Interface between air and earth

The earth surface corresponds to the line `y = 0` in the frame of the system. Cables above
this line are located in air, and cables below it in the earth. A cable may touch
the interface but cannot cross it, and [`LineParametersProblem`](@ref) throws an
`ArgumentError` for a cable that crosses it. The unified earth-return formula and the
finite-element model also reject a cable that touches the interface.
