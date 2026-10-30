# API reference

```@meta
CurrentModule = LineCableModels
```

The public names of the package, module by module. The formula families of the engine, the
reductions, the scan and the functions that the families extend are on the
[theory pages](theory/contents.md).

## LineCableModels

```@autodocs
Modules = [LineCableModels]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

## Commons

```@autodocs
Modules = [LineCableModels.Commons]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

## Materials

```@autodocs
Modules = [LineCableModels.Materials]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

## Earth

```@autodocs
Modules = [LineCableModels.Earth]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

## DataModel

```@autodocs
Modules = [LineCableModels.DataModel]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@docs
build(::Type{CableDesign}, id::AbstractString, parts::AbstractCablePart...)
build(::Type{CableSystem}, designs::AbstractVector{<:CableDesign}, placements::AbstractVector{<:Placement})
Base.:*(parent::Placement, child::Placement)
```

## Engine

```@autodocs
Modules = [LineCableModels.Engine]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
Filter = t -> !(t in (LineCableModels.Engine.impedance!, LineCableModels.Engine.potential!, LineCableModels.Engine.to_buffer!, LineCableModels.Engine.from_buffer, LineCableModels.Engine.Λ, LineCableModels.Engine.initialize_buffers, LineCableModels.Engine.admittivity, LineCableModels.Engine.reductions, LineCableModels.Engine.bundle!, LineCableModels.Engine.eliminate!, LineCableModels.Engine.transposition!, LineCableModels.Engine.solve!, LineCableModels.Engine.SpectralIntegral, LineCableModels.Engine.integrate))
```

```@docs
compute(problem::LineParametersProblem, formulation::AbstractFormulation)
compute(backend::AbstractBackend, problem::LineParametersProblem, formulation::AbstractFormulation)
compute(::LCMBackend, problem::LineParametersProblem{T}, formulation::CoaxialFormulation, options::ComputationOptions{LCMBackend}) where {T <: Real}
validate(problem::LineParametersProblem, ::AbstractFormulation)
validate(problem::LineParametersProblem, formulation::Union{CoaxialFormulation, FEM})
validate(earth::Tuple{Engine.EarthReturn.Formula{ID}, Engine.EarthReturn.Formula{ID}}, ::Type{CoaxialFormulation}) where {ID}
```

## CAD

```@autodocs
Modules = [LineCableModels.CAD]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

## Gmsh extension

```@autodocs
Modules = [Main.GmshExt]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@docs
defaults(::Type{GmshBackend})
compute(backend::GmshBackend, problem::LineParametersProblem, formulation::FEM)
compute(::GmshBackend, problem::LineParametersProblem, formulation::FEM, options::ComputationOptions{GmshBackend})
```

## Index

```@index
```
