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

```@autodocs
Modules = [LineCableModels]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## Commons

```@autodocs
Modules = [LineCableModels.Commons]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@autodocs
Modules = [LineCableModels.Commons]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## Materials

```@autodocs
Modules = [LineCableModels.Materials]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@autodocs
Modules = [LineCableModels.Materials]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## Earth

```@autodocs
Modules = [LineCableModels.Earth]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@autodocs
Modules = [LineCableModels.Earth]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## BaseGeometry

```@autodocs
Modules = [LineCableModels.BaseGeometry]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@autodocs
Modules = [LineCableModels.BaseGeometry]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## DataModel

```@autodocs
Modules = [LineCableModels.DataModel]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@autodocs
Modules = [LineCableModels.DataModel]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## Engine

```@autodocs
Modules = [LineCableModels.Engine]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
Filter = t -> !(t in (LineCableModels.Engine.impedance!, LineCableModels.Engine.potential!, LineCableModels.Engine.to_buffer!, LineCableModels.Engine.from_buffer, LineCableModels.Engine.Λ, LineCableModels.Engine.initialize_buffers, LineCableModels.Engine.admittivity, LineCableModels.Engine.reductions, LineCableModels.Engine.bundle!, LineCableModels.Engine.eliminate!, LineCableModels.Engine.factorize!, LineCableModels.Engine.transposition!, LineCableModels.Engine.solve!, LineCableModels.Engine.SpectralIntegral, LineCableModels.Engine.integrate))
```

```@autodocs
Modules = [LineCableModels.Engine]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## Gmsh extension

```@autodocs
Modules = [Main.GmshExt]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

```@autodocs
Modules = [Main.GmshExt]
Order = [:type, :function]
Public = false
Private = true
Filter = f -> Base.ispublic(parentmodule(f), nameof(f))
```

## Index

```@index
```
