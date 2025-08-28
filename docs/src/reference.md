# API reference

This page provides a comprehensive API reference for the [`LineCableModels.jl`](@ref) package. It documents all public modules, types, functions, and constants, organized by functional area. Each section corresponds to a major module in the package, with detailed information about parameters, return values, and usage examples.

## Contents
```@contents
Pages = ["reference.md"]
Depth = 3
```

---

## Data model
```@autodocs
Modules = [LineCableModels.DataModel]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

### Base parameters (R, L, C, G)
```@autodocs
Modules = [LineCableModels.DataModel.BaseParams]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

---

## Earth properties
```@autodocs
Modules = [LineCableModels.EarthProps]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

---

## Import & export
```@autodocs
Modules = [LineCableModels.ImportExport]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

---

## Materials library
```@autodocs
Modules = [LineCableModels.Materials]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

---

## Utilities
```@autodocs
Modules = [LineCableModels.Utils]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```

---

## Private API

#### Data model
```@autodocs
Modules = [LineCableModels.DataModel]
Order = [:module, :constant, :type, :function, :macro]
Public = false
Private = true
```

#### Earth properties
```@autodocs
Modules = [LineCableModels.EarthProps]
Order = [:module, :constant, :type, :function, :macro]
Public = false
Private = true
```

#### Materials library
```@autodocs
Modules = [LineCableModels.Materials]
Order = [:module, :constant, :type, :function, :macro]
Public = false
Private = true
```

#### Utilities
```@autodocs
Modules = [LineCableModels.Utils]
Order = [:module, :constant, :type, :function, :macro]
Public = false
Private = true
```

---

## Index
```@index
Pages   = ["reference.md"]
Order   = [:module, :constant, :type, :function, :macro]
```