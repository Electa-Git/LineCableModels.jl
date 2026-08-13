# Package conventions

---

## Function and method names

### Multi-dispatch resolution pattern

The codebase employs a consistent three-tier resolution pattern for handling user input processing through multi-dispatch. This standardized approach allows for predictable code organization and improved maintainability.

#### Resolution pattern naming structure

```julia
_resolve_<entity>           # Primary resolution function
├─ _parse_inputs_<entity>   # Type conversion & normalization 
└─ _do_resolve_<entity>     # Core implementation logic
```

The naming pattern has been carefully selected to reflect the purpose of each dispatch layer:

1. **`_resolve_<entity>`**: The primary entry point that coordinates the resolution process. The term "resolve" indicates that the function must determine an appropriate course of action based on input types that are not known in advance. This function validates inputs and delegates to specialized implementations.

2. **`_parse_inputs_<entity>`**: The intermediate layer responsible for converting diverse input types into standardized forms that can be processed by the implementation layer. This function normalizes inputs through type-specific conversions.

3. **`_do_resolve_<entity>`**: The implementation layer that performs the actual computations or transformations once inputs have been standardized. This function contains the core logic specific to each component type.

#### Function scope and visibility

All functions in this pattern are prefixed with an underscore (`_`) to indicate they are internal implementation details not intended for direct use by the package consumers.

### Calculation vs. computation pattern

The codebase distinguishes between calculation and computation methods through a clear naming convention:

#### Calculation methods (`calc_`)

Functions prefixed with `calc_` handle intermediate steps within a broader computational framework. These methods:

- Perform specific mathematical operations on well-defined inputs.
- Typically represent a single conceptual step in a larger process.
- Return intermediate results that will be used by higher-level functions.
- Are often associated with specific physical or mathematical formulations.

#### Computation methods (`comp_`)

Functions prefixed with `comp_` represent higher-level operations that perform multiple calculation steps to achieve a complete analysis. These methods:

- Coordinate multiple calculation steps toward a final result.
- Often work with complex objects rather than primitive types.
- Represent the primary technical capabilities of the toolbox.
- May store results in appropriate data structures for further processing.

This distinction reflects the hierarchical nature of the [`LineCableModels.jl`](@ref) package, where individual calculations support the broader computational objectives of modeling transmission lines and cables.

### Library management pattern

The codebase implements a consistent pattern for managing libraries of models and components through standardized naming conventions. This pattern facilitates the storage, retrieval, and management of reusable objects within the framework.

#### Library operations naming structure

```julia
store_<library>!       # Add or update an object in a library
remove_<library>!      # Remove an object from a library
save_<library>         # Writes the entire library contents to a file
list_<library>         # Display contents of a library
```

Each library operation follows a predictable naming convention:

1. **`store_<library>!`**: Adds or updates objects in the specified library. The exclamation mark indicates that this operation modifies the library state.

2. **`remove_<library>!`**: Removes an object from the specified library. The exclamation mark indicates that this operation modifies the library state.

3. **`save_<library>`**: Persists the current state of the library to external storage, typically a file. This operation does not modify the library itself.

4. **`list_<library>`**: Displays the contents of the library for inspection without modifying its state.

### Object modification pattern

For operations that modify components within larger structures (e.g., [`AbstractConductorPart`](@ref) within a parent [`ConductorGroup`](@ref)), the codebase employs the `addto_` prefix:

```julia
addto_<entity>!        # Add or modify a subcomponent within a larger component
```

The `addto_<entity>!` pattern:

- Indicates that a subcomponent is being added to or modified within a parent component.
- Always includes an exclamation mark to denote state modification.
- Typically invokes `calc_` methods to update derived properties.

This pattern allows for hierarchical composition of components while maintaining a clear distinction from library management operations.

### DataFrame view pattern

The codebase implements a consistent pattern for generating `DataFrame` views of complex objects through a standardized naming convention:

```julia
<entity>_todf        # Convert entity to a `DataFrame` representation
```

This pattern:

- Takes an entity object as input and creates a `DataFrame` visualization.
- Uses the suffix `_todf` to clearly indicate the conversion operation.
- Produces non-mutating transformations (no exclamation mark needed).
- Facilitates analysis, visualization, and reporting of complex data structures.

The `_todf` suffix provides a concise and immediately recognizable identifier for functions that expose object data in tabular format.

---

## Module organization

### Exports and visibility

- Public API functions have no underscore prefix and are explicitly exported.
- Internal functions use leading underscore prefix (`_function_name`) and are not exported.
- Modules use `@reexport` to propagate exports from submodules to parent modules.
- Exports are placed at the top of each module for visibility.

### Module structure

- Main module (`LineCableModels`) reexports from submodules.
- Submodules (`DataModel`, `Materials`, etc.) handle their own exports.
- Maximum nesting depth is 3 levels (parent → child → grandchild).
- Documentation is maintained at all levels using `DocStringExtensions`.

### Code navigation

- Module docstrings should list key exported functions and types.
- All exported items require docstrings.
- Internal implementation details use minimal documentation.

---

## Framework design pattern

### Core architecture

The `LineCableModels.jl` package implements a consistent framework pattern across all modules, designed to provide flexibility while maintaining type stability and performance. The architecture is built around three key components:

1. **Problem definition**: Defines the physics/mathematical approach
2. **Solver**: Controls execution parameters
3. **Workspace**: Centralizes all state during computation

This pattern separates *what* is being calculated (problem definition, formulation) from *how* calculations are executed (engine used, solver) and *where* data is stored (workspace).

### The workspace pattern

At the core of the framework is the Workspace pattern, exemplified by `FEMWorkspace` in the `FEMTools` module. This pattern can be replicated across other modules (e.g., `EMTWorkspace`).

A workspace:

- Acts as a centralized state container.
- Stores all intermediate computation data.
- Manages entity tracking and lookup tables.
- Provides consistent interfaces for data access.
- Maintains configuration and results.

### Type hierarchy

Each module follows a consistent type hierarchy:

```julia
AbstractProblemFormulation
  ├── FEMFormulation :> {Darwin, Electrodynamics, ...}
  ├── AbstractFDEMFormulation :> {CPEarth, CIGRE, ...}
  ├── AbstractEHEMFormulation :> {EnforceLayer, EquivalentSigma, ...}
  ├── ...
  └── [Other specialized formulations, concrete or abstract]

AbstractWorkspace
  ├── FEMWorkspace
  ├── EMTWorkspace
  ├── ...
  └── [Other specialized workspaces]

AbstractEntityData
  └── [Domain-specific entity types]
```

This hierarchy enables both specialization and shared interfaces.

### Data flow

Data flows through the system in a consistent pattern:

1. System definition (LineCableSystem instance).
2. Problem definition (physics parameters, formulations to employ).
3. Solver configuration (execution parameters).
4. Workspace initialization (state container).
5. Execution (multi-phase processing).
6. Result extraction (from workspace).

This pattern applies regardless of the specific module or calculation type.

### Multi-phase processing

All modules implement a multi-phase execution pattern with clear separation between phases. For example, the `FEMTools` module follows this pattern:

1. **Initialization phase**: Setup workspace, load configurations.
2. **Construction phase**: Create entities based on system definition (may include specific preliminary tasks, e.g. fragments/synchronization steps FEM simulations).
3. **Processing phase**: Execute main computation loops, store raw results in workspace container.
4. **Post processing phase**: Assign properties to processed entities.
5. **Result phase**: Extract and format results.

This pattern ensures clean separation of concerns, making the code more maintainable.

### State management

State is managed exclusively through the Workspace, which contains:

1. **Configuration state**: Original system, formulation, and opts.
2. **Entity state**: Collections of typed entities.
3. **Lookup maps**: Efficient mappings between entities and properties.
4. **Processing state**: Temporary calculation state.
5. **Result state**: Final calculation outputs.

This centralized approach eliminates global state and ensures thread safety.

### Implementation example: FEMTools.jl

The FEMTools module exemplifies this pattern:

- `FEMFormulation`: Physics parameters for FEM simulation.
- `FEMSolver`: Execution parameters for meshing and solving.
- `FEMWorkspace`: Central state container for all FEM operations.
- Entity types: Typed data containers for different geometric elements.
- Multi-phase workflow: Creation → Fragmentation → Identification → Assignment → Meshing → Solving → Post-processing.

### Extension to new modules

When creating new modules, the following patterns should be followed:

1. Define problem & formulation type (physics parameters).
2. Define solver type (execution parameters).
3. Define workspace type (state container).
4. Implement entity types specific to the domain.
5. Implement multi-phase workflow with clear separation.
6. Use the workspace pattern for state management.
7. Follow the standard data flow.

This framework ensures consistency, maintainability, and performance across all modules within the package.
