# Docstrings

The following docstring standards are generally adopted across the codebase.

## General Docstring principles

1. **Placement:** Docstrings must immediately precede the code entity (struct, function, module, constant) they describe.
2. **Delimiter:** Use triple double quotes (`"""Docstring content"""`) for all docstrings, *except* for individual struct field documentation.
3. **Conciseness:** Avoid redundancy. Information should be presented clearly and concisely in the appropriate section.
4. **Tone:** Use formal, precise scientific language suitable for technical documentation. Avoid contractions, colloquialisms, and ambiguous phrasing.

## Physical unit formatting

All variables corresponding to physical quantities must be annotated with their SI units and according to the following rules:

1. **Mandatory units:** ALL arguments, return values, struct fields, and constants representing **physical quantities** MUST have their SI units specified.
2. **Dimensionless quantities:** Physical quantities that are dimensionless MUST be explicitly marked as `\\[dimensionless\\]`.
3. **Non-physical quantities:** Do *not* add unit annotations to arguments, variables, or fields that do not represent physical quantities (e.g., counters, flags, indices).
4. **Standard format:** Units MUST be enclosed in double-backslash escaped square brackets: `\\[unit\\]`.
    - **Correct:** `\\[m\\]`, `\\[Hz\\]`, `\\[Ω\\]`, `\\[H/m\\]`, `\\[dimensionless\\]`
    - **Incorrect:** `[m]`, `\[m]`, `m` (as a standalone unit identifier)
5. **Exception for example comments:** Inside ` ```julia` code blocks within the `# Examples` section, use *regular* (non-escaped) square brackets for units within comments.
    - **Correct:** ````julia result = calculation(10.0) # Output in \\[m\\]````
    - **Incorrect:** ````julia result = calculation(10.0) # Output in \\\\[m\\\\]````
6. **Common units:** Use standard SI abbreviations (e.g., `m`, `s`, `kg`, `A`, `K`, `mol`, `cd`, `Hz`, `N`, `Pa`, `J`, `W`, `C`, `V`, `F`, `Ω`, `S`, `T`, `H`, `lm`, `lx`, `Bq`, `Gy`, `Sv`, `°C`). Use the Unicode middle dot `·` for multiplication where appropriate (e.g., `\\[Ω·m\\]`).

## Mathematical formulation formatting

1. **Requirement:** Mathematical formulas rendered using LaTeX are MANDATORY *only* for functions/methods whose names start with the prefix `calc_`.
2. **Location:** For `calc_` functions, the LaTeX formula MUST be placed within a ````math ...```` block inside the `# Notes` section.
3. **Forbidden:** Do NOT include ````math``` blocks or LaTeX formulations for any functions or methods *not* prefixed with `calc_`.
4. **LaTeX escaping:** Within documentation text AND inside ````math``` blocks, all LaTeX commands (like `\frac`, `\mu`) MUST have their backslashes escaped (`\\`).
    - **Correct:** `\\mu_r`, ```math \\frac{a}{b}```
    - **Incorrect:** `\mu_r`, ```math \frac{a}{b}```

## Documentation templates

The subsections below contain templates for different types of code elements.

### Structs

- **Main docstring:** Use `$(TYPEDEF)` for the signature and `$(TYPEDFIELDS)` to list the fields automatically. Provide a concise description of the struct purpose.

    ```julia
    """
    $(TYPEDEF)

    Represents a physical entity with specific properties...

    $(TYPEDFIELDS)
    """
    struct StructName
        # Field definitions follow
    end
    ```

- **Field documentation:**
  - Place *directly above* each field definition.
  - Use single-line double quotes: `"Description with units \\[unit\\] or \\[dimensionless\\] if applicable."`
  - Do NOT use `""" """` block quotes or inline comments (`#`) for documenting struct fields.

### Constructors (inside or outside structs)

- ALL constructors MUST be documented using the `@doc` macro placed immediately before the `function` keyword or the compact assignment form (`TypeName(...) = ...`). This applies even to default constructors if explicitly defined.
- **Format:** Use `$(TYPEDSIGNATURES)`. Include standard sections (`Arguments`, `Returns`, `Examples`).

    ````julia
    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`StructName`](@ref) instance.

    # Arguments

    - `arg_name`: Description including units `\\[unit\\]` if physical.

    # Returns

    - A [`StructName`](@ref) object. [Optionally add details about initialization].

    # Examples

    ```julia
    instance = $(FUNCTIONNAME)(...) # Provide meaningful example values
    ```

    """
    function StructName(...)
        # Implementation
    end
    ````

### Functions / methods

- **Format:** Start with `$(TYPEDSIGNATURES)`. Follow the section order described.
  
    ````julia
    """
    $(TYPEDSIGNATURES)

    Concise description of the function's purpose.

    # Arguments

    - `arg1`: Description, units `\\[unit\\]` if physical. Specify `Default: value` if applicable.
    - `arg2`: Description, `\\[dimensionless\\]` if physical and dimensionless.

    # Returns

    - Description of the return value, including units `\\[unit\\]` if physical. Document multiple return values individually if using tuples.

    # Notes  (OPTIONAL - MANDATORY ONLY for `calc_` functions for the formula)

    [For `calc_` functions: Explanation and formula]
    ```math
    \\LaTeX... \\escaped... \\formula...
    ```

    # Errors (OPTIONAL)

    - Describes potential errors or exceptions thrown.

    # Examples

    ```julia
    result = $(FUNCTIONNAME)(...) # Use representative values. Add expected output comment.
    # Example: result = $(FUNCTIONNAME)(0.02, 0.01, 1.0) # Expected output: ~0.0135 [m]
    ```

    # See also (OPTIONAL)

    - [`related_package_function`](@ref)
    """
    function function_name(...)
        # Implementation
    end
    ````
  
- **Section order:**
    1. Description (no heading)
    2. `# Arguments`
    3. `# Returns`
    4. `# Notes` (Only if needed; mandatory for `calc_` functions)
    5. `# Errors` (Only if needed)
    6. `# Examples`
    7. `# See also` (Only if needed)
- **Spacing:** Ensure exactly one blank line separates the description from `# Arguments` and precedes every subsequent section heading.
- **Examples:** Use the `$(FUNCTIONNAME)` macro instead of hardcoding the function name. Use meaningful, realistic input values. Include expected output or behavior in a comment, using *non-escaped* brackets for units (`[unit]`).
- **See also:** Only link to other functions *within this package* using `[`function_name`](@ref)`. Do not link to Base Julia functions or functions from external packages unless absolutely necessary for context. Only include if the linked function provides relevant context or alternatives.

### Modules

- **Format:** The first line must be the module name indented by four spaces. Use `$(IMPORTS)` and `$(EXPORTS)` literals.

    ````julia
    """
        ModuleName

    Brief description of the module purpose within the broader package (e.g., for [`Package.jl`](index.md)).

    # Overview

    - Bullet points describing key capabilities or features provided by the module.

    # Dependencies

    $(IMPORTS)

    # Exports

    $(EXPORTS)
    """
    module ModuleName
        # Contents
    end
    ````

### Constants

- **Format:** Use a single-line docstring with double quotes (`"..."`). Include a brief description, the symbol of the constant if standard (e.g., `μ₀`), its value, and its units using the `\\[unit\\]` format.

    ```julia
    "Magnetic constant (vacuum permeability), μ₀ = 4π * 1e-7 \\[H/m\\]]."
    const μ₀ = 4π * 1e-7
    ```

## Common mistakes to avoid

Double-check the docstrings to avoid these common errors:

- **Missing `@doc` for constructors:** ALL constructors require the `@doc` macro before their definition.
- **Incorrect struct field docstrings:** Use single-line `"..."` *above* the field, not block `"""..."""` quotes or inline `#` comments.
- **Incorrect section order:** Follow the specified order for function docstring sections precisely.
- **Hard-coding function names in examples:** Always use `$(FUNCTIONNAME)`.
- **Incorrect unit formatting:** Ensure `\\[unit\\]` syntax is used everywhere except comments within `Examples` blocks (`[unit]`). Double-check escaping (`\\`) for LaTeX.
- **Adding math formulas to non-`calc_` functions:** Math blocks are *only* for functions prefixed with `calc_`.
