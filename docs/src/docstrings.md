# Docstrings

Docstrings describe supported behavior, physical meaning, units, and limitations.
This release retains the established package style; it does not require a
wholesale signature or macro rewrite.

## Signatures and sections

Place a docstring immediately before its binding. Use an explicit Julia
signature when it is clearer to a user; existing DocStringExtensions signatures
may remain where they render correctly.

Use only the sections that add information:

- `# Arguments` for non-obvious inputs;
- `# Returns` for the returned type, shape, and units;
- `# Notes` for implemented physical mathematics and material assumptions;
- `# Errors` for deliberate exceptions;
- `# Examples` for concise supported syntax.

Dedicated cross-reference lists are not maintained.

## Physical quantities

State SI units for physically meaningful arguments, fields, and returns. State
when a physical value is dimensionless. Keep the package's existing escaped
unit notation, such as `\\[m\\]`, `\\[Ω·m\\]`, and
`\\[dimensionless\\]`.

## Mathematical notes

A method that implements a physically meaningful equation has a `# Notes`
section containing the equation and definitions used by the implementation.

For example:

````julia
"""
    resistance_per_length(resistivity, area)

Return conductor resistance per unit length.

# Arguments

- `resistivity`: Electrical resistivity \\[Ω·m\\].
- `area`: Conducting cross-sectional area \\[m²\\].

# Returns

- Resistance per unit length \\[Ω/m\\].

# Notes

```math
R' = \\frac{\\rho}{A}.
```
"""
resistance_per_length(resistivity, area) = resistivity / area
````

Preserve correct existing mathematical Notes. When the implementation and
documentation disagree, inspect and test the implementation before changing
either.

## Examples

Prefer `jldoctest` for self-contained public examples. When textual output is
not part of the expected output, place assertions before the doctest's `# output`
separator and leave the expected output empty. An example that needs an
external executable, repository fixture, network access, or a graphical
interaction belongs in integration documentation instead of a docstring.

Use explicit function names in new examples so copied code is readable.
Existing DocStringExtensions example macros may remain when they render valid
Julia syntax.
