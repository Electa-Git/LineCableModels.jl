"""
$(TYPEDSIGNATURES)

Calculates the coordinates of three cables arranged in a trifoil pattern.

# Arguments

- `x0`: X-coordinate of the center point \\[m\\].
- `y0`: Y-coordinate of the center point \\[m\\].
- `r_ext`: External radius of the circular layout \\[m\\].

# Returns

- A tuple containing:
  - `xa`, `ya`: Coordinates of the top cable \\[m\\].
  - `xb`, `yb`: Coordinates of the bottom-left cable \\[m\\].
  - `xc`, `yc`: Coordinates of the bottom-right cable \\[m\\].

# Examples

```julia
xa, ya, xb, yb, xc, yc = $(FUNCTIONNAME)(0.0, 0.0, 0.035)
println((xa, ya))  # Coordinates of top cable
println((xb, yb))  # Coordinates of bottom-left cable
println((xc, yc))  # Coordinates of bottom-right cable
```
"""
function trifoil_formation(x0::T, y0::T, r_ext::T) where {T<:REALSCALAR}
  @assert r_ext > 0 "External radius must be positive"

  d = r_ext / cos(deg2rad(30))
  xa = x0
  ya = y0 + d * sin(deg2rad(90))

  xb = x0 + d * cos(deg2rad(210))
  yb = y0 + d * sin(deg2rad(210))

  xc = x0 + d * cos(deg2rad(330))
  yc = y0 + d * sin(deg2rad(330))

  return xa, ya, xb, yb, xc, yc
end

"""
$(TYPEDSIGNATURES)

Calculates the coordinates of three conductors arranged in a flat (horizontal or vertical) formation.

# Arguments

- `xc`: X-coordinate of the reference point \\[m\\].
- `yc`: Y-coordinate of the reference point \\[m\\].
- `s`: Spacing between adjacent conductors \\[m\\].
- `vertical`: Boolean flag indicating whether the formation is vertical.

# Returns

- A tuple containing:
  - `xa`, `ya`: Coordinates of the first conductor \\[m\\].
  - `xb`, `yb`: Coordinates of the second conductor \\[m\\].
  - `xc`, `yc`: Coordinates of the third conductor \\[m\\].

# Examples

```julia
# Horizontal formation
xa, ya, xb, yb, xc, yc = $(FUNCTIONNAME)(0.0, 0.0, 0.1)
println((xa, ya))  # First conductor coordinates
println((xb, yb))  # Second conductor coordinates
println((xc, yc))  # Third conductor coordinates

# Vertical formation
xa, ya, xb, yb, xc, yc = $(FUNCTIONNAME)(0.0, 0.0, 0.1, vertical=true)
```
"""
function flat_formation(xc::T, yc::T, s::T; vertical=false) where {T<:REALSCALAR}
  if vertical
    # Layout is vertical; adjust only y-coordinates
    xa, ya = xc, yc
    xb, yb = xc, yc - s
    xc, yc = xc, yc - 2s
  else
    # Layout is horizontal; adjust only x-coordinates
    xa, ya = xc, yc
    xb, yb = xc + s, yc
    xc, yc = xc + 2s, yc
  end

  return xa, ya, xb, yb, xc, yc
end