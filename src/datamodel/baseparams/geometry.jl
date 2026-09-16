"""
$(TYPEDSIGNATURES)

Return `(mean_diameter, pitch_length, overlength)` for a helical radial layer:

```math
D_e=r_{in}+r_{ex},\\qquad L_p=\\lambda D_e,\\qquad
k=\\sqrt{1+(\\pi D_e/L_p)^2}.
```

`r_in` and `r_ex` are the inner and outer layer radii in meters, with
``0\\le r_{in}\\le r_{ex}``. `lay_ratio` is the nonnegative, dimensionless
ratio ``\\lambda=L_p/D_e`` used in EN 50182. All inputs must be finite.
The returned diameter and pitch are in meters; overlength ``k`` is dimensionless.
A zero pitch returns unit overlength. In particular, `lay_ratio=0` represents
a straight layer and returns a pitch of zero.
"""
function helix(r_in::Real, r_ex::Real, lay_ratio::Real)
    rin, rex, ratio = promote(float(r_in), float(r_ex), float(lay_ratio))
    all(isfinite, (rin, rex, ratio)) || throw(DomainError(
        (rin, rex, ratio), "helix geometry must be finite"
    ))
    rin >= zero(rin) || throw(DomainError(rin, "inner radius must be nonnegative"))
    rex >= rin ||
        throw(DomainError(rex, "outer radius must not be smaller than inner radius"))
    ratio >= zero(ratio) || throw(DomainError(ratio, "lay ratio must be nonnegative"))
    mean_diameter = rin + rex
    pitch_length = ratio * mean_diameter
    overlength = iszero(pitch_length) ? one(pitch_length) :
                 sqrt(one(pitch_length) +
                      (oftype(pitch_length, π) * mean_diameter / pitch_length)^2)
    return mean_diameter, pitch_length, overlength
end

"""
$(TYPEDSIGNATURES)

Return the centers of `num_wires` equally spaced circular wires. For
``N>1`` wires, the center of wire ``i=0,\\ldots,N-1`` is

```math
C+r_l(\\cos(2\\pi i/N),\\sin(2\\pi i/N)),\\qquad
r_l=r_{in}+r_w.
```

`num_wires` is the nonnegative integer count ``N``. `radius_wire` is the
positive wire radius ``r_w`` and `r_in` is the nonnegative inner radius of
the wire layer, both in meters. `C` is the layer center `(x, y)` in meters
and defaults to `(0, 0)`. Coordinates and radii must be finite.

The result is a vector of `(x, y)` tuples in meters, ordered counterclockwise
from the positive x direction. Zero wires return an empty vector; one wire
is placed at `C`.
"""
function wire_coordinates(
        num_wires::Integer,
        radius_wire::Real,
        r_in::Real;
        C = nothing
)
    num_wires >= 0 || throw(DomainError(num_wires, "number of wires must be nonnegative"))
    center = C === nothing ? (zero(radius_wire), zero(radius_wire)) : C
    rw, rin, cx, cy = promote(
        float(radius_wire), float(r_in), float(center[1]), float(center[2])
    )
    all(isfinite, (rw, rin, cx, cy)) || throw(DomainError(
        (rw, rin, cx, cy), "wire geometry must be finite"
    ))
    rw > zero(rw) || throw(DomainError(rw, "wire radius must be positive"))
    rin >= zero(rin) || throw(DomainError(rin, "inner radius must be nonnegative"))
    radius = num_wires == 1 ? zero(rin) : rin + rw
    num_wires == 0 && return Tuple{typeof(rin), typeof(rin)}[]
    step = 2 * (one(rin) * π) / num_wires
    return [(cx + radius * cos(index * step), cy + radius * sin(index * step))
            for index in 0:(num_wires - 1)]
end

function wire_coordinates(num_wires::Integer, radius_wire::Real, r_in::Real, C::Tuple)
    wire_coordinates(num_wires, radius_wire, r_in; C)
end

"""
$(TYPEDSIGNATURES)

Return the dimensionless permeability factor for a helical solenoid:

```math
k_\\mu=1+\\frac{2\\pi^2N^2(r_i^2-r_c^2)}{\\log(r_i/r_c)}.
```

`num_turns` is the nonnegative number of turns per meter ``N``.
`r_con` is the conductor radius ``r_c`` and `r_ins` is the outer insulation
radius ``r_i``, both in meters, with ``0\\le r_c<r_i``.
`num_turns=NaN` denotes an unspecified winding and returns `1` before checking
the radii.
"""
function solenoid_factor(num_turns::Real, r_con::Real, r_ins::Real)
    turns, conductor, insulator = promote(
        float(num_turns), float(r_con), float(r_ins)
    )
    isnan(turns) && return one(turns)
    turns >= zero(turns) || throw(DomainError(turns, "number of turns must be nonnegative"))
    conductor >= zero(conductor) || throw(DomainError(
        conductor, "conductor radius must be nonnegative"
    ))
    insulator > conductor || throw(DomainError(
        insulator, "insulator radius must exceed conductor radius"
    ))
    pi_value = one(turns) * π
    return one(turns) +
           2 * pi_value^2 * turns^2 *
           (insulator^2 - conductor^2) / log(insulator / conductor)
end
