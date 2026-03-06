# CableBuilder Developer Guide

## Extending the Cable Modeling DSL

This document describes the **internal grammar of the CableBuilder API** and the **required steps to add new cable elements**.

The architecture separates responsibilities into four layers:

1. **Shape payloads** — geometric primitives
2. **Cable parts** — attach electrical meaning
3. **Builders** — materialization functors
4. **Specs** — lazy blueprints supporting combinatorics (`Grid`)

All cable layers follow the same extension pattern.

---

# 1. Core Concepts

## 1.1 Cable Parts

Cable parts attach **electrical role** to **geometric shapes**.

```julia
abstract type AbstractCablePart end

@inline r_ex(p::AbstractCablePart) = r_ex(p.shape)
```

Two built-in roles exist:

```julia
struct ConductorPart{L,T,S<:AbstractShape{L,T}} <: AbstractCablePart
    tag::Symbol
    shape::S
    material::Material{T}
end

struct InsulatorPart{L,T,S<:AbstractShape{L,T}} <: AbstractCablePart
    tag::Symbol
    shape::S
    material::Material{T}
end
```

Both rely on identical promotion logic:

```julia
function ConductorPart(tag::Symbol,
                       shape::AbstractShape{L,Tshape},
                       mat::Material{Tmat}) where {L,Tshape<:Real,Tmat<:Real}

    T = promote_type(Tshape,Tmat)

    s = convert(AbstractShape{L,T}, shape)
    m = convert(Material{T}, mat)

    ConductorPart{L,T,typeof(s)}(tag,s,m)
end
```

Identical pattern applies to `InsulatorPart`.

### Contract

A valid cable part requires:

* `tag::Symbol`
* `shape <: AbstractShape`
* `material::Material`

No logic is allowed inside the struct.

Promotion must occur in outer constructors.

---

# 2. Layout and Shape System

Shapes encode **pure geometry**.

Electrical semantics are handled by cable parts.

## Layout markers

```julia
abstract type AbstractLayout end

struct Concentric <: AbstractLayout end
struct SectorShaped <: AbstractLayout end
```

## Shape interface

```julia
abstract type AbstractShape{L<:AbstractLayout,T<:Real} end
```

Every shape must implement:

```julia
r_in(shape)
r_ex(shape)
```

These accessors are **mandatory**.

All geometry logic in downstream physics relies on them.

---

# 3. Shape Implementation Pattern

Every shape file contains **three components**.

### 3.1 The Vault (struct)

Strict parametric storage.

No logic.

Example:

```julia
struct SolidCore{L,T<:Real} <: AbstractShape{L,T}
    r_ex::T
end
```

---

### 3.2 The Janitor (outer constructor)

Handles promotion.

```julia
SolidCore{L}(r_ex::T) where {L,T<:Real} = SolidCore{L,T}(r_ex)
```

---

### 3.3 The Diplomat (`convert`)

Allows upgrading precision.

Required for compatibility with:

* `Measurements`
* grid sampling
* solver promotion

```julia
function Base.convert(::Type{<:AbstractShape{L,T}},
                      s::SolidCore{L}) where {L,T<:Real}

    SolidCore{L,T}(convert(T,r_ex(s)))
end
```

---

### 3.4 Accessors

```julia
@inline r_in(s::SolidCore) = zero(typeof(s.r_ex))
@inline r_ex(s::SolidCore) = s.r_ex
```

These must always exist.

Never access fields directly outside the shape file.

---

# 4. Builders

Builders are **functors that materialize cable parts**.

They do not perform combinatorics.

They receive the current stacking radius.

---

## Example: Solid core

```julia
struct SolidCoreBuilder{P,Tgeom<:Real,Tmat<:Real}
    tag::Symbol
    r_ex::Tgeom
    mat::Material{Tmat}
end
```

Constructor:

```julia
@inline function SolidCoreBuilder{P}(tag::Symbol,
                                     r_ex::Tgeom,
                                     mat::Material{Tmat}) where {P,Tgeom,Tmat}

    SolidCoreBuilder{P,Tgeom,Tmat}(tag,r_ex,mat)
end
```

Materialization:

```julia
@inline function (b::SolidCoreBuilder{P})(current_r::T) where {P,T<:Real}

    current_r != zero(T) &&
        error("Topological violation: Solid core must be at r=0.")

    P(b.tag, SolidCore{Concentric}(b.r_ex), b.mat)
end
```

---

## Example: Tubular layer

```julia
struct TubularBuilder{P,Tgeom<:Real,Tmat<:Real}
    tag::Symbol
    t::Tgeom
    mat::Material{Tmat}
end
```

Materialization:

```julia
@inline function (b::TubularBuilder{P})(current_r::T) where {P,T<:Real}

    r_ex = current_r + b.t

    P(b.tag,
      TubularShape{Concentric}(current_r,r_ex),
      b.mat)
end
```

---

# 5. Specs (Blueprint Layer)

Specs represent **lazy parameter spaces**.

They expand into builders during iteration.

All specs subtype:

```julia
abstract type AbstractSpec{Target} end
```

---

## SolidCoreSpec

```julia
struct SolidCoreSpec{P,Ttag,Tr,M<:AbstractSpec{Material}} <:
       AbstractSpec{SolidCoreBuilder{P}}

    tag::Ttag
    r_ex::Tr
    mat::M
end
```

Constructor:

```julia
SolidCoreSpec(::Type{P},
              tag::Ttag,
              r_ex::Tr,
              mat::M) where {P,Ttag,Tr,M<:AbstractSpec{Material}}
```

---

## TubularPartSpec

```julia
struct TubularPartSpec{P,Ttag,Tt,M<:AbstractSpec{Material}} <:
       AbstractSpec{TubularBuilder{P}}

    tag::Ttag
    t::Tt
    mat::M
end
```

Constructor:

```julia
TubularPartSpec(::Type{P},
                tag::Ttag,
                t::Tt,
                mat::M)
```

---

# 6. Grid System

`Grid` represents **deterministic parameter variation**.

Example:

```julia
Grid([0.02,0.025,0.03])
```

Specs store parameters as grids.

Iteration is performed through:

```julia
Iterators.product(grid_args(spec)...)
```

Materialization occurs lazily.

---

# 7. Builder Materialization

Builders are executed by the stacking engine:

```julia
build_layer(current_r, builders)
```

Each builder:

```
(current_r) -> CablePart
```

Stacking radius evolves via:

```
current_r → r_ex(part)
```

---

# 8. API Sugar

User-facing constructors exist in DSL modules:

Example:

```julia
Conductor.Solid(tag,material;r)
Conductor.Tubular(tag,material;t)
```

They construct specs:

```
SolidCoreSpec(...)
TubularPartSpec(...)
```

---

# 9. Extending the API

To add a new cable element **all steps below are mandatory**.

---

# Step 1 — Implement Shape

File: `newshape.jl`

Required components:

```
struct NewShape{L,T} <: AbstractShape{L,T}
    fields...
end
```

Janitor:

```
function NewShape{L}(args...) where {L}
    T = promote_type(typeof(arg1), typeof(arg2)...)
    return NewShape{L, T}(convert(T, arg1), convert(T, arg2)...)
end
```

Diplomat:

```
function Base.convert(::Type{<:AbstractShape{L,T}}, s::NewShape{L}) where {L, T <: Real}
    return NewShape{L, T}(convert(T, ...), convert(T, ...), ...)
end
```

Accessors:

```
r_in(s)
r_ex(s)
```

---

# Step 2 — Implement Builder

```
struct NewShapeBuilder{P,Tgeom,Tmat}
```

Constructor:

```
NewShapeBuilder{P}(args...)
```

Functor:

```
(b::NewShapeBuilder)(current_r)
```

Must return:

```
P(tag, NewShape(...), material)
```

---

# Step 3 — Implement Spec

```
struct NewShapeSpec{P,...} <: AbstractSpec{NewShapeBuilder{P}}
```

Constructor:

```
NewShapeSpec(::Type{P}, ...)
```

---

# Step 4 — Define Grid Arguments

If parameters are gridable, ensure they propagate through `grid_args`.

---

# Step 5 — Implement DSL Constructor

Example:

```julia
function Conductor.NewShape(tag::Symbol, mat; parameters...)
```

Must return:

```
NewShapeSpec(ConductorPart,...)
```

---

# 10. Stacking Engine

Stacking is implemented through tuple recursion.

```julia
build_layer(r,builders)
```

Each layer returns:

```
(part, build_layer(...))
```

Final object:

```
CableDesign{Tuple}(parts)
```

---

# 11. Required Guarantees

Every new element must satisfy:

* No abstract fields in builders
* No closures
* Pure struct functors
* Explicit promotion
* Correct `convert` implementation
* Valid `r_in` / `r_ex` accessors

Failure to implement any of these breaks stacking or solver compatibility.

---

# End of Guide
