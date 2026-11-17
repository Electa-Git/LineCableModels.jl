"""
    touching(a, b)

Whether `a` and `b` are equal within `isapprox` at a relative tolerance of `1e-8`.

```julia
@assert touching(1.0, 1.0 + 1e-12)
@assert 1.0 ≈ 1.0 + 1e-12
```
"""
touching(a, b) = isapprox(a, b; rtol = 1e-8)
tied(a, b) = isapprox(a, b, atol = 1e-8)
named(a, b; rtol, atol) = isapprox(a, b; rtol, atol)
qualified(a, b) = Base.isapprox(a, b; rtol = 0, atol = 1e-8)
elementwise(as, bs) = isapprox.(as, bs; rtol = 1e-8)
prefixed(a, b) = ≈(a, b; atol = 1e-8)
partial(a, bs) = all(isapprox(a; rtol = 1e-8), bs)
text = "isapprox(a, b), a ≈ b" # isapprox(a, b), a ≈ b
