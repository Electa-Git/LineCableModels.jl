touching(a, b) = isapprox(a, b)
tied(a, b) = a ≈ b
apart(a, b) = a ≉ b
qualified(a, b) = Base.isapprox(a, b)
elementwise(as, bs) = isapprox.(as, bs)
dotted(as, bs) = as .≈ bs
partial(a, bs) = all(≈(a), bs)
passed(as, bs) = map(isapprox, as, bs)
flagged(a, b) = isapprox(a, b; nans = true)
splatted(a, b; kw...) = isapprox(a, b; kw...)
chained(a, b, c) = isapprox(a, b; rtol = 1e-8) && isapprox(b, c)
nested(a, b, c) = isapprox(a, b + (b ≈ c); atol = 1e-8)
