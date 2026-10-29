probe(f, x) = applicable(f, x) || Base.applicable(f, x)
@eval generated() = 1
Base.@eval generated_too() = 2
switch(x) = x.kind === :a || x.kind == :b || x.kind in (:c, :d) || in(x.kind, [:e])
negated(x) = x.kind !== :a && x.kind != :b && x.kind ∉ (:c, :d)
moded(mode) = mode === :fast || :slow == mode
typed(x) = x isa Int || isa(x, Float64) || Core.isa(x, Char)
subtyped(T) = T <: Real
spelled(x) = typeof(x) <: Real || typeof(x) == Int || Int !== typeof(x)
tested(xs) = any(x -> x isa Int, xs)
mapped(xs) =
    map(xs) do x
        x isa Int
    end
