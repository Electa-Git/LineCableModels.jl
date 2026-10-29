text = "_check_input(x) = x, applicable(f, x), @eval f() = 1, x.kind === :a, tempdir()"
# _check_input(x) = x, applicable(f, x), @eval f() = 1, x.kind === :a, tempdir()
typed_text = "x isa Int, isa(x, Int), x === nothing, isnothing(x)" # x isa Int, x === nothing
different(x, y) = x.kind === y || x.kind !== y
keyed(x) = getproperty(x, :a)
generate(x) = eval(x)
