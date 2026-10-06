# The formula families of the package, found from the registry instead of a hand-kept list.
@testmodule FormulaFamilies begin
    using LineCableModels
    const C = LineCableModels.Commons

    # Each package module that owns a `Formula` type registered with `Commons.formulas`,
    # sorted by module name.
    function families(root::Module = LineCableModels, found = Module[])
        for name in names(root; all = true)
            isdefined(root, name) || continue
            value = getfield(root, name)
            value isa Module && parentmodule(value) === root && value !== root || continue
            families(value, found)
        end
        isdefined(root, :Formula) || return found
        F = getfield(root, :Formula)
        F isa Type && parentmodule(F) === root &&
            hasmethod(C.formulas, Tuple{Type{F}}) && push!(found, root)
        return sort!(found; by = string)
    end
end
