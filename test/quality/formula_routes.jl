@testitem "Quality / equation ownership and absence of cross-author delegation" tags=[:quality] begin
    const E=LineCableModels.Engine
    const EP=LineCableModels.Earth
    const FM=LineCableModels.FormulaMethod
    owners=(E.InternalImpedance, E.InsulationImpedance, E.EarthImpedance,
        E.InsulationAdmittance, E.SemiconAdmittance, E.EarthAdmittance,
        EP.FrequencyDependent, EP.EquivalentHomogeneous, LineCableModels.Transforms,
        LineCableModels.Materials.TemperatureDependent)
    function inspect(node, identifier, violations)
        node isa Expr || return
        if node.head === :call &&
           first(node.args) in (:earth_impedance, :earth_potential_coefficient,
            :internal_impedance, :insulation_impedance, :insulation_material, :semicon_material,
            :earth_material, :equivalent_material, :modal_operators, :temperature_resistivity)
            arg=length(node.args)>1 ? node.args[2] : nothing
            if arg isa Expr && arg.head === :call && first(arg.args) === :Val &&
               last(arg.args) isa QuoteNode
                target=last(arg.args).value
                target===identifier || push!(violations, (identifier, target))
            end
        end
        foreach(child->inspect(child, identifier, violations), node.args)
    end
    for owner in owners, identifier in owner.formulas()

        selected=owner.Formula(identifier)
        routes = if owner in (E.EarthImpedance, E.EarthAdmittance)
            bindings = FM[]
            for kind in (:self, :mutual), source in 1:2, target in 1:2
                kind === :self && source != target && continue
                heights = (source == 1 ? 1.0 : -1.0, target == 1 ? 1.0 : -1.0)
                pair = E.EarthPair(1, kind === :self ? 1 : 2, heights,
                    kind === :self ? 0.0 : 1.0, (source, target);
                    radius = kind === :self ? 0.01 : nothing)
                binding = FM(selected, pair)
                signature = Tuple{
                    Val{identifier}, typeof.(binding.arguments)..., Any, Any, Any}
                which(binding.method, signature) === owner.EQUATION_FALLBACK && continue
                push!(bindings, validate(selected, pair).equation)
            end
            @test !isempty(bindings)
            bindings
        elseif owner === E.InternalImpedance
            values(selected.binding)
        elseif owner === EP.EquivalentHomogeneous
            (validate(selected, E.EarthPair(1, 1, (1.0, 1.0), 0.0, (1, 1); radius = 0.01)).equation,)
        else
            (selected.binding,)
        end
        for route in routes
            @test route isa FM
            @test typeof(route).parameters[1] === identifier
            @test parentmodule(route.method) === owner
            @test all(arg->arg isa Val, route.arguments)
        end
        family=owner===LineCableModels.Transforms ? "transforms" :
               lowercase(String(nameof(owner)))
        location=owner in (EP.FrequencyDependent, EP.EquivalentHomogeneous) ? "earth" :
                 owner === LineCableModels.Materials.TemperatureDependent ? "materials" : "engine"
        directory=owner===LineCableModels.Transforms ?
                  joinpath(pkgdir(LineCableModels), "src", "transforms", "formulas") :
                  joinpath(pkgdir(LineCableModels), "src", location, family, "formulas")
        syntax=Meta.parseall(read(joinpath(directory, lowercase(String(identifier))*".jl"), String))
        violations=Tuple{Symbol, Symbol}[]
        inspect(syntax, identifier, violations)
        @test isempty(violations)
        custom=(args...)->args
        key=owner===E.InternalImpedance ? :inner : :contribution
        for route in routes
            @eval LineCableModels.computation_options(::$(typeof(route)), ::$(typeof(custom))) = (;)
        end
        modified=owner.Formula(identifier; hooks = NamedTuple{(key,)}((custom,)))
        @test modified.hooks[key] === custom
        @test formula_id(modified) === identifier
        if owner in (E.EarthImpedance, E.EarthAdmittance)
            unsupported = owner.Formula(identifier; hooks = (unknown = custom,))
            kind, source,
            target = map(value -> typeof(value).parameters[1], first(routes).arguments)
            pair = E.EarthPair(1, kind === :self ? 1 : 2,
                (source == 1 ? 1.0 : -1.0, target == 1 ? 1.0 : -1.0),
                kind === :self ? 0.0 : 1.0, (source, target);
                radius = kind === :self ? 0.01 : nothing)
            @test_throws ArgumentError validate(unsupported, pair)
        else
            @test_throws ArgumentError owner.Formula(identifier; hooks = (unknown = custom,))
        end
    end
    # Guard must detect semantic delegation even when arguments are modified;
    # an exact-text forwarding check would miss this call.
    violations=Tuple{Symbol, Symbol}[]
    inspect(:(earth_impedance(Val(:OtherAuthor), Val(:mutual), modified_state, pair)),
        :Carson1926, violations)
    @test violations == [(:Carson1926, :OtherAuthor)]
end
