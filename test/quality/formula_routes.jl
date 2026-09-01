@testitem "Quality / formula catalogues / bound route ownership" tags=[:quality] begin
    using LineCableModels

    const EN=LineCableModels.Engine
    const EP=LineCableModels.EarthProps
    const FM=LineCableModels.FormulaMethod

    catalogues=(
        (
            owner = EN.InternalImpedance,
            path = ("engine", "internalimpedance", "formulas"),
            template = :internal_impedance,
            methods = (
                :internal_impedance,
                :surface_impedance_state,
                :surface_impedances
            )
        ),
        (
            owner = EN.InsulationImpedance,
            path = ("engine", "insulationimpedance", "formulas"),
            template = :insulation_impedance,
            methods = (:insulation_impedance,)
        ),
        (
            owner = EN.InsulationAdmittance,
            path = ("engine", "insulationadmittance", "formulas"),
            template = :insulation_material,
            methods = (:insulation_material,)
        ),
        (
            owner = EN.SemiconAdmittance,
            path = ("engine", "semiconadmittance", "formulas"),
            template = :semicon_material,
            methods = (:semicon_material,)
        ),
        (
            owner = EP.FD,
            path = ("earthprops", "fd", "formulas"),
            template = :earth_material,
            methods = (:earth_material,)
        ),
        (
            owner = EP.EHEM,
            path = ("earthprops", "ehem", "formulas"),
            template = :equivalent_material,
            methods = (
                :equivalent_material,
                :transverse_propagation_constant,
                :equivalent_conductivity_step
            )
        ),
        (
            owner = LineCableModels.Transforms,
            path = ("transforms", "formulas"),
            template = :modal_operators,
            methods = (
                :modal_operators,
                :modal_basis,
                :newton_workspace,
                :newton_step!,
                :levenberg_marquardt_workspace,
                :levenberg_marquardt_step!,
                :levenberg_marquardt_residual!,
                :levenberg_marquardt_jacobian!
            )
        )
    )

    binder=(identity, selector, value)->(identity, selector, value)
    bound=FM(Val(:BoundFormula), binder, Val(:interaction))
    @test bound(:runtime) == (
        Val(:BoundFormula), Val(:interaction), :runtime
    )
    @test_throws ArgumentError FM(Val(:BoundFormula), binder, :interaction)

    root=pkgdir(LineCableModels)
    common_methods=(:assumptions, :description, :routes)
    named_method=r"(?ms)^\s*(?:@inline\s+)?function\s+([A-Za-z_][A-Za-z0-9_!]*)\s*\(\s*(.*?)\)"

    for catalogue in catalogues
        owner=catalogue.owner
        identifiers=owner.formulas()
        directory=joinpath(root, "src", catalogue.path...)
        files=sort(filter(endswith(".jl"), readdir(directory)))
        @test length(files) == length(identifiers)

        for (identifier, file) in zip(identifiers, files)
            formula=owner.Formula(identifier)
            @test formula_id(formula) === identifier
            @test !isempty(description(formula))
            @test occursin(r"\(\d{4}(?:/\d{4})?\)$", description(formula))

            selected=isdefined(owner, :routes) ?
                     values(owner.routes(formula)) : (formula.route,)
            for route in selected
                @test route isa FM
                @test typeof(route).parameters[1] === identifier
                @test nameof(route.method) === catalogue.template
                @test all(argument->argument isa Val, route.arguments)
            end

            @test_throws ArgumentError owner.Formula(identifier; unknown = true)

            source=read(joinpath(directory, file), String)
            for matched in eachmatch(named_method, source)
                method_name=Symbol(matched.captures[1])
                method_name in common_methods && continue
                @test method_name in catalogue.methods
                arguments=strip(matched.captures[2])
                identity=Regex(
                    "^(?:[A-Za-z_][A-Za-z0-9_]*)?::Val\\{:$identifier\\}"
                )
                @test occursin(identity, arguments)
                @test isnothing(match(r"\d{4}", String(method_name)))
            end
        end

        first_identifier=first(identifiers)
        custom=(arguments...)->arguments
        overridden=owner === EN.InternalImpedance ?
                   owner.Formula(first_identifier; inner = custom) :
                   owner.Formula(first_identifier; route = custom)
        @test formula_id(overridden) === first_identifier
        if owner === EN.InternalImpedance
            @test owner.routes(overridden).inner === custom
        else
            @test overridden.route === custom
        end

        @test_throws ArgumentError owner.Formula(:UnregisteredFormula)
    end

    for owner in (
            EN.InsulationImpedance,
            EN.InsulationAdmittance,
            EN.SemiconAdmittance,
            EP.FD,
            EP.EHEM,
            LineCableModels.Transforms
    )
        @test !isdefined(owner, :Functor)
    end

    insulation=EN.InsulationAdmittance.Formula(:Ametani2004)
    semicon=EN.SemiconAdmittance.Formula(:Ametani2004)
    @test formula_id(insulation) === formula_id(semicon) === :Ametani2004
    @test typeof(insulation) !== typeof(semicon)
end
